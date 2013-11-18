import time
from collections import defaultdict
from string import strip
import cPickle
    
from Bio import SeqIO
from Bio.SeqFeature import CompoundLocation

from bottle import Bottle, run, get, post, request, route, response    

HOST_NAME = 'localhost' # !!!REMEMBER TO CHANGE THIS!!!
PORT_NUMBER = 9000

#class MyHandler(BaseHTTPServer.BaseHTTPRequestHandler):
#    def do_HEAD(s):
#        s.send_response(200)
#        s.send_header("Content-type", "text")
#        s.end_headers()
#        
#    def do_GET(s):
#        """Respond to a GET request."""
#        s.send_response(200)
#        s.send_header("Content-type", "text")
#        s.end_headers()
#        s.wfile.write("<html><head><title>Title goes here.</title></head>")
#        s.wfile.write("<body><p>This is a test.</p>")
#        # If someone went to "http://something.somewhere.net/foo/bar/",
#        # then s.path equals "/foo/bar/".
#        s.wfile.write("<p>You accessed path: %s</p>" % s.path)
#        s.wfile.write("</body></html>")


def genome_pos(f, scaffolds):
    start = f.location.start + f.parent.scaffold_start
    end = f.location.end + f.parent.scaffold_start

    # Test that local and genome region are actually the same seq
    if isinstance(f.location, CompoundLocation):
        for loc in f.location.parts:
            substart = loc.start + f.parent.scaffold_start
            subend = loc.end + f.parent.scaffold_start
            seq2 =  scaffolds[f.parent.scaffold][substart:subend]
            if loc.strand == -1:
                seq1 = loc.extract(f.parent.seq).reverse_complement()
            else:
                seq1 = loc.extract(f.parent.seq)
            if str(seq1) != str(seq2):
                print seq1
                print seq2
                raw_input()
    else:
        seq2 =  scaffolds[f.parent.scaffold][start:end]
        if f.location.strand == -1:
            seq1 = f.extract(f.parent.seq).reverse_complement()
        else:
            seq1 = f.extract(f.parent.seq)
    
        if str(seq1) != str(seq2):
            print seq1
            print seq2
            raw_input()
    return start, end

def read_and_index_genome():
    gbrecords = {}
    features = defaultdict(list)
    scaffolds = {}
    #Reads genome 
    for sca in SeqIO.parse(open("C_thermophilum.scaffolds.fa"), "fasta"):
        if sca.name in scaffolds:
            raise ValueError("Duplicated scaffold %s" %sca.name)
        scaffolds[sca.name] = sca.seq
    for sca in SeqIO.parse(open("C_thermophilum.mitochondrial.fa"), "fasta"):
        if sca.name in scaffolds:
            raise ValueError("Duplicated scaffold %s" %sca.name)
        scaffolds[sca.name] = sca.seq
    for sca in SeqIO.parse(open("C_thermophilum.rrn.fa"), "fasta"):
        if sca.name in scaffolds:
            raise ValueError("Duplicated scaffold %s" %sca.name)
        scaffolds[sca.name] = sca.seq

    #Map genebank regions and features into the genome assembly
    for e in SeqIO.parse(open("C_thermophilum.annotation.gbf", "rU"), "genbank"):
        gbrecords[e.id] = e
        startpos = None
        
        # Find position in scaffolds
        for sca, seq in scaffolds.iteritems():
            try:
                startpos = str(seq).index(str(e.seq))
            except ValueError, x:
                pass
            else:
                target_sca = sca
                
        if startpos == None:
            print e.id
            continue
            
        e.scaffold_start = startpos
        e.scaffold = target_sca
        
        # Index features by type and global genome position
        for f in e.features:
            f.parent = e
            genome_start, genome_end =  genome_pos(f, scaffolds)
            features[f.type].append([target_sca, genome_start, genome_end, f])
   
    return gbrecords, features, scaffolds

    
@get('/sequence')
def sequence():
    region = request.GET["region"]
    seqs = []
    #for reg in region.split(","):
    #    ch, start, end = parse_region(reg)
    #    seqs.append(SCAFFOLDS[ch][start:end+1])
    ch, start, end = parse_region(region)
    seq = str(SCAFFOLDS[ch][start:end+1])

    return {"chromosome":ch,
            "start":start,
            "end":end,
            "strand":1,
            "sequence":seq}

@get('/gene')
def gene():
    response.content_type = "application/json"
    
    region = request.GET["region"]
    ch, start, end = parse_region(region)
    print FEATURES["gene"][0]
    target_g = [f for fch, fstart, fend, f in FEATURES["gene"] if fch == ch and (between(fstart, start, end) or between(fend, start, end))]
    
    target_p = [f for fch, fstart, fend, f in FEATURES["CDS"] if fch == ch and (between(fstart, start, end) or between(fend, start, end))]

    genes = {}
    #create new gene instances
    for gene in target_g:
        gstart = gene.location.start
        gend = gene.location.end
        gstrand = gene.location.strand
        gname = gene.qualifiers["locus_tag"][0]
        genedict = {"id": f.qualifiers["locus_tag"],
                    "name": f.qualifiers["locus_tag"],
                    "biotype":"",
                    "status":"",
                    "chromosome":ch,
                    "start":gstart,
                    "end":gend,
                    "strand":gstrand,
                    "source":"CT",
                    "description":"",
                    "transcripts":[]}
        genes[gname] = genedict
    # populate transcripts
    for prot in target_p:
        pstart = prot.location.start
        pend = prot.location.end
        pstrand = prot.location.strand
        gname = prot.qualifiers["locus_tag"][0]
        pname = prot.qualifiers["protein_id"][0]
        transdict = {"id":pname,
                     "name":pname,
                     "biotype":"",
                     "status":"",
                     "chromosome":ch,
                     "start":pstart,
                     "end":pend,
                     "strand": pstrand,
                     "genomicCodingStart":0,
                     "genomicCodingEnd":0,
                     "cdnaCodingStart":0,
                     "cdnaCodingEnd":0,
                     "cdsLength":0,
                     "proteinID":"",
                     "description":"",
                     "xrefs":[{"id":"","dbNameShort":"","dbName":"","description":""},],
                     "exons":[]
                     }
        genes[gname]["transcripts"].append(transdict)
        # populate exons
        exonnumber = 1
        for estart, eend, estrand, eseq in iter_exons(prot):
            exondict = {"id": "%s_exon%s" %(pname, exonnumber),
                         "chromosome":ch,
                         "start": estart,
                         "end": eend,
                         "strand": estrand,
                         "genomicCodingStart":0,
                         "genomicCodingEnd":0,
                         "cdnaCodingStart":0,
                         "cdnaCodingEnd":0,
                         "cdsStart":0,"cdsEnd":0,
                         "phase":-1,"exonNumber":exonnumber,
                         "sequence": str(eseq)}
            transdict["exons"].append(exondict)
            exonnumber += 1
            
    return str(genes.values())

def iter_exons(f):
    if isinstance(f.location, CompoundLocation):
        for loc in f.location.parts:
            substart = loc.start + f.parent.scaffold_start
            subend = loc.end + f.parent.scaffold_start
            seq = SCAFFOLDS[f.parent.scaffold][substart:subend]
            yield substart, subend, loc.strand, seq
    else:
        substart = f.location.start + f.parent.scaffold_start
        subend = f.location.end + f.parent.scaffold_start
        seq =  SCAFFOLDS[f.parent.scaffold][start:end]
        yield substart, subend, f.location.strand, seq

 
def between(pos, start, end):
    return True if pos>=start and pos<=end else False
      
    
    
def parse_region(region):
    ch, pos = map(strip, region.split(":"))
    start, end = map(int, pos.split("-"))
    return ch, start, end



    
if __name__ == '__main__':
    try:
        GBRECORDS, FEATURES, SCAFFOLDS = cPickle.load(open("genomedb.cache.pkl"))
    except Exception:
        GBRECORDS, FEATURES, SCAFFOLDS = read_and_index_genome()
        cPickle.dump([GBRECORDS, FEATURES, SCAFFOLDS], open("genomedb.cache.pkl", "wb"))
        
    print '\n'.join(FEATURES.keys())
        
    run(host='localhost', port=9000)    



    
    
    
    
#    
#    server_class = BaseHTTPServer.HTTPServer
#    httpd = server_class((HOST_NAME, PORT_NUMBER), MyHandler)
#    print time.asctime(), "Server Starts - %s:%s" % (HOST_NAME, PORT_NUMBER)
#    try:
#        httpd.serve_forever()
#    except KeyboardInterrupt:
#        pass
#    httpd.server_close()
#    print time.asctime(), "Server Stops - %s:%s" % (HOST_NAME, PORT_NUMBER)
