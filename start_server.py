import sys
import os
import time
from collections import defaultdict
from string import strip
import cPickle
import argparse
import re

from Bio import SeqIO
from Bio.SeqFeature import CompoundLocation
from bottle import Bottle, run, get, post, request, route, response    

BASEPATH = os.path.split(os.path.realpath(__file__))[0]
sys.path.insert(0, os.path.join(BASEPATH, "sphinx/api/"))
try:
    import sphinxapi as sphinx
except ImportError:
    print "Sphinx API could not be imported. Searches will be disabled!"


FASTA_FILES = [os.path.join(BASEPATH, "data/v1/C_thermophilum.scaffolds.v1.fa"),
               os.path.join(BASEPATH, "data/v1/C_thermophilum.mitochondrial.v1.fa"),
               os.path.join(BASEPATH, "data/v1/C_thermophilum.rrn.v1.fa")]
GBF_FILE = os.path.join(BASEPATH, "data/v1/C_thermophilum.annotation.v1.gbf")
EXPR_FILE = os.path.join(BASEPATH, "data/v1/C_thermophilum.expressed.proteins.v1.txt")
DESC_FILE = os.path.join(BASEPATH, "data/v1/C_thermophilum.additional.descr.v1.txt")
DBFILE = os.path.join(BASEPATH, "genome.db.pkl")
WEBSERVICE_PORT = 9000 # needs to be externally open or bypassed with apache proxy

# Sphinx config
SPHINX_PORT = 9001
SPHINX_CONFIG_FILE = os.path.join(BASEPATH,"ct.sphinx.conf")
SPHINX_ANNOTATION_FILE = os.path.join(BASEPATH,"ct_annotations.tab")
SEARCHD = os.path.join(BASEPATH, "sphinx/src/searchd")
INDEXER = os.path.join(BASEPATH, "sphinx/src/indexer")
SEARCHER = os.path.join(BASEPATH, "sphinx/src/search")

# The following data from the gbf file is not retrieved by web-services
AVOIDED_QUALIFIERS = set(["translation"])

def system(cmd):
    print cmd
    return os.system(cmd)

# indexing functions
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
                raw_input("Inconsistency in seqs")
    else:
        seq2 =  scaffolds[f.parent.scaffold][start:end]
        if f.location.strand == -1:
            seq1 = f.extract(f.parent.seq).reverse_complement()
        else:
            seq1 = f.extract(f.parent.seq)
    
        if str(seq1) != str(seq2):
            print seq1
            print seq2
            raw_input("Inconsistency in seqs")
    return start, end

def read_and_index_extra_info():
    gene2expr = defaultdict(dict)
    for line in open(EXPR_FILE, "rU"):
        if not line.strip() or line.startswith("#"):
            continue
        try:
            #['CTHT_0043330', 'PH domain (AA 28-163)', 'SEC3', 'Saccharomyces cerevisiae',
            # 'Ivana Vonkova',
            # 'PH domain binding specificity',
            # 'now expressed with N-terminal HisSUMO-tag and C-terminal sfGFP; could be cleaved out and re-cloned using BamHI/HindIII',
            # 'Escherichia coli BL21 (DE3)',
            # 'GGATCCGACGGTTCTGTGCCAGAAACCTACATCACTCACATCCGCATCACCGAATGGCAGAACTACCCGTCTTCCCCGCCACCACCGTCTGCACGTGCTCCGCAGTACGAAAAACCGCGTGTTATCATTGTAGCTGTGCGTAAAAGCGGTCGTCTGCGCGTTCACAAATCCAAAGAAAACGCGAACGGCACCTTCAGCATTGGCAAAACTTGGTGGCTGGACGATCTGCAGAGCATCGAATCTTTCACGTCTCCGTCTGCAAACCCAAACCTGCGCGAATGGGCTCGTGATGTCGGCTTTATTGTAACCCTCGGTAAACCGTACTATTGGGAGGCTCACTCCGACAAAGAGAAAAAATTCTTCATCGCGTCCCTGATCAAAATCTTCAACCGTTACACCGGTGGTCGTACTCCAAAGCTT',
            # 'DGSVPETYITHIRITEWQNYPSSPPPPSARAPQYEKPRVIIVAVRKSGRLRVHKSKENANGTFSIGKTWWLDDLQSIESFTSPSANPNLREWARDVGFIVTLGKPYYWEAHSDKEKKFFIASLIKIFNRYTGGRTP\n']
            #        
            (locus, prj, qgene, qorg, researcher, variant,
             functag, ex_system, ntseq, aaseq) = map(strip, line.split("\t"))
        except ValueError:
            print "Skipped line:", line
            pass
        else:
            gene2expr[locus] = {
                "Project":prj,
                "Query gene":qgene,
                "Query organism":qorg,
                "Researcher":researcher,
                "Variant":variant,
                "Functional TAG":functag,
                "Expressed system": ex_system,
                "DNA seq": ntseq,
                "Protein seq": aaseq,
            }
            
    gene2eggnog = {}
    for line in open(DESC_FILE, "rU"):
        if not line.strip() or line.startswith("LOCUS"):
            continue
        try:
            #CTHT_0000070	KOG1515	28	907	V	Arylacetamide deacetylase	
            locus, eggname, x, y, z, gtype, desc = map(strip, line.split("\t"))

        except ValueError:
            print "Skipped line:", line

        else:
            gene2eggnog[locus] = {
                "EggNOG name": eggname,
                "short description": gtype,
                "Additional description": desc,
            }
    return gene2expr, gene2eggnog
             
def read_and_index_genome():
    gbrecords = {}
    features = defaultdict(list)
    scaffolds = {}
    scaffold_genes = defaultdict(set)
    #Reads genome
    for fasta_file in FASTA_FILES:
        print "Reading assembly sequences from", fasta_file
        for sca in SeqIO.parse(open(fasta_file, "rU"), "fasta"):
            if sca.name in scaffolds:
                raise ValueError("Duplicated scaffold %s" %sca.name)
            scaffolds[sca.name] = sca.seq
    print len(scaffolds), "total scaffolds"
    
    #Map genebank regions and features into the genome assembly
    print "Reading annotations from ", GBF_FILE
    for e in SeqIO.parse(open(GBF_FILE, "rU"), "genbank"):
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
            print "Skipping GBF region without a match in FASTA files", e.id
            continue
            
        e.scaffold_start = startpos
        e.scaffold = target_sca
        
        # Index features by type and global genome position
        for f in e.features:
            f.parent = e
            genome_start, genome_end =  genome_pos(f, scaffolds)
            features[f.type].append([target_sca, genome_start, genome_end, f])
            if f.type == "gene":
                scaffold_genes[target_sca].add(f)
            
    return gbrecords, features, scaffolds, scaffold_genes

# Web services
    
@get('/all')
def all():
    response.set_header("Access-Control-Allow-Origin","*")
    response.content_type = "application/json"
    chrdict = {"species":"CT","chromosomes":[]}
    for sca in SCAFFOLDS:
        chr_info = {"cytobands":[],
                    "numberGenes":0,"name":sca,"isCircular":0,"size":len(SCAFFOLDS[sca]),"end":len(SCAFFOLDS[sca]),"start":1}
        chrdict["chromosomes"].append(chr_info)

    return {"result":chrdict}

@get('/info')
def info():
    response.set_header("Access-Control-Allow-Origin","*")
    response.content_type = "application/json"
    chrdict = {"species":"CT","chromosomes":[]}
    for sca in SCAFFOLDS:
        chr_info = {"cytobands":[],
                    "numberGenes":0,"name":sca,"isCircular":0,"size":len(SCAFFOLDS[sca]),"end":len(SCAFFOLDS[sca]),"start":1}
        
        chrdict["chromosomes"].append(chr_info)
        
    return {"result":chrdict}
    
@get('/sequence')
def sequence():
    t1 = time.time()
    response.set_header("Access-Control-Allow-Origin","*")
    response.content_type = "application/json"

    by_region = []
    region = request.GET["region"]

    for reg in region.split(","):
        ch, start, end = parse_region(reg)
        by_region.append({"id":reg, "result": {"chromosome":ch,
                                               "start":start,
                                               "end":end,
                                               "strand":"+",
                                               "sequence": str(SCAFFOLDS[ch][start:end+1])}})
    print "**Returning SEQ json in", time.time() -t1
    return as_json(by_region)


@get('/gc_content')
def gc_content():
    t1 = time.time()
    response.set_header("Access-Control-Allow-Origin","*")
    response.content_type = "application/json"
    by_region = []
    region = request.GET["region"]
    histogram = request.GET.get("histogram", None)
    interval = int(request.GET.get("interval", 100))
    if histogram:
        for reg in region.split(","):
            ch, start, end = parse_region(reg)
            seq = str(SCAFFOLDS[ch][start:end+1]).upper()
            reg_intervals = []
            for i in xrange(0, len(seq), interval):
                #gc_content = sum([1 for nt in seq[i:i+interval] if nt == "G" or nt == "C"]) / float(interval)
                gc_content =  (seq[i:i+interval].count("C") + seq[i:i+interval].count("G")) / float(interval)
                reg_intervals.append({"chromosome":ch,
                                      "start":start+i,
                                      "length":interval,
                                       "end":start+i+interval-1,
                                       "strand":"+",
                                       "features_count": gc_content})
            by_region.append({"id":reg,
                              "resultType": "frequencies",
                              "result": reg_intervals})
            
    print "**Returning GC-content json in", time.time() -t1
    return as_json(by_region)
    
@get('/gene')
def gene():
    t1 = time.time()
    response.set_header("Access-Control-Allow-Origin","*")
    response.content_type = "application/json"
    
    by_region = []
    region = request.GET["region"]

    histogram = request.GET.get("histogram", None)
    interval = int(request.GET.get("interval", 1000))
    biotypes = map(strip, request.GET.get("biotypes", "CDS,tRNA,rRNA").split(","))
    excludes = set(map(strip, request.GET.get("exclude", "").split(",")))
    exclude_transcripts = "transcripts" in excludes
    if histogram:
        for reg in region.split(","):
            ch, start, end = parse_region(reg)
            total_genes = float(len(SCAFFOLD_GENES[ch]))
            all_genes = get_exons(reg, biotypes, exclude_transcripts)
            reg_intervals = []
            for i in xrange(start, end, interval):
                features = 0 
                for g in all_genes:
                    if g["start"] >= i and g["end"] <= i+interval:
                        features += 1
                    elif g["start"] > i+interval:
                        break
                reg_intervals.append({"chromosome":ch,
                                      "start":i,
                                      "end":i+interval-1,
                                      "strand":"+",
                                      "features_count": features/total_genes})

                   
            by_region.append({"id":reg, "resultType": "frequencies", "result": reg_intervals})
    else:
        for reg in region.split(","):
            by_region.append({"id":reg, "result":get_exons(reg, biotypes, exclude_transcripts)})

    print "**Returning GENE json in", time.time() -t1        
    return as_json(by_region)

@route('/search/<q>')
def search(q=None):
    if not q:
        return
    cl = sphinx.SphinxClient()
    cl.SetServer ("localhost", SPHINX_PORT)
    cl.SetMatchMode(sphinx.SPH_MATCH_ALL)
    cl.SetLimits (0, 100, 100)
    res = cl.Query("*%s*" %q, '*')
    if not res:
        print 'query failed: %s' % cl.GetLastError()
        return

    
    if cl.GetLastWarning():
        print 'WARNING: %s\n' % cl.GetLastWarning()

    print 'Query \'%s\' retrieved %d of %d matches in %s sec' % (q, res['total'], res['total_found'], res['time'])
    print 'Query stats:'

    result = {}
    result["words"] = []
    if res.has_key('words'):
        for info in res['words']:
            txt = '"%s" found %d times in %d entries' % (info['word'], info['hits'], info['docs'])
            result["words"].append(txt)
            print txt
            
    result["matches"] = []
    MATCHER = re.compile("(%s)"%q, re.IGNORECASE)
    if res.has_key('matches'):
        n = 1
        for match in res['matches']:
                attrsdump = ''
                for attr in res['attrs']:
                        attrname = attr[0]
                        attrtype = attr[1]
                        value = match['attrs'][attrname]
                        if attrtype==sphinx.SPH_ATTR_TIMESTAMP:
                                value = time.strftime ( '%Y-%m-%d %H:%M:%S', time.localtime(value) )
                        attrsdump = '%s, %s=%s' % ( attrsdump, attrname, value )

                #print '%d. doc_id=%s, weight=%d%s' % (n, match['id'], match['weight'], attrsdump, INDEX2REGION.get(int(match['id']), ""))
                #print '%d. doc_id=%s, weight=%d%s' % (n, match['id'], match['weight'], attrsdump)
                locus, biotype, region = INDEX2REGION[int(match['id'])]
                desc = ', '.join(INDEX2DESC[int(match['id'])])
                desc = re.sub(MATCHER, "<b>\\1</b>", desc)
                result["matches"].append([region, locus, biotype, desc])
                n += 1
    return as_json(result)
    
def as_json(obj):
    return {"response":obj}

def expr_data(locus_name):
    html = ""
    if locus_name in GENE2EXPR:
        html += "<div><ul>"
        for key, value in GENE2EXPR[locus_name].iteritems():
            html += "<li>%s: %s</li>" %(key, value)
        html += "</div></ul>"
    return html

def additional_description(locus_name):
    html = ""
    if locus_name in GENE2DESC:
        html += "<div><ul>"
        for key, value in GENE2DESC[locus_name].iteritems():
            html += "<li>%s: %s</li>" %(key, value)
        html += "</div></ul>"
    return html
    
def get_exons(region, biotypes, exclude_transcripts=False):
    t1 = time.time()
    ch, start, end = parse_region(region)
    genes = {}
    
    if start > len(SCAFFOLDS[ch]) or end < 0:
        return []
    
    #create new gene instances
    for fch, fstart, fend, gene in FEATURES["gene"]:
        if str(fch) == ch and (between(fstart, start, end) or between(fend, start, end)):
            gstrand = parse_strand(gene.location.strand) 
            gname = gene.qualifiers["locus_tag"][0]
            genedict = {"id": gname,
                        "name": gname,
                        "biotype": "gene",
                        "featureType": "gene",
                        "chromosome":ch,
                        "start":fstart,
                        "end":fend,
                        "strand":gstrand,
                        "source":"",
                        "description": additional_description(gname),
                        "expression": expr_data(gname),
                        "transcripts":[]}
            genes[gname] = genedict

    print "   Number of candidate genes", len(genes)
    if len(genes) == 0:
        return []
    
    # populate transcripts
    if not exclude_transcripts:
        for ftype in biotypes: 
            for fch, fstart, fend, prot in FEATURES[ftype]:
                if str(fch) == ch and (between(fstart, start, end) or between(fend, start, end)):
                    pstrand = parse_strand(gene.location.strand)
                    gname = prot.qualifiers["locus_tag"][0]
                    if gname not in genes:
                        print gname, "gname not found!!!"
                        continue
                    transdict = {
                        "id": gname,
                        "biotype": ftype,
                        "chromosome":ch,
                        "start":fstart,
                        "end":fend,
                        "strand": pstrand,
                        "exons":[]
                        }
                    for key, value in prot.qualifiers.iteritems():
                        if key in AVOIDED_QUALIFIERS: continue
                        if isinstance(value, list) and len(value) == 1:
                            transdict[key] = value[0]
                        else:
                            transdict[key] ='; '.join(map(str, value))

                    genes[gname]["transcripts"].append(transdict)
                    genes[gname]["biotype"] = ftype
                    # populate exons
                    exonnumber = 1
                    for estart, eend, estrand, eseq in iter_exons(prot):
                        estrand = parse_strand(gene.location.strand)
                        exondict = {"id": "exon_%s" %(exonnumber),
                                     "chromosome":ch,
                                     "start": int(estart),
                                     "end": int(eend),
                                     "strand": estrand,
                                    "exonNumber":exonnumber,
                                    #"sequence": str(eseq)
                        }
                        transdict["exons"].append(exondict)
                        for key, value in prot.qualifiers.iteritems():
                            if key in AVOIDED_QUALIFIERS: continue
                            if isinstance(value, list) and len(value) == 1:
                                exondict[key] = value[0]
                            else:
                                exondict[key] ='; '.join(map(str, value))
                        exonnumber += 1

        for g, gv in genes.items():
            if not gv["transcripts"]:
                del genes[g]

    result = sorted(genes.values(), lambda a,b: cmp(a["start"], b["start"]))             
    print "   Number of final genes", len(genes), "completed in ", time.time() - t1, "secs"
    return result
    
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
        seq =  SCAFFOLDS[f.parent.scaffold][substart:subend]
        yield substart, subend, f.location.strand, seq
 
def between(pos, start, end):
    pos, start, end = map(int, [pos, start, end])
    if pos>=start and pos<=end:
        return True
    return False

def parse_strand(strand):
    if int(strand) == -1:
        return "-"
    return "+"
       
def parse_region(region):
    ch, pos = map(strip, region.split(":"))
    start, end = map(int, pos.split("-"))
    return ch, start, end
    
# Command line options
    
def refresh_db():
    gbrecords, features, scaffolds, scaffold_genes = read_and_index_genome()
    gene2expr, gene2desc = read_and_index_extra_info()
    return gbrecords, features, scaffolds, scaffold_genes, gene2expr, gene2desc

def refresh_sphinx():
    index = []
    id2region = {}
    id2desc = {}
    i = 1
    for ftype in ['gene', 'CDS', 'tRNA', 'rRNA', 'misc_RNA']:
        for fch, fstart, fend, f in FEATURES[ftype]:
            gname = f.qualifiers.get("locus_tag", ["NoName"])[0]
            txt = []
            for key, value in f.qualifiers.iteritems():
                if key in AVOIDED_QUALIFIERS: continue
                if isinstance(value, list) and len(value) == 1:
                    txt.append(value[0])
                else:
                    txt.append('; '.join(map(str, value)))
            index.append([i, gname, ', '.join(txt)])
            id2region[i] = [gname, ftype, "%s:%d-%d" %(fch, fstart, fend)]
            id2desc[i] = txt
            i += 1
            
    open(SPHINX_ANNOTATION_FILE, "w").write('\n'.join(map(lambda x: "%s\t%s\t%s" %(x[0], x[1], x[2]), index))+"\n")
    
    stop_sphinx()
    print 'Updating sphinx index... '
    system('cd %s && %s --config %s annotations' %(BASEPATH, INDEXER, SPHINX_CONFIG_FILE))
   
    return id2region, id2desc
        
def start_sphinx():
    stop_sphinx()
    print 'Starting sphinx sever... '
    system('cd %s && %s --config %s' %(BASEPATH, SEARCHD, SPHINX_CONFIG_FILE))

def stop_sphinx():
    print 'Stopping sphinx sever... '
    system('cd %s && %s --config %s --stop' %(BASEPATH, SEARCHD, SPHINX_CONFIG_FILE))
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--refresh', dest='refresh', action='store_true',
                        help='upgrades all databases and start the daemon')
    parser.add_argument('--host', dest='host', type=str,
                        help='upgrades all databases and start the daemon')

    parser.add_argument('-s', dest='search', type=str,
                        help='test a query search')

    parser.add_argument('--refresh_sphinx', dest='sphinx', action="store_true",
                        help='refresh sphinx database')
   

    args = parser.parse_args()
    if args.refresh:
        GBRECORDS, FEATURES, SCAFFOLDS, SCAFFOLD_GENES, GENE2EXPR, GENE2DESC = refresh_db()
        INDEX2REGION, INDEX2DESC = refresh_sphinx()
        cPickle.dump([GBRECORDS, FEATURES, SCAFFOLDS, SCAFFOLD_GENES, GENE2EXPR, GENE2DESC, INDEX2REGION, INDEX2DESC], open(DBFILE, "wb"))
        
    else:
        try:
            GBRECORDS, FEATURES, SCAFFOLDS, SCAFFOLD_GENES, GENE2EXPR, GENE2DESC, INDEX2REGION, INDEX2DESC = cPickle.load(open(DBFILE))
        except Exception:
            print "Unable to load cached data"
            sys.exit(1)


    if args.search:
        search(args.search)
        sys.exit(0)
    
            
    print GENE2EXPR.values()[0]
    print GENE2DESC.values()[0]
    sorted_sca = sorted(SCAFFOLDS.keys(), lambda a,b: cmp(len(SCAFFOLD_GENES[a]), len(SCAFFOLD_GENES[b])), reverse=True)
    print "Serving sequences and annotations for the following scaffolds:\n", '\n'.join(map(lambda x: '"%s with %d genes"'%(x, len(SCAFFOLD_GENES[x])), sorted_sca))
    print "The following annotations were found:\n", '\n'.join(map(lambda x: "  %s %d" %(x[0], len(x[1])), FEATURES.items()))
    
    # starts 
    try:
        HOST = open("HOST").readline().strip()
    except Exception:
        HOST = "localhost"
    if args.host:
        HOST = args.host
        
    start_sphinx()
    run(host=HOST, port=WEBSERVICE_PORT)
    stop_sphinx()
   
