import sys
import os
import time
from collections import defaultdict
from string import strip
import cPickle
import argparse
import re
import logging
from StringIO import StringIO
import gzip
import json


from Bio import SeqIO
from Bio.SeqFeature import CompoundLocation
from bottle import Bottle, run, get, post, request, route, response    
from daemonize import Daemonize

BASEPATH = os.path.split(os.path.realpath(__file__))[0]
sys.path.insert(0, os.path.join(BASEPATH, "sphinx/api/"))
os.chdir(BASEPATH)

class WhitespaceRemovingFormatter(logging.Formatter):
    def format(self, record):
        record.msg = record.msg.strip()
        return super(WhitespaceRemovingFormatter, self).format(record)
    
LOG = logging.getLogger('ctbrowser')
LOG.setLevel(logging.DEBUG)

try:
    import sphinxapi as sphinx
except ImportError:
    print "Sphinx API could not be imported. Searches will be disabled!"

AN_VERSION = 2
    
if AN_VERSION == 1:
    # Blast DB will be based/upgraded based on these fasta files
    FASTA_FILES = [os.path.join(BASEPATH, "data/v1/C_thermophilum.scaffolds.v1.fa"),
                   os.path.join(BASEPATH, "data/v1/C_thermophilum.mitochondrial.v1.fa"),
                   os.path.join(BASEPATH, "data/v1/C_thermophilum.rrn.v1.fa")]

    # Annotations CDS, rRNA, tRNA and misc_RNA are extracted from the GBF file
    GBF_FILE = os.path.join(BASEPATH, "data/v1/C_thermophilum.annotation.v1.gbf")

    # Additional annotations by locus tag are loaded from here
    EXPR_FILE = os.path.join(BASEPATH, "data/v1/C_thermophilum.expressed.proteins.v1.txt")
    DESC_FILE = os.path.join(BASEPATH, "data/v1/C_thermophilum.additional.descr.v1.txt")
elif AN_VERSION == 2:
    # Blast DB will be based/upgraded based on these fasta files
    FASTA_FILES = [os.path.join(BASEPATH, "data/v2/C_thermophilum.scaffolds.v2.fa")]

    # Annotations CDS, rRNA, tRNA and misc_RNA are extracted from the GBF file
    GBF_FILE = os.path.join(BASEPATH, "data/v2/C_thermophilum.annotation.v2.gbf")

    # Additional annotations by locus tag are loaded from here
    EXPR_FILE = os.path.join(BASEPATH, "data/v2/C_thermophilum.expressed.proteins.v1.txt")
    DESC_FILE = os.path.join(BASEPATH, "data/v2/C_thermophilum.additional.descr.v1.txt")
    UNIPROT_CONVERSION = os.path.join(BASEPATH, "data/v2/uniprot2geneid.tsv")
    UNIPROT2GENE = dict([map(strip, line.split('\t')) for line in open(UNIPROT_CONVERSION) if line.strip()])
    GENE2UNIPROT = dict([v, k] for k,v in UNIPROT2GENE.iteritems())
    REPEATS_GFF_FILE = os.path.join(BASEPATH, "data/v2/Predicted_repeats_C_thermophilum.scaffolds.fa.out.gff")
    
PEP_IMGS_DIR = os.path.join(BASEPATH, "ctbrowser/pep_img/")
    
# Where all parsed information is sotored
DBFILE = os.path.join(BASEPATH, "cache/genome.db.pkl")
JS_VARS_FILE = os.path.join(BASEPATH, "ctbrowser/ct_genome_info.js") # init vars for the js browsers
BLAST_DB_PATH = os.path.join(BASEPATH, "blastDB/")
WEBSERVICE_PORT = 9000 # needs to be externally open or bypassed with apache
                       # proxy. Currently bound to ctbrowser.embl.de.
PID_FILE=os.path.join(BASEPATH, "ct_daemon.pid")

# Sphinx configuration (used to index and fast search in annotations and names)
SEARCHD = os.path.join(BASEPATH, "sphinx/src/searchd")
INDEXER = os.path.join(BASEPATH, "sphinx/src/indexer")
SEARCHER = os.path.join(BASEPATH, "sphinx/src/search")
SPHINX_PORT = 9001
# this file has all de words indexed for each locus
SPHINX_ANNOTATION_FILE = os.path.join(BASEPATH,"cache/ct_annotations.tab")
# Define how to index all the words
SPHINX_CONFIG_FILE = os.path.join(BASEPATH,"ct.sphinx.conf")

# The following data (parsed as GBF entry qualifiers) are not retrieved by
# web-services 
AVOIDED_QUALIFIERS = set(["translation"])

COMPRESS_DATA = True

def system(cmd):
    print color(cmd, "blue")
    return os.system(cmd)

def web_return(html, response, min_len=1000):
    if COMPRESS_DATA and len(html)>=min_len:
        chtmlF = StringIO()
        z = gzip.GzipFile(fileobj=chtmlF, mode='w')
        z.write(html)
        z.close()
        chtmlF.seek(0)
        html = chtmlF.read()
        response.set_header( 'Content-encoding', 'gzip')
        response.set_header( 'Content-length', len(html))
    return html


# INDEXING FUNCTIONS (CALLED ONLY WHEN DATABASE NEEDS TO BE UPGRADED)
def genome_pos(f, scaffolds):
    ''' Checks and returns the genomic coordinates of a feature'''
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
    ''' Read custom files including info about eggnog mappings and expression experiments'''
    gene2expr = defaultdict(dict)
    if EXPR_FILE:
        for line in open(EXPR_FILE, "rU"):
            if not line.strip() or line.startswith("#"):
                continue
            try:
                # Gene_name	UniProt_ID	Domain type	e-value score (SMART)	Domain borders (AA)	Domain sequence (AA)	S
                # ource of the sequence	Expression level in E.coli (1=lowest, 5=highest)	Solubility (1=poor, 3=very good)	 Interaction with PIP containing liposome
                # CTHT_0049260
                # G0SB85
                # PH
                # 1.05E-16
                # 360-472
                # NEVVKSGYLSKCGKRNPKYNRYWFRLKGVVLSYYRDPQDLYFPSGHIDLRYGISASITDKDKEGINFTIETHNRTYYFRADSAQSAKEWVKCIQRVIFRSHNDGDSVKISLPI
                # synthetic gene
                # 5
                # 3
                # Yes
                (locus, uniprot, domtype, evalue, dom_coords, peptide, source,
                 exp_level, solubility, inter_pip) = map(strip, line.split("\t"))
            except ValueError:
                print "Skipped line:", line
                pass
            else:
                gene2expr[locus] = {
                    "Gene Id":locus,
                    "Domain type":domtype,
                    "E-value (SMART)":evalue,
                    "Domain position":dom_coords,
                    "Domain sequence":peptide,
                    "source": source,
                    "Expression level in E.coli (1=lowest, 5=highest)":exp_level,
                    "Solubility (1=lowest, 3=highest)": solubility,
                    "Interaction with PIP containing liposome": inter_pip,
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
                "eggNOG": eggname,
                "eggNOG description": gtype,
                #"GO terms": desc, # A long list of GO terms
            }
    return gene2expr, gene2eggnog
             
def read_and_index_genome():
    '''
    Uses BioPython to parse the GBF and FASTA files, extract all their genes and
    features and store them in memory. Note that GBF includes many entries, for
    which gene are in local coordinates. Also, FASTA files provide the final
    genome assembly, so every GBF entrie is required to be mapped to each
    scaffold in the FASTA files, and positions converted to the global genomic
    coordinate system. '''
    
    gbrecords = {}
    features = defaultdict(list)
    scaffolds = {}
    scaffold_genes = defaultdict(set)
    # Reads genome assembly (Scaffolds) 
    for fasta_file in FASTA_FILES:
        print color("Reading assembly sequences from", "green"), fasta_file
        for sca in SeqIO.parse(open(fasta_file, "rU"), "fasta"):
            if sca.name in scaffolds:
                raise ValueError("Duplicated scaffold %s" %sca.name)
            scaffolds[sca.name] = sca.seq
    print len(scaffolds), "total scaffolds"
    
    # Map genebank regions and features into the genome assembly. Every GBF
    # entry should find a match in a scaffold.
    print color("Reading annotations from ", "green"), GBF_FILE
    for e in SeqIO.parse(open(GBF_FILE, "rU"), "genbank"):
        gbrecords[e.id] = e
        startpos = None
        
        # Find position of this region in scaffolds
        for sca, seq in scaffolds.iteritems():
            try:
                startpos = str(seq).index(str(e.seq))
            except ValueError, x:
                pass
            else:
                target_sca = sca

        # Orphan GBF region! If this occurs, manual intervention may be
        # necessary
        if startpos == None:
            print "Skipping GBF region without a match in FASTA files", e.id
            continue
            
        e.scaffold_start = startpos
        e.scaffold = target_sca
        
        # Index the global position of every feature in this GBF entry.
        for f in e.features:
            f.parent = e
            genome_start, genome_end =  genome_pos(f, scaffolds)
            features[f.type].append([target_sca, genome_start, genome_end, f])
            if f.type == "gene":
                scaffold_genes[target_sca].add(f)
            
    return gbrecords, features, scaffolds, scaffold_genes

# THESE ARE WEB SERVICES PROVIDING JSON DATA TO THE GENOME VIEWER APPLICATION
    
@get('/sequence')
def sequence():
    t1 = time.time()
    response.set_header("Access-Control-Allow-Origin","*")
    response.content_type = "application/json"

    by_region = []
    region = request.GET["region"]

    for reg in region.split(","):
        ch, start, end = parse_region(reg)
        #translate to string coordinates
        if end < 1:
            return as_json([])
        if start < 1:
            start = 1
        start -= 1
        
        seq = str(SCAFFOLDS[ch][start:end]).upper()
        print seq, start, end
        by_region.append({"id":reg, "result": {"chromosome":ch,
                                               "start":start,
                                               "end":end,
                                               "strand":"+",
                                               "sequence": seq}})
    print "**Returning SEQ json in", time.time() -t1
    return web_return(as_json(by_region), response)


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
            if end < 1:
                return as_json([])
            if start < 1:
                start = 1
            start -= 1
            
            seq = str(SCAFFOLDS[ch][start:end]).upper()

            reg_intervals = []
            for i in xrange(0, len(seq), interval):
                gc_content =  (seq[i:i+interval].count("C") + seq[i:i+interval].count("G")) / float(interval)
                reg_intervals.append({"chromosome":ch,
                                      "start":start+i+1,
                                      "length":interval,
                                       "end":start+i+interval,
                                       "strand":"+",
                                       "features_count": gc_content})
            by_region.append({"id":reg,
                              "resultType": "frequencies",
                              "result": reg_intervals})
            
    print "**Returning GC-content json in", time.time() -t1
    return web_return(as_json(by_region), response)
    
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
    return web_return(as_json(by_region), response)


@get('/repeat_region')
def repeat_regions():
    t1 = time.time()
    response.set_header("Access-Control-Allow-Origin","*")
    response.content_type = "application/json"
    
    by_region = []
    region = request.GET["region"]
    histogram = request.GET.get("histogram", None)
    interval = int(request.GET.get("interval", 1000))
    exclude_transcripts = True
    if histogram:
        print 'Histogram requested, but not available for repeat_regions'
    else:
        for reg in region.split(","):
            by_region.append({"id":reg, "result":get_repeats(reg)})

    print "**Returning GENE json in", time.time() -t1        
    return web_return(as_json(by_region), response)

@post('/search/')
def search(q=None):
    t1 = time.time()
    response.set_header("Access-Control-Allow-Origin","*")
    response.content_type = "application/json"
    seqid = request.POST.get('seqid', '')
    cl = sphinx.SphinxClient()
    cl.SetServer ("localhost", SPHINX_PORT)
    cl.SetMatchMode(sphinx.SPH_MATCH_ALL)
    cl.SetRankingMode(sphinx.SPH_RANK_SPH04)
    cl.SetSortMode(sphinx.SPH_SORT_RELEVANCE)
    cl.SetLimits (0, 50, 100)
    res = cl.Query("*%s*" %seqid)
    if not res:
        print 'query failed: %s' % cl.GetLastError()
        return ''
    
    if cl.GetLastWarning():
        print 'WARNING: %s\n' % cl.GetLastWarning()

    print 'Query \'%s\' retrieved %d of %d matches in %s sec' % (q, res['total'], res['total_found'], res['time'])
    print 'Query stats:'

    # result = {}
    # result["words"] = []
    # if res.has_key('words'):
    #     for info in res['words']:
    #         txt = '"%s" found %d times in %d entries' % (info['word'], info['hits'], info['docs'])
    #         result["words"].append(txt)
    #         print txt
            
    matches = {"results":[]}
    MATCHER = re.compile("(%s)"%q, re.IGNORECASE)
    if res.has_key('matches'):
        n = 1
        for match in res['matches']:
                locus, biotype, region = INDEX2REGION[int(match['id'])]
                desc = ', '.join(INDEX2DESC[int(match['id'])])
                chrom, _pos = region.split(':')
                start, end = _pos.split('-')
                matches['results'].append({"id":locus, "text":locus, "desc":desc, 'biotype':biotype,
                                           'chr':chrom,
                                           'start':start, 
                                           'end':end,
                                           'reg':region,
                                           'uniprot':GENE2UNIPROT.get(locus, ''),
                                           })
                n += 1
    print time.time() - t1, "search query" 
    return web_return(json.dumps(matches), response)

def get_repeats(region):
    t1 = time.time()
    ch, start, end = parse_region(region)
    genes = {}
    if start > len(SCAFFOLDS[ch]) or end < 0:
        return []

    rep_entries = FEATURES.get('repeat_region', []) + FEATURES.get('repeat_region2', [])
             
    for fch, fstart, fend, gene in rep_entries:
        if str(fch) == ch and (between(fstart, start, end) or between(fend, start, end)):
            gstrand = parse_strand(gene.location.strand)
            try:
                note = gene.qualifiers["note"][0]
                biotype = 'repeat_region'
            except KeyError:
                note = gene.qualifiers["Target"][0]
                biotype = 'repeat_region2'
                
            m = re.match("RepeatMasker,\s+Target\s+'([^']+)'\s+(\d+)\s(\d+)", note)
            if m:
                gname = "%s (%s-%s)" %(m.groups()[0], m.groups()[1], m.groups()[2])
            else:
                gname = note
            genedict = {"id": gname,
                        "biotype": biotype,
                        "featureType": "repeat",
                        "chromosome":ch,
                        "start":fstart,
                        "end":fend,
                        "strand":gstrand,
                        "transcripts":[]}
            genes[gname] = genedict
           
    result = sorted(genes.values(), lambda a,b: cmp(a["start"], b["start"]))             
    return result

def get_exons(region, biotypes, exclude_transcripts=False):
    ''' Returns all genes, transcripts and exons within a genomic region '''
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
                        "uniprot": GENE2UNIPROT.get(gname, ''),
                        "biotype": "gene",
                        "pep_img": int(os.path.exists(os.path.join(PEP_IMGS_DIR, "svg", gname+".svg"))),
                        "featureType": "gene",
                        "chromosome":ch,
                        "start":fstart,
                        "end":fend,
                        "strand":gstrand,
                        "source":"",
                        "annotations": additional_description(gname),
                        "expression": expr_info(gname),
                        "transcripts":[]}
            genes[gname] = genedict
            
    if len(genes) == 0:
        return []
    
    # populate transcripts
    exclude_biotypes = set(['repeat_region', 'misc_feature'])
    if not exclude_transcripts:
        for ftype in biotypes:
            if ftype in exclude_biotypes:
                continue
            for fch, fstart, fend, prot in FEATURES[ftype]:
                if str(fch) == ch and (between(fstart, start, end) or between(fend, start, end)):
                    pstrand = parse_strand(gene.location.strand)
                    try:
                        gname = prot.qualifiers["locus_tag"][0]
                    except KeyError:
                        print ' not link to LOCUS_TAG for biotype', ftype
                        continue
                    
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
    #print "   Number of final genes", len(genes), "completed in ", time.time() - t1, "secs"
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


def as_json(obj):
    return json.dumps({"response":obj})

def expr_info(locus_name):
    an = {}
    if locus_name in GENE2EXPR:
        for key, value in GENE2EXPR[locus_name].iteritems():
            an[key] = value

def additional_description(locus_name):
    an = {}
    if locus_name in GENE2DESC:
        for key, value in GENE2DESC[locus_name].iteritems():
            an[key] = value
    return an
 
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
    
# COMMAND LINE RELATED OPTIONS
    
def refresh_all_dbs():
    gbrecords, features, scaffolds, scaffold_genes = read_and_index_genome()
    gene2expr, gene2desc = read_and_index_extra_info()

    
    return gbrecords, features, scaffolds, scaffold_genes, gene2expr, gene2desc

def refresh_sphinx():
    index = []
    id2region = {}
    id2desc = {}
    i = 1
    for ftype in ['gene', 'CDS', 'tRNA', 'rRNA', 'misc_RNA',
                  'repeat_region', 'misc_feature']:
        for fch, fstart, fend, f in FEATURES[ftype]:
            gname = f.qualifiers.get("locus_tag", ["NoName"])[0]
            txt = []
            for key, value in f.qualifiers.iteritems():
                if key in AVOIDED_QUALIFIERS: continue
                if isinstance(value, list) and len(value) == 1:
                    txt.append(value[0])
                else:
                    txt.append('; '.join(map(str, value)))
            index.append([i, gname, GENE2UNIPROT.get(gname, ''), ', '.join(txt)])
            id2region[i] = [gname, ftype, "%s:%d-%d" %(fch, fstart, fend)]
            id2desc[i] = txt
            i += 1
            
    open(SPHINX_ANNOTATION_FILE, "w").write(
        '\n'.join(map(lambda x: '\t'.join(map(str, x)), index))+"\n")
    
    stop_sphinx()
    print 'Updating sphinx index... '
    system('cd %s && %s --config %s annotations' %(BASEPATH, INDEXER, SPHINX_CONFIG_FILE))
   
    return id2region, id2desc
        
def start_sphinx():
    stop_sphinx()
    print 'Starting sphinx sever... '
    s = system('cd %s && %s --config %s' %(BASEPATH, SEARCHD, SPHINX_CONFIG_FILE))
    if s:
        print color("Sphinx Start action Failed!", "red")
    else:
        print color("Sphinx Start action OK!", "green")

def stop_sphinx():
    print 'Stopping sphinx sever... '
    s = system('cd %s && %s --config %s --stop' %(BASEPATH, SEARCHD, SPHINX_CONFIG_FILE))
    if s:
        print color("Sphinx Stop action Failed!", "red")
    else:
        print color("Sphinx Stop action OK!", "green")

    
def start_webservices():
    # This starts the web service on the host and port specified. Currently it
    # should listen in "ct.bork.embl.de" and use the port 9000.
    run(host=HOST, port=WEBSERVICE_PORT, server='cherrypy')

def parse_gff(fname, ftype='unknown'):
    from BCBio import GFF
    entries = []
    for e in GFF.parse(fname):
        for f in e.features:
            entries.append([e.id, f.location.start, f.location.end, f])
    print len(entries), "entries read from GFF file", fname
    return entries
            
# Utils 
SHELL_COLORS = {
    "wr": '\033[1;37;41m', # white on red
    "wo": '\033[1;37;43m', # white on orange
    "wm": '\033[1;37;45m', # white on magenta
    "wb": '\033[1;37;46m', # white on blue
    "bw": '\033[1;37;40m', # black on white
    "lblue": '\033[1;34m', # light blue
    "lred": '\033[1;31m', # light red
    "lgreen": '\033[1;32m', # light green
    "yellow": '\033[1;33m', # yellow
    "cyan": '\033[36m', # cyan
    "blue": '\033[34m', # blue
    "green": '\033[32m', # green
    "orange": '\033[33m', # orange
    "red": '\033[31m', # red
    "magenta": "\033[35m", # magenta
    "white": "\033[0m", # white
    None: "\033[0m", # end
}
def color(string, color):
    return "%s%s%s" %(SHELL_COLORS[color], string, SHELL_COLORS[None])

class LogStdout(StringIO):
    def write(self, s):
        LOG.info(s.rstrip("\n"))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('action', metavar='start|stop|status', type=str, nargs=1, choices=["start", "stop", "status"],
                        help='')

    parser.add_argument('--refresh', dest='refresh', action='store_true',
                        help=('Refreshes all databases to account for changes in'
                        ' FASTA or GBF files, and start the daemon. This option will'
                        ' cause blast DB, search sphinx DB, and GBF annotation index'
                        ' to be updated.'))
    
    parser.add_argument('--host', dest='host', type=str,
                        help=('This option overwrites the HOST file information and'
                              ' starts the daemon listening in a different address.'))

    parser.add_argument('--search', dest='search', type=str,
                        help='test a query search')

    parser.add_argument('--daemon', dest='daemon', action = "store_true",
                        help='Start the service as a daemon')


    
    args = parser.parse_args()

    if args.search:
        search(args.search)
        sys.exit(0)
               
   
    # starts 
    try:
        HOST = open("HOST").readline().strip()
    except Exception:
        HOST = "localhost"
    if args.host:
        HOST = args.host

    action = args.action[0]
    if action == "stop":
        if os.path.exists(PID_FILE):
            pid = strip(open(PID_FILE).readline())
            system("kill %s" %pid)
        stop_sphinx()
                
    elif action == "start":
        if os.path.exists(PID_FILE):
            pid = strip(open(PID_FILE).readline())
            system("kill %s" %pid)
        stop_sphinx()
        
        if args.refresh:
            (GBRECORDS, FEATURES, SCAFFOLDS, SCAFFOLD_GENES, GENE2EXPR,
             GENE2DESC) = refresh_all_dbs()

            INDEX2REGION, INDEX2DESC = refresh_sphinx()
            
            # Saves the database for faster loading in the fugture
            cPickle.dump([GBRECORDS, FEATURES, SCAFFOLDS, SCAFFOLD_GENES, GENE2EXPR,
                          GENE2DESC, INDEX2REGION, INDEX2DESC],
                         open(DBFILE, "wb"))

        else:
            try:
                (GBRECORDS, FEATURES, SCAFFOLDS, SCAFFOLD_GENES, GENE2EXPR, GENE2DESC,
                 INDEX2REGION, INDEX2DESC) = cPickle.load(open(DBFILE))

            except Exception, e:
                print e
                print "Unable to load cached data. Try with --refresh"
                sys.exit(1)
        
        print FEATURES['repeat_region'][0]
        FEATURES['repeat_region2'] = parse_gff(REPEATS_GFF_FILE)
        sorted_sca = sorted(SCAFFOLDS.keys(),
                            lambda a,b: cmp(len(SCAFFOLD_GENES[a]), len(SCAFFOLD_GENES[b])),
                            reverse=True)


        # Sphinx will not start if another instance is running in the same
        # port. This call will try to stop any other running instance, but it might
        # fail if the service was started from a different user or directory. Just
        # 'pkill searchd' to kill any running process.
        start_sphinx()
        
        js = ['var scaffolds = [%s];\n' %(','.join(map(lambda x: '"%s"'%x, sorted_sca)))]
        js.append('var scaffolds_data = [')
        for sca_name in sorted_sca:
            size = len(SCAFFOLDS[sca_name])
            ngenes = len(SCAFFOLD_GENES[sca_name])
            js.append('{name:"%s",  cytobands: [], isCircular: 0, start: 1, end: %d, size: %d, numberGenes: %d },' \
                          %(sca_name, size, size, ngenes))
        js.append('];')
        open(JS_VARS_FILE, "w").write("/* Auto-generated file. Do not modify manually. */\n" + ' '.join(js))
        print color("Serving sequences and annotations for the following scaffolds:", 'green')
        print '\n'.join(map(lambda x: '  % 20s (%d sites, %d genes)'%(x, len(SCAFFOLDS[x]), len(SCAFFOLD_GENES[x])), sorted_sca))
        print color("The following annotations were found:", 'green')
        print '\n'.join(map(lambda x: "  % 20s %d" %(x[0], len(x[1])), sorted(FEATURES.items(), lambda a,b: cmp(len(a[1]), len(b[1])), reverse=True)))
        print color("Additional annotations:", 'green')
        print "  genes with expression data: ", len(GENE2EXPR)
        print "  genes with description data:", len(GENE2DESC)
        
        if args.daemon:
            # create logger with 'spam_application'
            fh = logging.FileHandler('log/ctbrowser_daemon.log')
            LOG.addHandler(fh)
            fh.setFormatter(WhitespaceRemovingFormatter())
            daemon = Daemonize(app="ctbrowser",
                               pid=PID_FILE,
                               action=start_webservices,
                               keep_fds=[3])
            sys.stdout = LogStdout()
            sys.stderr = LogStdout()
            print "Start!"
            
            daemon.start()
        else:
            print color("\n*** Running in debug mode. If you want to start the service in a production system, please use the --daemon flag.", "red")
            start_webservices()
            stop_sphinx()
            
    elif action == "status":
        if os.path.exists(PID_FILE):
            pid = strip(open(PID_FILE).readline())
            print "CT Browser daemon seems to be UP with PID %s" %pid
        else:
            print "CT Browser daemon seems to be DOWN"
            
        try:
            pid = strip(open(os.path.join(BASE_PATH, "searchd.pid")).readline())
        except Exception:
            print "CT Sphinx daemon seems to be DOWN"
        else:
            print "CT Sphinx daemon seems to be UP with PID %s" %pid
        
            
    
    
    
