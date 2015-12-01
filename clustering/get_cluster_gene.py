import gzip
import sys

def get_feature(fname, feature = "exon"):
    ss2gene = {}
    for ln in gzip.open(fname):
        if ln[0] == "#": continue
        ln = ln.split('\t')
        gID = ln[-1].split('gene_id "')[1].split('"')[0]
        if ln[2] != feature:
            continue
        ss2gene[(ln[0], int(ln[3]))] = gID
        ss2gene[(ln[0], int(ln[4]))] = gID
    return ss2gene


ss2gene = get_feature("gencode.v19.annotation.gtf.gz", "exon")

import numpy as np
W = file("%s.clu2gene.txt"%sys.argv[1].split("_perind")[0],'w')
for ln in gzip.open(sys.argv[1]):
    if "chrom" in ln: continue
    chrom, A, B, clu = ln.split()[0].split(":")
    gs = []
    vals = [int(x.split('/')[0])/(1+float(x.split('/')[1])) for x in ln.split()[1:] if x != "NA"]
    mean, median, minAS, maxAS = np.mean(vals), np.median(vals), min(vals), max(vals)
    if (chrom, int(A)) in ss2gene:
        gs.append(ss2gene[(chrom, int(A))])
    if (chrom, int(B)) in ss2gene:
        gs.append(ss2gene[(chrom, int(B))])

    
    if len(gs) > 0:
        W.write("%s %s %s %s %s %.2f %.2f %.2f %.2f\n"%(clu,chrom, A,B, gs[0], mean, median, minAS, maxAS))
    else:
        W.write("%s %s %s %s %s %.2f %.2f %.2f %.2f\n"%(clu,chrom, A,B, "?", mean, median, minAS, maxAS))
W.close()
