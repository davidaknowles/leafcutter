#!/usr/bin/env python

import gzip
import sys


def getmean(lst):
    return sum(lst)/float(len(lst))

def getmedian(lst):
    lst.sort()
    n = len(lst)
    if n < 1:
        return None
    if n % 2 == 1:
        return lst[n//2]
    else:
        return sum(lst[n//2-1:n//2+1])/2.0

def get_feature(fname, feature = "exon"):
    ss2gene = {}
    if ".gz" in fname:
        F = gzip.open(fname)
    else:
        F = open(fname)
    for ln in F:
        if ln[0] == "#": continue
        ln = ln.split('\t')
        gID = ln[-1].split('gene_id "')[1].split('"')[0]
        if ln[2] != feature:
            continue
        ss2gene[(ln[0], int(ln[3]))] = gID
        ss2gene[(ln[0], int(ln[4]))] = gID
    return ss2gene


ss2gene = get_feature(sys.argv[1], "exon")

W = file("%s.clu2gene.txt"%sys.argv[2].split("_perind")[0],'w')
for ln in gzip.open(sys.argv[2]):
    if "chrom" in ln: continue
    
    if len(ln.split()[0].split(":")) == 5:
        chrom, A, B, clu, strand = ln.split()[0].split(":")
    else:
        chrom, A, B, clu = ln.split()[0].split(":")
    gs = []
    vals = [int(x.split('/')[0])/(1+float(x.split('/')[1])) for x in ln.split()[1:] if x != "NA"]
    mean, median, minAS, maxAS = getmean(vals), getmedian(vals), min(vals), max(vals)
    if (chrom, int(A)) in ss2gene:
        gs.append(ss2gene[(chrom, int(A))])
    if (chrom, int(B)) in ss2gene:
        gs.append(ss2gene[(chrom, int(B))])

    
    if len(gs) > 0:
        W.write("%s %s %s %s %s %.2f %.2f %.2f %.2f\n"%(clu,chrom, A,B, gs[0], mean, median, minAS, maxAS))
    else:
        W.write("%s %s %s %s %s %.2f %.2f %.2f %.2f\n"%(clu,chrom, A,B, "?", mean, median, minAS, maxAS))
W.close()
