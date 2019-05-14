#!/usr/bin/env python

import sys
import gzip
import numpy as np
import scipy as sc
import pickle

from optparse import OptionParser
    
from sklearn.decomposition import PCA
from sklearn import preprocessing
from sklearn import linear_model

from scipy.stats import rankdata
from scipy.stats import norm

def qqnorm(x):
    n=len(x)
    a=3.0/8.0 if n<=10 else 0.5
    return(norm.ppf( (rankdata(x)-a)/(n+1.0-2.0*a) ))

def stream_table(f, ss = ''):
    fc = '#'
    while fc[0] == "#":
        fc = f.readline().strip()
        head = fc.split(ss)

    for ln in f:
        ln = ln.strip().split(ss)
        attr = {}

        for i in range(len(head)):
            try: attr[head[i]] = ln[i]
            except: break
        yield attr

def get_chromosomes(ratio_file):
    """Get chromosomes from table. Returns set of chromosome names"""
    try: open(ratio_file)
    except:
        sys.stderr.write("Can't find %s..exiting\n"%(ratio_file))
        return
    sys.stderr.write("Parsing chromosome names...\n")
    chromosomes = set()
    with gzip.open(ratio_file, 'rt') as f:
            f.readline()
            for line in f:
                chromosomes.add(line.split(":")[0])
    return(chromosomes)

def get_blacklist_chromosomes(chromosome_blacklist_file):
    """
    Get list of chromosomes to ignore from a file with one blacklisted
    chromosome per line. Returns list. eg. ['X', 'Y', 'MT']
    """
    if chromosome_blacklist_file:
        with open(chromosome_blacklist_file, 'r') as f:
            return(f.read().splitlines())
    else:
        return(["X", "Y"])



def main(ratio_file, chroms, blacklist_chroms, pcs=50):
    
    dic_pop, fout = {}, {}
    try: open(ratio_file)
    except: 
        sys.stderr.write("Can't find %s..exiting\n"%(ratio_file))
        return

    sys.stderr.write("Starting...\n")
    for i in chroms:
        fout[i] = file(ratio_file+".phen_chr"+i,'w')
        fout_ave = file(ratio_file+".ave",'w')
    valRows, valRowsnn, geneRows = [], [], []
    finished = False
    header = gzip.open(ratio_file).readline().split()[1:]


    for i in fout:
        fout[i].write("\t".join(["#Chr","start", "end", "ID"]+header)+'\n')

    for dic in stream_table(gzip.open(ratio_file),' '):

        chrom = dic['chrom'].replace("chr",'')
        chr_ = chrom.split(":")[0]
        if chr_ in blacklist_chroms: continue
        NA_indices, valRow, aveReads = [], [], []
        tmpvalRow = []

        i = 0
        for sample in header:

            try: count = dic[sample]
            except: print chrom, len(dic)
            num, denom = count.split('/')
            if float(denom) < 1:
                count = "NA"
                tmpvalRow.append("NA")
                NA_indices.append(i)
            else:
                # add a 0.5 pseudocount
                count = (float(num)+0.5)/((float(denom))+0.5)
                tmpvalRow.append(count) 
                aveReads.append(count)


        # If ratio is missing for over 40% of the samples, skip
        if tmpvalRow.count("NA") > len(tmpvalRow)*0.4:
            continue

        ave = np.mean(aveReads)

        # Set missing values as the mean of all values
        for c in tmpvalRow:
            if c == "NA": valRow.append(ave)
            else: valRow.append(c)

        # If there is too little variation, skip (there is a bug in fastqtl which doesn't handle cases with no variation)
        if np.std(valRow) < 0.005: continue

        chr_, s, e, clu = chrom.split(":")
        if len(valRow) > 0:                
            fout[chr_].write("\t".join([chr_,s,e,chrom]+[str(x) for x in valRow])+'\n')
            fout_ave.write(" ".join(["%s"%chrom]+[str(min(aveReads)), str(max(aveReads)), str(np.mean(aveReads))])+'\n')

            # scale normalize
            valRowsnn.append(valRow)                
            valRow = preprocessing.scale(valRow)

            valRows.append(valRow)
            geneRows.append("\t".join([chr_,s,e,chrom]))
            if len(geneRows) % 1000 == 0:
                sys.stderr.write("Parsed %s introns...\n"%len(geneRows))
                
    for i in fout:
        fout[i].close()

    # qqnorms on the columns
    matrix = np.array(valRows)
    for i in xrange(len(matrix[0,:])):
        matrix[:,i] = qqnorm(matrix[:,i])
        
    # write the corrected tables
    fout = {}
    for i in chroms:
        fn="%s.qqnorm_chr%s"%(ratio_file,i)
        print("Outputting: " + fn)
        fout[i] = file(fn,'w')
        fout[i].write("\t".join(['#Chr','start','end','ID'] + header)+'\n')
    lst = []
    for i in xrange(len(matrix)):
        chrom, s = geneRows[i].split()[:2]
        
        lst.append((chrom.replace("chr",""), int(s), "\t".join([geneRows[i]] + [str(x) for x in  matrix[i]])+'\n'))

    lst.sort()
    for ln in lst:
        fout[ln[0]].write(ln[2])
        
    fout_run = file("%s_prepare.sh"%ratio_file,'w')

    for i in fout:
        fout[i].close()
        fout_run.write("bgzip -f %s.qqnorm_chr%s\n"%(ratio_file, i))
        fout_run.write("tabix -p bed %s.qqnorm_chr%s.gz\n"%(ratio_file, i))
    fout_run.close()

    sys.stdout.write("Use `sh %s_prepare.sh' to create index for fastQTL (requires tabix and bgzip).\n"%ratio_file)

    if pcs>0:
        #matrix = np.transpose(matrix) # important bug fix (removed as of Jan 1 2018)
        pcs = min([len(header), pcs])
        pca = PCA(n_components=pcs)                                                                                                                                                                            
        pca.fit(matrix)  
        pca_fn=ratio_file+".PCs"
        print("Outputting PCs: " + pca_fn)
        pcafile = file(pca_fn,'w')  
        pcafile.write("\t".join(['id']+header)+'\n')
        pcacomp = list(pca.components_)
    
        for i in range(len(pcacomp)):
            pcafile.write("\t".join([str(i+1)]+[str(x) for x in pcacomp[i]])+'\n')

        pcafile.close()

if __name__ == "__main__":

    parser = OptionParser(usage="usage: %prog [-p num_PCs] input_perind.counts.gz")
    parser.add_option("-p", "--pcs", dest="npcs", default = 50, help="number of PCs output")
    parser.add_option("--ChromosomeBlackList", dest="cbl", default="", help="file of blacklisted chromosomes to exclude from analysis, one per line. If none is provided, will default to blacklisting X and Y")
    (options, args) = parser.parse_args()
    if len(args)==0:
        sys.stderr.write("Error: no ratio file provided... (e.g. python leafcutter/scripts/prepare_phenotype_table.py input_perind.counts.gz\n")
        exit(0)
    main(args[0], get_chromosomes(args[0]), get_blacklist_chromosomes(options.cbl), int(options.npcs) )
    
