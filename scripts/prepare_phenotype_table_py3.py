#!/usr/bin/env python
# written by Yang Li
# ported to python 3 by Jack Humphrey
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

def main(ratio_file, pcs=50):
    
    dic_pop, fout = {}, {}
    try: open(ratio_file)
    except: 
        sys.stderr.write("Can't find %s..exiting\n"%(ratio_file))
        return

    sys.stderr.write("Starting...\n")
    for i in range(1,23):
        fout[i] = open(ratio_file+".phen_chr%d"%i,'w')
        fout_ave = open(ratio_file+".ave",'w')
    valRows, valRowsnn, geneRows = [], [], []
    finished = False
    # have to now specify text mode
    header = gzip.open(ratio_file, mode = 'rt').readline().split()[1:]

    # problematic line - something about py2->py3 and handling gzipped files
    for i in fout:
        header_string = "\t".join(["#Chr","start", "end", "ID"]+header)+'\n'
        fout[i].write(header_string) 
	#fout[i].write("\t".join(["#Chr","start", "end", "ID"]+header)+'\n')
	
    for dic in stream_table(gzip.open(ratio_file, mode = 'rt'),' '):

        chrom = dic['chrom'].replace("chr",'')
        chr_ = chrom.split(":")[0]
        if chr_ in 'XY': continue
        NA_indices, valRow, aveReads = [], [], []
        tmpvalRow = []

        i = 0
        for sample in header:

            try: count = dic[sample]
            except: print(chrom, len(dic))
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
            chrom_int = int(chr_)
            fout[chrom_int].write("\t".join([chr_,s,e,chrom]+[str(x) for x in valRow])+'\n')
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
    for i in range(len(matrix[0,:])):
        matrix[:,i] = qqnorm(matrix[:,i])
        
    # write the corrected tables
    fout = {}
    for i in range(1,23):
        fn="%s.qqnorm_chr%d"%(ratio_file,i)
        print(("Outputting: " + fn))
        fout[i] = open(fn,'w')
        fout[i].write("\t".join(['#Chr','start','end','ID'] + header)+'\n')
    lst = []
    for i in range(len(matrix)):
        chrom, s = geneRows[i].split()[:2]
        
        lst.append((int(chrom.replace("chr","")), int(s), "\t".join([geneRows[i]] + [str(x) for x in  matrix[i]])+'\n'))

    lst.sort()
    for ln in lst:
        fout[ln[0]].write(ln[2])
        
    fout_run = open("%s_prepare.sh"%ratio_file,'w')

    for i in fout:
        fout[i].close()
        fout_run.write("bgzip -f %s.qqnorm_chr%d\n"%(ratio_file, i))
        fout_run.write("tabix -p bed %s.qqnorm_chr%d.gz\n"%(ratio_file, i))
    fout_run.close()

    sys.stdout.write("Use `sh %s_prepare.sh' to create index for fastQTL (requires tabix and bgzip).\n"%ratio_file)

    if pcs>0:
        #matrix = np.transpose(matrix) # important bug fix (removed as of Jan 1 2018)
        pcs = min([len(header), pcs])
        pca = PCA(n_components=pcs)                                                                                                                                                                            
        pca.fit(matrix)  
        pca_fn=ratio_file+".PCs"
        print(("Outputting PCs: " + pca_fn))
        pcafile = open(pca_fn,'w')  
        pcafile.write("\t".join(['id']+header)+'\n')
        pcacomp = list(pca.components_)
    
        for i in range(len(pcacomp)):
            pcafile.write("\t".join([str(i+1)]+[str(x) for x in pcacomp[i]])+'\n')

        pcafile.close()

if __name__ == "__main__":

    parser = OptionParser(usage="usage: %prog [-p num_PCs] input_perind.counts.gz")
    parser.add_option("-p", "--pcs", dest="npcs", default = 50, help="number of PCs output")
    (options, args) = parser.parse_args()
    if len(args)==0:
        sys.stderr.write("Error: no ratio file provided... (e.g. python leafcutter/scripts/prepare_phenotype_table.py input_perind.counts.gz\n")
        exit(0)
    main(args[0], int(options.npcs) )
    
