import sys
import gzip
import copy
import numpy as np


def merge_dicts(x, y):
    '''Given two dicts, merge them into a new dict as a shallow copy.'''
    z = x.copy()
    z.update(y)
    return z

def overlaps(A,B):
    '''
    Checks if A and B overlaps
    '''

    if A[1] < B[0] or B[1] < A[0]:
        return False
    else: return True

def cluster_intervals(E):
    '''                                                                                                                                                                                                     
    Clusters intervals together.                                                                                                                                                                            
    '''
    E.sort()
    current = E[0]
    Eclusters = []
    
    cluster = []

    i = 0
    while i < len(E):

        if overlaps(E[i], current):
            cluster.append(E[i])
        else:
            Eclusters.append(cluster)
            cluster = [E[i]]
        current = (E[i][0], max([current[1], E[i][1]]))
        i += 1

    if len(cluster) > 0:
        Eclusters.append(cluster)

    return Eclusters


# get splicing events that matter

allClu, exons, cluMaxEnd, cluChrom, cluExons = {}, {}, {}, {}, {}
cluN = 0


outPrefix = sys.argv[1]
refined_cluster = outPrefix + "_refined"
outName = outPrefix + "_perind"
libFileNames = [lib.strip() for lib in open(outPrefix+'_sortedlibs').readlines()]

for ln in open(refined_cluster):
    chrom = ln.split()[0]
    cluN += 1
    for exon in ln.split()[1:]:
        A, B, count = exon.split(":")
        if chrom not in exons:
            exons[chrom] = {}
            cluChrom[chrom] = []
        exons[chrom][(int(A),int(B))] = cluN
        if cluN not in cluExons:
            cluExons[cluN] = []
        cluExons[cluN].append((chrom, A, B))
        if cluN not in cluMaxEnd:
            cluChrom[chrom].append(cluN)
            cluMaxEnd[cluN] = (int(A),int(B))
        else:
            cluMaxEnd[cluN] = (min([int(A), cluMaxEnd[cluN][0]]), max([int(B), cluMaxEnd[cluN][1]]))



libs = []
by_chrom = {}


ks = exons.keys()
ks.sort()


libs = []
libFiles = {}
for lib in libFileNames:
    libFiles[lib] = gzip.open(lib)
    libs.append(lib)


# Now cluster clusters to determine the ends of scanning phase

chromCuts = {}
skip = 150
for chrom in cluChrom:
    
    clusters = []
    for clu in cluChrom[chrom]:
        clusters.append(cluMaxEnd[clu])
    N = 0
    for clu in cluster_intervals(clusters):
        N += 1
        if N % skip == 0:
            if chrom not in chromCuts:
                chromCuts[chrom] = []
            chromCuts[chrom].append(max([intv[1] for intv in clu]))

# Chromcuts tells us where we can stop reading

fout_counts = gzip.open(outName+'.counts.gz','w')
fout_counts.write("chrom "+" ".join([name.split("/")[-1].split('.')[0] for name in libs])+'\n')

index, cChrom = -1, None
finished = False

libChrom = {}

dic = {}
dic2 = {}
dics_leftover = []
finished = {}
while len(finished) < len(libs):
    
    sys.stderr.write(".")
    for lib in libs:
        # For each library read until readTo and then start writing
        if lib in libChrom:
            # Do not start reading a file until all files are at the same cChrom
            if libChrom[lib] != cChrom: continue
        while True:
            ln = libFiles[lib].readline()
            
            if len(ln) == 0:
                finished[lib] = ''
                break
            
            chrom, start, end, dot, count, strand = ln.split()
            libChrom[lib] = chrom
            
            end = str(int(end)+1)
            if cChrom == None:
                readTo = chromCuts[chrom][0]
                index = 0
                cChrom = chrom
            else:
                if index >= len(chromCuts[chrom]):
                    readTo = np.inf
                else:
                    readTo = chromCuts[chrom][index]
            
            if int(start) > readTo or index == -1 or chrom != cChrom:
                # This exon's data must be stored for later use 
                if (int(start), int(end)) in exons[chrom]: 
                    cluN = exons[chrom][(int(start), int(end))]
                    if cluN not in dic2:
                        dic2[cluN] = {}
                    if lib not in dic2[cluN]:
                        dic2[cluN][lib] = {}
                    dic2[cluN][lib][(chrom, start, end)] = count
                break
            
            if (int(start), int(end)) not in exons[chrom]: continue

            cluN = exons[chrom][(int(start), int(end))]
            if cluN not in dic:
                dic[cluN] = {}
            if lib not in dic[cluN]:
                dic[cluN][lib] = {}
            dic[cluN][lib][(chrom, start, end)] = count
            
    chroms = libChrom.values()

    for clu in dic:
    
        libCounts = {}
        allClu[clu] = ''
        for lib in dic[clu]:
            libCounts[lib] = 0
            for intron in cluExons[clu]:
                if intron in dic[clu][lib]:
                    libCounts[lib] += int(dic[clu][lib][intron])        
        
        for intron in cluExons[clu]:
            buf = ["%s:%s:%s:"%intron+"clu_%s"%clu]
            for lib in libs:
                
                if lib not in dic[clu]:
                    buf.append("0/0")
                elif intron not in dic[clu][lib]:
                    buf.append("0/%d"%libCounts[lib])
                else:
                    buf.append("%s/%d"%(dic[clu][lib][intron],libCounts[lib]))
            fout_counts.write(" ".join(buf)+'\n')

    if (len(set(chroms)) == 1 and cChrom != chroms[0]) or len(finished) == len(libs):
        sys.stderr.write("%s done.\n"%cChrom)
        cChrom = chroms[0]
        index = -1
        
        dic = copy.deepcopy(dic2)
        #break
        for d in dics_leftover:
            dic = merge_dicts(dic, d)
        dic2 = {}
        dics_leftover = []
    elif len(set(chroms)) >= 2:
        # some junction files ran out things
        dics_leftover.append(dic2)
    else:
        dic = copy.deepcopy(dic2)
        dic2 = {}
        dics_leftover = []
    index += 1

if len(dic) > 0:
    for clu in dic:
        libCounts = {}
        allClu[clu] = ''
        for lib in dic[clu]:
            libCounts[lib] = 0
            for intron in cluExons[clu]:
                if intron in dic[clu][lib]:
                    libCounts[lib] += int(dic[clu][lib][intron])

        for intron in cluExons[clu]:
            buf = ["%s:%s:%s:"%intron+"clu_%s"%clu]
            for lib in libs:
                if lib not in dic[clu]:
                    buf.append("0/0")
                elif intron not in dic[clu][lib]:    
                    buf.append("0/%d"%libCounts[lib])
                else:
                    buf.append("%s/%d"%(dic[clu][lib][intron],libCounts[lib]))
            fout_counts.write(" ".join(buf)+'\n')

sys.stderr.write("Summarized %d clusters\n"% len(allClu))

fout_counts.close()
