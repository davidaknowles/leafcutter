import sys
import gzip


runName = sys.argv[1]
chromLst = [chrom.strip() for chrom in  open(runName+"_chrom").readlines()]
refined_cluster = runName + "_refined"

exons, cluExons = {}, {}
cluN = 0

for ln in open(refined_cluster):
    chrom = ln.split()[0]
    cluN += 1
    for exon in ln.split()[1:]:
        A, B, count = exon.split(":")
        if chrom not in exons:
            exons[chrom] = {}

        exons[chrom][(int(A),int(B))] = cluN
        if cluN not in cluExons:
            cluExons[cluN] = []
        cluExons[cluN].append((chrom, A, B))

merges = {}
for lib in sys.argv[2:]:
    libN = lib.replace("_1.",".").replace("_2.",".")
    if libN not in merges:
        merges[libN] = []
    merges[libN].append(lib)


fout_runlibs = file(runName+"_sortedlibs",'w')

for libN in merges:

    by_chrom = {}
    foutName = libN+'.%s.sorted.gz'%(runName.split("/")[-1])

    fout_runlibs.write(foutName+'\n')

    try: gzip.open(foutName)
    except:
        pass
    else:
        #continue
        pass
    sys.stderr.write("Sorting %s..\n"%libN)
    if len(merges[libN]) > 1:
        sys.stderr.write("merging %s...\n"%(" ".join(merges[libN])))
    else:
        #sys.stderr.write("")
        pass
    fout = gzip.open(foutName,'w')

    fout.write("chrom %s\n"%libN.split("/")[-1].split(".junc")[0])

    for lib in merges[libN]:
        
        for ln in open(lib):
            chrom, start, end, dot, count, strand = ln.split()
            if chrom not in chromLst: continue
            if chrom not in by_chrom:
                by_chrom[chrom] = {}
            intron = (int(start), int(end)+1)
            if intron in by_chrom[chrom]:
                by_chrom[chrom][intron] += int(count)
            else:
                by_chrom[chrom][intron] = int(count)

    for clu in cluExons:
        buf = []
        ks = cluExons[clu]
        ks.sort()
        tot = 0
        for exon in ks:
            chrom, start, end = exon
            start, end = int(start), int(end)
            if (start,end) in by_chrom[chrom]:
                tot += by_chrom[chrom][(start,end)]
        for exon in ks:
            
            chrom, start, end = exon
            start, end = int(start), int(end)
            if (start,end) in by_chrom[chrom]:
                #print chrom, start, end, by_chrom[chrom][(start,end)]
                buf.append("%s:%d:%d:clu_%d %d/%d\n"%(chrom,start, end, clu, by_chrom[chrom][(start,end)], tot))
            else:
                buf.append("%s:%d:%d:clu_%d 0/%d\n"%(chrom,start, end,clu, tot))
        #print chrom, " ".join(buf)
        fout.write("".join(buf))
    fout.close()
fout_runlibs.close()
