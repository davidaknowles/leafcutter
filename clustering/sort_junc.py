import sys
import gzip


runName = sys.argv[1]
chromLst = [chrom.strip() for chrom in  open(runName+"_chrom").readlines()]

merges = {}
for lib in sys.argv[2:]:
    libN = lib.replace("_1.",".").replace("_2.",".")
    if libN not in merges:
        merges[libN] = []
    merges[libN].append(lib)

fout_runlibs = file(runName+"_sortedlibs",'w')

for libN in merges:

    by_chrom = {}
    foutName = libN+'.sorted.gz'

    fout_runlibs.write(foutName+'\n')

    try: gzip.open(foutName)
    except:
        sys.stderr.write("Sorting %s.."%libN)
    else:
        continue

    if len(merges[libN]) > 1:
        sys.stderr.write("merging %s...\n"%(" ".join(merges[libN])))
    else:
        sys.stderr.write("\n")

    for lib in merges[libN]:
        
        for ln in open(lib):
        

            lnsplit=ln.split()
            if len(lnsplit)<6: 
                sys.stderr.write("Error in %s \n" % lib)
                continue
            chrom, A, B, dot, counts, strand = lnsplit
            
            if chrom not in chromLst: continue
            if chrom not in by_chrom:
                by_chrom[chrom] = {}
            intron = (int(start), int(end), strand)
            if intron in by_chrom[chrom]:
                by_chrom[chrom][intron] += int(count)
            else:
                by_chrom[chrom][intron] = int(count)

    ks = by_chrom.keys()
    ks.sort()
    
    fout = gzip.open(foutName,'w')
    for k in ks:
        sortedIntrons = by_chrom[k].keys()
        sortedIntrons.sort()
        for start, end, strand in sortedIntrons:
            count = by_chrom[k][(start,end,strand)]
            fout.write("\t".join([k, str(start), str(end), '.', str(count), strand])+'\n')
    fout.close()

fout_runlibs.close()
