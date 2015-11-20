def main(outPrefix, maxIntronLen, flist, ov_cutoff = 6):

    outFile = outPrefix+"_pooled"
    
    chromLst = [chrom.strip() for chrom in  open(outPrefix+"_chrom").readlines()]
    
    by_chrom = {}
    for lib in flist:
        
        lib = lib.strip()
        sys.stderr.write("scanning %s...\n"%lib)

        for ln in open(lib):
            chrom, A, B, dot, counts, strand = ln.split()
            if chrom not in chromLst: continue
            A, B = int(A), int(B)+1
            if B-A > int(maxIntronLen): continue
            try: by_chrom[chrom][(A,B)] = int(counts) + by_chrom[chrom][(A,B)]
            except: 
                try: by_chrom[chrom][(A,B)] = int(counts)
                except: by_chrom[chrom] = {(A,B):int(counts)}

    W = file(outFile, 'w')
    N = 0
    sys.stderr.write("Parsing ")
    for chrom in by_chrom:
        read_ks = by_chrom[chrom].keys()
        read_ks.sort()
        sys.stderr.write("%s.."%chrom)
        clu = cluster_intervals(read_ks)
        for cl in clu:
            if len(cl) > 1:
                buf = '%s '%chrom
                
                for c in refine_linked(cl):
                    if len(c) > 1:
                        buf = '%s '%chrom
                        for interval, count in [(x, by_chrom[chrom][x]) for x in c]:
                            buf += "%d:%d" % interval + ":%d"%count+ " "
                        W.write(buf+'\n')
                N += 1
    sys.stderr.write("\nWrote %d clusters..\n"%N)
    W.close()

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

def overlaps(A,B):
    '''
    Checks if A and B overlaps
    '''

    if A[1] < B[0] or B[1] < A[0]:
        return False
    else: return True


def refine_linked(clusters):
    all_clus = []
    while len(clusters) > 0:
        start = clusters[0]
        clusters = clusters[1:]
        current_pos = [start[0],start[1]]
        current_clu = [start]
        while True:
            torm = []
            for intron in clusters:
                if intron[0] in current_pos or intron[1] in current_pos:
                    current_clu.append(intron)
                    current_pos += [intron[0],intron[1]]
                    torm.append(intron)
            for intron in torm:
                clusters.remove(intron)
            if len(torm) == 0: break
            
                
        all_clus.append(current_clu)
    return all_clus


if __name__ == "__main__":
    import sys
    
    main(sys.argv[1], sys.argv[2], sys.argv[3:])
