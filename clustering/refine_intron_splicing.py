




def main(outPrefix, minReads=50):

    inFile = outPrefix + "_pooled"
    outFile = outPrefix + "_refined"

    W = file(outFile,'w')
    Ncl = 0
    for ln in open(inFile):
        clu = []
        totN = 0
        chrom = ln.split()[0]
        for ex in ln.split()[1:]:
            A, B, N = ex.split(":")
            clu.append(((int(A),int(B)), int(N)))
            totN += int(N)
        if totN < minReads: continue
        rc = refine_cluster(clu)
        if len(rc) > 0:
            #print "clu", len(rc)
            for clu in rc:
                buf = '%s ' % chrom
                #print clu
                for interval, count in clu:
                    #print interval, count
                    buf += "%d:%d" % interval + ":%d"%(count)+ " "
                Ncl += 1
                W.write(buf+'\n')
    sys.stderr.write("Split into %s clusters...\n"%Ncl)
    W.close()

def cluster_intervals(E):
    '''                                                                                                                                                                                                                                                                                                     
    Clusters intervals together.                                                                                                                                                                                                                                                                            
    '''
    E.sort()
    current = E[0]
    Eclusters = []
    Edic = {}
    cluster = []

    i = 0
    while i < len(E):

        if overlaps(E[i], current):
            cluster.append(E[i])
        else:
            Eclusters.append(cluster)
            for k in set(cluster):
                Edic[k] = set(cluster)

            cluster = [E[i]]
        current = (E[i][0], max([current[1], E[i][1]]))
        i += 1

    if len(cluster) > 0:
        Eclusters.append(cluster)
        for k in set(cluster):
            Edic[k] = set(cluster)

    return Eclusters, Edic


def overlaps(A,B):
    '''                                                                                                                                                                                                                                                                                                    
    Checks if A and B overlaps                                                                                                                                                                                                                                                                             
    '''

    if A[1] < B[0] or B[1] < A[0]:
        return False
    else: return True

def refine_cluster(clu, cutoff = 0.01):
    ''' for each exon in the cluster compute the ratio of reads, if smaller than 0.1,
    remove and recluster '''
    
    remove = []
    dic = {}
    intervals = []

    totN = 0
    for inter, count in clu:
        totN += count
    for inter, count in clu:
        if (count/float(totN) >= cutoff and count >= 10) or count >= 20:
            intervals.append(inter)
            dic[inter] = count
            
    if len(intervals) == 0: return []
    A, B = cluster_intervals(intervals)
    
    if len(A) == 1:
        rc = [(x, dic[x]) for x in A[0]]
        if len(rc) > 1:
            return [[(x, dic[x]) for x in A[0]]]
        else:
            return []
    NCs = []
    for c in A:
        if len(c) > 1:
            NC = refine_cluster([(x, dic[x]) for x in c])
            NCs += NC
    return NCs

if __name__ == "__main__":
    import sys
    main(sys.argv[1],int(sys.argv[2]))
