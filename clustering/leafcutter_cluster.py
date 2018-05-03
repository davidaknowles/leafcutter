#!/usr/bin/env python

import sys
import tempfile
import os
import gzip
import shutil

def main(options,libl):
    
    pool_junc_reads(libl, options)
    refine_clusters(options)
    sort_junctions(libl, options)
    merge_junctions(options)
    get_numers(options)

def pool_junc_reads(flist, options):

    outPrefix = options.outprefix
    rundir = options.rundir
    maxIntronLen = int(options.maxintronlen)
    checkchrom = options.checkchrom
    useStrand = options.strand

    outFile = "%s/%s_pooled"%(rundir,outPrefix)
    
    chromLst = ["chr%d"%x for x in range(1,23)]+['chrX','chrY']+["%d"%x for x in range(1,23)]+['X','Y']
    by_chrom = {}
    for libl in flist:
        
        lib = libl.strip()
        if not os.path.isfile(lib):
            continue

        if options.verbose:
            sys.stderr.write("scanning %s...\n"%lib)

        if lib[-3:] == ".gz": F = gzip.open(lib)
        else: F = open(lib)
        for ln in F:
            
            lnsplit=ln.split()
            if len(lnsplit)<6: 
                sys.stderr.write("Error in %s \n" % lib)
                continue
            chrom, A, B, dot, counts, strand = lnsplit
            
            if not useStrand:
                strand = "NA"
            if checkchrom and (chrom not in chromLst): continue
            A, B = int(A), int(B)+1
            if B-A > int(maxIntronLen): continue
            try: by_chrom[(chrom,strand)][(A,B)] = int(counts) + by_chrom[(chrom, strand)][(A,B)]
            except: 
                try: by_chrom[(chrom,strand)][(A,B)] = int(counts)
                except: by_chrom[(chrom, strand)] = {(A,B):int(counts)}

    fout = file(outFile, 'w')
    Ncluster = 0
    sys.stderr.write("Parsing...\n")
    for chrom in by_chrom:
        read_ks = [k for k,v in by_chrom[chrom].items() if v >= 3] # a junction must have at least 3 reads
        read_ks.sort()
        sys.stderr.write("%s:%s.."%chrom)
        clu = cluster_intervals(read_ks)[0]
        for cl in clu:
            if len(cl) > 1: # if cluster has more than one intron  
                buf = '%s:%s '%chrom
                for interval, count in [(x, by_chrom[chrom][x]) for x in cl]:
                    buf += "%d:%d" % interval + ":%d"%count+ " "
                fout.write(buf+'\n')
            Ncluster += 1
    sys.stderr.write("\nWrote %d clusters..\n"%Ncluster)
    fout.close()


def sort_junctions(libl, options):

    chromLst = ["chr%d"%x for x in range(1,23)]+['chrX','chrY']+["%d"%x for x in range(1,23)]+['X','Y'] 
    outPrefix = options.outprefix
    rundir = options.rundir
    refined_cluster = "%s/%s_refined"%(rundir,outPrefix)
    checkchrom = options.checkchrom
    useStrand = options.strand
    runName = "%s/%s"%(rundir, outPrefix)


    exons, cluExons = {}, {}
    cluN = 0

    for ln in open(refined_cluster):
        chrom = ln.split()[0]
        chrom = tuple(chrom.split(":"))
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
    for ll in libl:
        lib=ll.rstrip()
        if not os.path.isfile(lib):
            continue
        libN = lib
        if libN not in merges:
            merges[libN] = []
        merges[libN].append(lib)

    fout_runlibs = file(runName+"_sortedlibs",'w')

    for libN in merges:
        libName = "%s/%s"%(rundir,libN.split('/')[-1])
        by_chrom = {}
        foutName = libName+'.%s.sorted.gz'%(runName.split("/")[-1])

        fout_runlibs.write(foutName+'\n')

        if options.verbose:   
            sys.stderr.write("Sorting %s..\n"%libN)
        if len(merges[libN]) > 1:
            if options.verbose:   
                sys.stderr.write("merging %s...\n"%(" ".join(merges[libN])))
        else:
            pass
        fout = gzip.open(foutName,'w')

        fout.write("chrom %s\n"%libN.split("/")[-1].split(".junc")[0])

        for lib in merges[libN]:
            if lib[-3:] == ".gz": F = gzip.open(lib)
            else: F = open(lib)
            for ln in F:

                lnsplit=ln.split()
                if len(lnsplit)<6: 
                    sys.stderr.write("Error in %s \n" % lib)
                    continue
                chrom, start, end, dot, count, strand = ln.split()
                if not useStrand:
                    strand = "NA"
                if checkchrom and (chrom not in chromLst): continue
                chrom = (chrom, strand)
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
                if chrom not in by_chrom:
                    pass
                elif (start,end) in by_chrom[chrom]:
                    tot += by_chrom[chrom][(start,end)]
            for exon in ks:
            
                chrom, start, end = exon
                start, end = int(start), int(end)
                chromID, strand = chrom
                if chrom not in by_chrom:
                    buf.append("%s:%d:%d:clu_%d_%s 0/%d\n"%(chromID,start, end,clu, strand, tot))
                elif (start,end) in by_chrom[chrom]:                
                    buf.append("%s:%d:%d:clu_%d_%s %d/%d\n"%(chromID,start, end, clu,strand, by_chrom[chrom][(start,end)], tot))
                else:
                    buf.append("%s:%d:%d:clu_%d_%s 0/%d\n"%(chromID,start, end,clu,strand, tot))
        
            fout.write("".join(buf))
        fout.close()
    fout_runlibs.close()

def refine_clusters(options):
    
    outPrefix = options.outprefix
    rundir = options.rundir
    minratio = float(options.mincluratio)
    minreads = int(options.minclureads)

    inFile = "%s/%s_pooled"%(rundir,outPrefix)
    outFile = "%s/%s_refined"%(rundir,outPrefix)

    fout = file(outFile,'w')
    Ncl = 0
    for ln in open(inFile):
        clu = []
        totN = 0
        chrom = ln.split()[0]
        for ex in ln.split()[1:]:
            A, B, N = ex.split(":")
            clu.append(((int(A),int(B)), int(N)))
            totN += int(N)
        if totN < minreads: continue
        #print "CLU",clu
        #print "linked",refine_linked(clu)
        #print '\n\n'
    
        for cl in refine_linked(clu):
            rc = refine_cluster(cl,minratio, minreads)
            if len(rc) > 0:
                for clu in rc:
                    buf = '%s ' % chrom
                    for interval, count in clu:
                        buf += "%d:%d" % interval + ":%d"%(count)+ " "
                    Ncl += 1
                    fout.write(buf+'\n')
    sys.stderr.write("Split into %s clusters...\n"%Ncl)
    fout.close()


def merge_junctions(options):    
    ''' function to merge junctions '''

    outPrefix = options.outprefix
    rundir = options.rundir
    fnameout = "%s/%s"%(rundir,outPrefix)
    flist = "%s/%s_sortedlibs"%(rundir, outPrefix)
    
    lsts = []
    for ln in open(flist):
        lsts.append(ln.strip())
    if options.verbose:
        sys.stderr.write("merging %d junction files...\n"%(len(lsts)))
    
    # Change 300 if max open file is < 300
    N = min([300, max([100, int(len(lsts)**(0.5))])])

    tmpfiles = []
    while len(lsts) > 1:    
        clst = []
        
        for i in range(0,(len(lsts)/N)+1): 
            lst = lsts[N*i:N*(i+1)]
            if len(lst) > 0:
                clst.append(lst)
        lsts = []
    
        for lst in clst:
            if len(lst) == 0: continue
            tmpfile = tempfile.mktemp()
            os.mkdir(tmpfile)
            foutname = tmpfile+"/tmpmerge.gz"
            fout = gzip.open(foutname,'w')
            
            merge_files(lst, fout, options)
            lsts.append(foutname)
            tmpfiles.append(foutname)
            fout.close()
    
    shutil.move(lsts[0], fnameout+"_perind.counts.gz")

def merge_files(fnames, fout, options):

    fopen = []
    for fname in fnames:
        if fname[-3:] == ".gz":
            fopen.append(gzip.open(fname))
        else:
            fopen.append(open(fname))

    finished = False
    N = 0
    while not finished:
        N += 1
        if N % 50000 == 0: 
            sys.stderr.write(".")
        buf = []
        for f in fopen:
            ln = f.readline().split()
            if len(ln) == 0: 
                finished = True
                break
            chrom = ln[0]
            data = ln[1:]
            if len(buf) == 0:
                buf.append(chrom)
            buf += data
        if len(buf) > 0:
            if buf[0] == "chrom":
                if options.verbose:
                    sys.stderr.write("merging %d files"%(len(buf)-1))
            fout.write(" ".join(buf)+'\n')
        else:
            break

    sys.stderr.write(" done.\n")
    for fin in fopen:
        fin.close()


def cluster_intervals(E):
    ''' Clusters intervals together. '''
    E.sort()
    if len(E) == 0:
        return [], []
    current = E[0]
    Eclusters, cluster = [], []

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

    
    return Eclusters, E

def overlaps(A,B):
    '''
    Checks if A and B overlaps
    '''

    if A[1] < B[0] or B[1] < A[0]:
        return False
    else: return True

def refine_linked(clusters):

    unassigned = [x for x in clusters[1:]]
    current = [clusters[0]]
    splicesites = set([current[0][0][0],current[0][0][1]])
    newClusters = []
    while len(unassigned) > 0:
        finished = False
    
        while not finished:
            finished = True
            torm = []
            for intron in unassigned:
                inter, count = intron
                start, end = inter
                if start in splicesites or end in splicesites:
                    current.append(intron)
                    splicesites.add(start)
                    splicesites.add(end)
                    finished = False
                    torm.append(intron)
            for intron in torm:
                unassigned.remove(intron)
        newClusters.append(current)
        current = []
        if len(unassigned) > 0:
            current = [unassigned[0]]
            splicesites = set([current[0][0][0],current[0][0][1]])
            unassigned = unassigned[1:]
    return newClusters


def refine_cluster(clu, cutoff, readcutoff):
    ''' for each exon in the cluster compute the ratio of reads, if smaller than cutoff,
    remove and recluster '''
    
    remove = []
    dic = {}
    intervals = []

    reCLU = False
    totN = 0

    for inter, count in clu:
        totN += count
    for inter, count in clu:
        if (count/float(totN) >= cutoff and count >= readcutoff):
            intervals.append(inter)
            dic[inter] = count
        else:
            reCLU = True
    if len(intervals) == 0: return []
    
    # This makes sure that after trimming, the clusters are still good
    Atmp, B = cluster_intervals(intervals)
    A = []
    for cl in Atmp:
        for c in refine_linked([(x,0) for x in cl]):
            if len(c) > 0:
                A.append([x[0] for x in c])
    
    if len(A) == 1:
        rc = [(x, dic[x]) for x in A[0]]
        if len(rc) > 1:
            if reCLU:
                return refine_cluster([(x, dic[x]) for x in A[0]], cutoff, readcutoff)
            else:
                return [[(x, dic[x]) for x in A[0]]]
        else:
            return []
    NCs = []
    for c in A:
        if len(c) > 1:
            NC = refine_cluster([(x, dic[x]) for x in c], cutoff, readcutoff)
            NCs += NC
    return NCs


def get_numers(options):
    outPrefix = options.outprefix
    rundir = options.rundir
    fname = "%s/%s_perind.counts.gz"%(rundir,outPrefix)
    fnameout = "%s/%s_perind_numers.counts.gz"%(rundir,outPrefix)
    input_file=gzip.open(fname, 'rb')
    fout = gzip.open(fnameout,'w')
    first_line=True
    
    for l in input_file:
        if first_line:
            fout.write(" ".join(l.strip().split(" ")[1:])+'\n') # print the sample names
            first_line=False
        else:
            l=l.strip()
            words=l.split(" ")
            
            fout.write(words[0]+ " "+ " ".join( [ g.split("/")[0] for g in words[1:] ] ) +'\n')

    input_file.close()
    fout.close()

if __name__ == "__main__":

    from optparse import OptionParser

    parser = OptionParser()

    parser.add_option("-j", "--juncfiles", dest="juncfiles",
                  help="text file with all junction files to be processed")

    parser.add_option("-o", "--outprefix", dest="outprefix", default = 'leafcutter',
                  help="output prefix (default leafcutter)")

    parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")

    parser.add_option("-r", "--rundir", dest="rundir", default='./',
                  help="write to directory (default ./)")
    
    parser.add_option("-l", "--maxintronlen", dest="maxintronlen", default = 100000,
                  help="maximum intron length in bp (default 100,000bp)")

    parser.add_option("-m", "--minclureads", dest="minclureads", default = 30,
                  help="minimum reads in a cluster (default 30 reads)")

    parser.add_option("-p", "--mincluratio", dest="mincluratio", default = 0.001,
                  help="minimum fraction of reads in a cluster that support a junction (default 0.001)")
    parser.add_option("-k", "--checkchrom", dest="checkchrom", default = False,
                  help="check that the chromosomes are well formated e.g. chr1, chr2, ..., or 1, 2, ... (default=False)")

    parser.add_option("-s", "--strand", dest="strand", default = False,
                  help="use strand info (default=False)")

    (options, args) = parser.parse_args()

    if options.juncfiles == None:
        sys.stderr.write("Error: no junction file provided...\n")
        exit(0)
    
    # Get the junction file list
    libl = []
    for junc in open(options.juncfiles):
        junc = junc.strip()
        try:
            open(junc)
        except: 
            sys.stderr.write("%s does not exist... check your junction files.\n"%junc)
            exit(0)
        libl.append(junc)

    main(options, libl)
