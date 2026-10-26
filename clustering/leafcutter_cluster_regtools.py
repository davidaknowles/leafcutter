#!/usr/bin/env python
# https://github.com/griffithlab/regtools
# /home/yangili1/tools/regtools/build/regtools junctions extract -a 8 -i 50 -I 500000 bamfile.bam -o outfile.junc
# Using regtools speeds up the junction extraction step by an order of magnitude or more

import sys
import tempfile
import os
import gzip
import shutil
import array
import codecs
import collections
import io
import itertools
import multiprocessing
import operator
import queue
import threading
import traceback

# The per-sample sorted files are intermediates and are stored uncompressed inside gzip.
SORTED_LEVEL = 0
OUTPUT_LEVEL = 1

SORT_WORKERS = 4      # forked processes of sort_junctions
POOL_WORKERS = 4      # blocks of files in pool_junc_reads, one parsed by the parent
COUNT_WORKERS = 3     # forked processes of sort_and_merge
PIPE_WORKERS = 3      # forked processes of pipeline_run, beside the parent
PIPE_CHUNK = 1 << 16  # bytes a pipeline_run worker reads at a time
PIPE_MAX_CHROMS = 2000  # chromosomes beyond which pool_junc_reads' path is faster
WRITE_ROWS = 2048     # rows formatted at a time

CHROMOSOMES = ["chr%d"%x for x in range(1,23)]+['chrX','chrY']+["%d"%x for x in range(1,23)]+['X','Y']
CHROMOSOME_BYTES = frozenset(c.encode('ascii') for c in CHROMOSOMES)

# ASCII bytes besides tab and newline that str.split() splits on
SPLIT_BYTES = (b' ', b'\r', b'\x0b', b'\x0c', b'\x1c', b'\x1d', b'\x1e', b'\x1f')

def main(options,libl):

    if options.cluster is None:
        if pipeline_run(libl, options):
            return
        pooled = pool_junc_reads(libl, options)
        refine_clusters(options)
        if sort_and_merge(libl, options, pooled):
            return

    sort_junctions(libl, options)
    merge_junctions(options)
    get_numers(options)

def pool_junc_reads(flist, options):

    outPrefix = options.outprefix
    rundir = options.rundir
    maxIntronLen = int(options.maxintronlen)
    nochromcheck = options.nochromcheck

    outFile = "%s/%s_pooled"%(rundir,outPrefix)

    # check that all files exist
    for libl in flist:
        if not os.path.isfile(libl.strip()):
            raise ValueError("File %s does not exist."%lib)

    pooled = PooledReads()
    ctx = fork_context()
    if ctx is not None and len(flist) > 1:
        pool_forked(ctx, flist, options, pooled)
    else:
        for k,libl in enumerate(flist, 1):
            lib = libl.strip()
            if options.verbose:
                sys.stderr.write("scanning %d/%d: %s ...\n"%(k, len(flist), lib))
            pooled.add_file(lib, maxIntronLen, nochromcheck)
    by_chrom = pooled.by_chrom()

    fout = open(outFile, 'w')
    Ncluster = 0
    sys.stderr.write("Parsing...\n")
    for chrom in by_chrom:
        clu = pool_clusters(by_chrom[chrom])
        sys.stderr.write("%s:%s.."%chrom)
        for cl in clu:
            if len(cl) > 1: # if cluster has more than one intron
                buf = '%s:%s '%chrom
                for interval, count in [(x, by_chrom[chrom][x]) for x in cl]:
                    buf += "%d:%d" % interval + ":%d"%count+ " "
                fout.write(buf+'\n')
            Ncluster += 1
    sys.stderr.write("\nWrote %d clusters...\n"%Ncluster)
    fout.close()
    return pooled


def pool_clusters(introns):
    ''' Clusters the introns {(A, B): reads} of one chromosome and strand that have at
    least 3 reads; returns every cluster, one-intron clusters included. '''

    read_ks = [k for k,v in list(introns.items()) if v >= 3] # a junction must have at least 3 reads
    read_ks.sort()
    if len(read_ks) == 0:
        return [] # weird test case for toy data with only 1 gene - two chroms but one is empty after filtering
    return cluster_intervals(read_ks)[0]


class PooledReads(object):
    ''' Reads per intron pooled over junction files, as pool_junc_reads counts them.

    introns maps (chrom, strand, A, B) to an index into tot, the pooled reads. files has,
    per file, (block, ids, reads) for its counted lines, where ids index the introns of
    that block and idmaps[block] maps them to this object's (None: they are this
    object's). A file is None instead when sort_junctions would not count it from
    those lines alone: a line it reports, prints or fails to parse. '''

    def __init__(self):
        self.introns = {}
        self.tot = []
        self.files = []
        self.idmaps = [None]

    def add_file(self, lib, maxIntronLen, nochromcheck):
        introns, tot = self.introns, self.tot
        ids, reads = array.array('i'), array.array('i')
        clean = True

        if lib.endswith('.gz'):
            opener = gzip.open(lib, 'rt')
        else:
            opener = open(lib, 'r')

        for ln in opener:
            lnsplit = ln.split()
            if len(lnsplit) < 6:
                sys.stderr.write("Error in %s\n" % lib)
                clean = False
                continue
            chrom, A, B, dot, counts, strand, rA, rb, rgb, blockCount, blockSize, blockStarts = lnsplit
            if int(blockCount) > 2:
                print(ln)
                clean = False
                continue
            # regtools -s 0 (unstranded) now puts "?" in strand field when strand is ambiguous
            if strand == "?":
                clean = clean and sorts_cleanly(chrom, A, B, counts, blockSize, nochromcheck)
                continue
            if not nochromcheck and (chrom not in CHROMOSOMES): continue  # excludes ALT contigs etc
            Aoff, Boff = blockSize.split(",")
            A, B = int(A)+int(Aoff), int(B)-int(Boff)+1

            if B-A > int(maxIntronLen):
                clean = clean and parses_as_int(counts)
                continue
            count = int(counts)
            key = (chrom, strand, A, B)
            i = introns.get(key)
            if i is None:
                i = introns[key] = len(tot)
                tot.append(count)
            else:
                tot[i] += count
            if clean:
                try:
                    ids.append(i)
                    reads.append(count)
                except OverflowError:
                    clean = False
        opener.close()
        self.files.append((0, ids, reads) if clean else None)

    def merge(self, introns, tot, files):
        ''' Adds the reads of a later block of files. '''

        own, own_tot = self.introns, self.tot
        idmap = []
        for key, count in zip(introns, tot):
            i = own.get(key)
            if i is None:
                i = own[key] = len(own_tot)
                own_tot.append(count)
            else:
                own_tot[i] += count
            idmap.append(i)
        block = len(self.idmaps)
        self.idmaps.append(idmap)
        self.files.extend(None if f is None else (block, f[1], f[2]) for f in files)

    def by_chrom(self):
        ''' {(chrom, strand): {(A, B): reads}} in the order the files first count them. '''

        by_chrom = {}
        tot = self.tot
        for (chrom, strand, A, B), i in self.introns.items():
            introns = by_chrom.get((chrom, strand))
            if introns is None:
                introns = by_chrom[(chrom, strand)] = {}
            introns[(A, B)] = tot[i]
        return by_chrom


def sorts_cleanly(chrom, A, B, counts, blockSize, nochromcheck):
    ''' Whether sort_junctions parses a line that pool_junc_reads skips without raising. '''

    if not nochromcheck and (chrom not in CHROMOSOMES):
        return True
    try:
        Aoff, Boff = blockSize.split(",")
        int(A)+int(Aoff), int(B)-int(Boff)+1
    except ValueError:
        return False
    return parses_as_int(counts)


def parses_as_int(text):
    try:
        int(text)
    except ValueError:
        return False
    return True


def pool_block(pooled, flist, start, end, options):
    maxIntronLen = int(options.maxintronlen)
    for k in range(start, end):
        lib = flist[k].strip()
        if options.verbose:
            sys.stderr.write("scanning %d/%d: %s ...\n"%(k + 1, len(flist), lib))
        pooled.add_file(lib, maxIntronLen, options.nochromcheck)


def pool_forked(ctx, flist, options, pooled):
    ''' pool_junc_reads' parse in POOL_WORKERS contiguous blocks of files: the first in
    this process, the others in forked workers, merged and replayed in list order. '''

    blocks = contiguous_blocks([lib.strip() for lib in flist], POOL_WORKERS)
    workers = []
    try:
        for start, end in blocks[1:]:
            workers.append(fork_worker(ctx, pool_worker, [w[1] for w in workers],
                                       flist, start, end, options))
        pool_block(pooled, flist, blocks[0][0], blocks[0][1], options)
        for worker, (start, end) in zip(workers, blocks[1:]):
            log, error, tb, block = receive(worker, "scanning %s" % flist[start].strip())
            replay(log)
            if error is not None:
                raise error from WorkerTraceback(tb)
            pooled.merge(*block)
    finally:
        stop_workers(workers)


def pool_worker(conn, flist, start, end, options):
    pooled = PooledReads()
    log, error, tb = capture(pool_block, pooled, flist, start, end, options)
    if error is None:
        send_result(conn, log, None, None, (pooled.introns, pooled.tot, pooled.files))
    else:
        send_result(conn, log, error, tb)


def sort_junctions(libl, options):

    outPrefix = options.outprefix
    rundir = options.rundir

    if options.cluster is None:
        refined_cluster = "%s/%s_refined"%(rundir,outPrefix)
    else:
        refined_cluster = options.cluster
    runName = "%s/%s"%(rundir, outPrefix)


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
    for ll in libl:
        lib=ll.rstrip()
        if not os.path.isfile(lib):
            continue
        libN = lib
        if libN not in merges:
            merges[libN] = []
        merges[libN].append(lib)

    jobs = []
    for libN in merges:
        libName = "%s/%s"%(rundir,libN.split('/')[-1])
        jobs.append((libN, merges[libN], libName+'.%s.sorted.gz'%(runName.split("/")[-1])))

    fout_runlibs = open(runName+"_sortedlibs",'w')
    ctx = fork_context()
    if ctx is not None and len(jobs) > 1:
        sort_forked(ctx, jobs, cluExons, options, fout_runlibs)
    else:
        for k, (libN, libs, foutName) in enumerate(jobs, 1):
            fout_runlibs.write(foutName+'\n')
            sorting_messages(k, len(jobs), libN, libs, options)
            sort_sample(libN, libs, foutName, cluExons, options)
    fout_runlibs.close()


def sorting_messages(k, n, libN, libs, options):
    if options.verbose:
        sys.stderr.write("Sorting %d/%d: %s ...\n"%(k, n, libN))
    if len(libs) > 1:
        if options.verbose:
            sys.stderr.write("merging %s...\n"%(" ".join(libs)))


def sort_sample(libN, libs, foutName, cluExons, options):
    ''' Writes to foutName one sample's reads for every intron of cluExons. '''

    nochromcheck = options.nochromcheck
    by_chrom = {}

    with gzip.open(foutName, 'wb', compresslevel=SORTED_LEVEL) as fout:
        header_string = "chrom %s\n"%libN.split("/")[-1].split(".junc")[0]
        fout.write(header_string.encode('utf-8'))

        for lib in libs:
            if lib.endswith('.gz'):
                opener = gzip.open(lib, 'rt')
            else:
                opener = open(lib, 'r')

            for ln in opener:
                lnsplit = ln.split()
                if len(lnsplit) < 6:
                    sys.stderr.write("Error in %s \n" % lib)
                    continue
                chrom, A, B, dot, count, strand, rA, rb, rgb, blockCount, blockSize, blockStarts = lnsplit
                if int(blockCount) > 2:
                    print(ln)
                    continue

                if not nochromcheck and (chrom not in CHROMOSOMES): continue
                Aoff, Boff = blockSize.split(",")
                A, B = int(A)+int(Aoff), int(B)-int(Boff)+1
                chrom = (chrom,strand)
                if chrom not in by_chrom:
                    by_chrom[chrom] = {}
                intron = (A, B)

                if intron in by_chrom[chrom]:
                    by_chrom[chrom][intron] += int(count)
                else:
                    by_chrom[chrom][intron] = int(count)
            opener.close()

        for clu in cluExons:
            buf = []
            ks = cluExons[clu]
            ks.sort()
            tot = 0
            for exon in ks:
                chrom, start, end = exon
                chrom = tuple(chrom.split(":"))
                start, end = int(start), int(end)
                if chrom not in by_chrom:
                    pass
                elif (start,end) in by_chrom[chrom]:
                    tot += by_chrom[chrom][(start,end)]
            for exon in ks:

                chrom, start, end = exon
                start, end = int(start), int(end)
                chrom = tuple(chrom.split(":"))
                chromID, strand = chrom
                if chrom not in by_chrom:
                    buf.append("%s:%d:%d:clu_%d_%s 0/%d\n"%(chromID,start, end,clu, strand, tot))
                elif (start,end) in by_chrom[chrom]:
                    buf.append("%s:%d:%d:clu_%d_%s %d/%d\n"%(chromID,start, end, clu,strand, by_chrom[chrom][(start,end)], tot))
                else:
                    buf.append("%s:%d:%d:clu_%d_%s 0/%d\n"%(chromID,start, end,clu,strand, tot))

            fout.write("".join(buf).encode('utf-8'))


def sort_forked(ctx, jobs, cluExons, options, fout_runlibs):
    ''' sort_sample for each job in SORT_WORKERS forked workers, each writing to a .part
    file; _sortedlibs, messages, renames and errors follow the serial order. '''

    parts = ["%s.%d.%d.part" % (job[2], os.getpid(), k) for k, job in enumerate(jobs)]
    workers = []
    try:
        for w in range(min(SORT_WORKERS, len(jobs))):
            mine = [(jobs[k], parts[k]) for k in range(w, len(jobs), SORT_WORKERS)]
            workers.append(fork_worker(ctx, sort_worker, [x[1] for x in workers],
                                       mine, cluExons, options))
        for k, (libN, libs, foutName) in enumerate(jobs):
            fout_runlibs.write(foutName+'\n')
            sorting_messages(k + 1, len(jobs), libN, libs, options)
            log, error, tb, _ = receive(workers[k % SORT_WORKERS], "sorting %s" % libN)
            replay(log)
            if os.path.exists(parts[k]):
                os.replace(parts[k], foutName)
            if error is not None:
                raise error from WorkerTraceback(tb)
    finally:
        stop_workers(workers)
        remove_files(parts)


def sort_worker(conn, mine, cluExons, options):
    for (libN, libs, foutName), part in mine:
        log, error, tb = capture(sort_sample, libN, libs, part, cluExons, options)
        send_result(conn, log, error, tb)
        if error is not None:
            return


def refine_clusters(options):

    outPrefix = options.outprefix
    rundir = options.rundir
    minratio = float(options.mincluratio)
    minreads = int(options.minclureads)

    inFile = "%s/%s_pooled"%(rundir,outPrefix)
    outFile = "%s/%s_refined"%(rundir,outPrefix)
    fout = open(outFile,'w')

    Ncl = 0
    for ln in open(inFile):
        clu = []
        chrom = ln.split()[0]
        for ex in ln.split()[1:]:
            A, B, N = ex.split(":")
            clu.append(((int(A),int(B)), int(N)))

        for cl in refine_pooled(clu, minratio, minreads, options.const):
            buf = '%s ' % chrom
            for interval, count in cl:
                buf += "%d:%d" % interval + ":%d"%(count)+ " "
            Ncl += 1
            fout.write(buf+'\n')

    sys.stderr.write("Split into %s clusters...\n"%Ncl)
    fout.close()


def refine_pooled(clu, minratio, minreads, const):
    ''' The clusters refine_clusters writes for one _pooled line, clu = [((A, B), reads)]. '''

    totN = 0
    for interval, count in clu:
        totN += count
    if totN < minreads: return []

    refined = []
    if const:
        if len(clu) == 1:
            refined.append(clu)

    for cl in refine_linked(clu):
        rc = refine_cluster(cl,minratio, minreads)
        refined.extend(rc)
    return refined


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

        for i in range(0,(len(lsts)//N)+1):
            lst = lsts[N*i:N*(i+1)]
            if len(lst) > 0:
                clst.append(lst)
        lsts = []

        for lst in clst:
            if len(lst) == 0: continue
            tmpfile = tempfile.mktemp()
            os.mkdir(tmpfile)
            foutname = tmpfile+"/tmpmerge.gz"
            fout = gzip.open(foutname,'w',compresslevel=OUTPUT_LEVEL)

            merge_files(lst, fout, options)
            lsts.append(foutname)
            tmpfiles.append(foutname)
            fout.close()

    if not options.const:
        shutil.move(lsts[0], fnameout+"_perind.counts.gz")
    else:
        shutil.move(lsts[0], fnameout+"_perind.constcounts.gz")


def merge_files(fnames, fout, options):

    fopen = []
    for fname in fnames:
        if fname.endswith(".gz"):
            fopen.append(gzip.open(fname, "rt"))
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

            out_string = " ".join(buf)+'\n'
            fout.write(out_string.encode('utf-8'))
        else:
            break

    sys.stderr.write(" done.\n")
    for fin in fopen:
        fin.close()


def merge_messages(options, tokens, nrows):
    ''' Prints what merge_junctions prints when it merges sorted files whose headers hold
    tokens[k] after "chrom" and nrows rows each. '''

    if options.verbose:
        sys.stderr.write("merging %d junction files...\n"%(len(tokens)))
    widths = [len(t) for t in tokens]
    N = min([300, max([100, int(len(widths)**(0.5))])])
    while len(widths) > 1:
        groups = [widths[i:i + N] for i in range(0, len(widths), N)]
        for group in groups:
            if options.verbose:
                sys.stderr.write("merging %d files"%(sum(group)))
            # merge_files prints a dot every 50,000 lines, counting the header and the end
            sys.stderr.write("." * ((nrows + 2) // 50000))
            sys.stderr.write(" done.\n")
        widths = [sum(group) for group in groups]


def cluster_intervals(E):
    ''' Clusters intervals together. '''
    E.sort()
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

    if not options.const:
        fname = "%s/%s_perind.counts.gz"%(rundir,outPrefix)
        fnameout = "%s/%s_perind_numers.counts.gz"%(rundir,outPrefix)
    else:
        fname = "%s/%s_perind.constcounts.gz"%(rundir,outPrefix)
        fnameout = "%s/%s_perind_numers.constcounts.gz"%(rundir,outPrefix)

    input_file=gzip.open(fname, 'rt')
    fout = gzip.open(fnameout,'w',compresslevel=OUTPUT_LEVEL)
    first_line=True

    for l in input_file:
        if first_line:
            header_string = " ".join(l.strip().split(" ")[1:])+'\n'
            fout.write(header_string.encode('utf-8')) # print the sample names
            first_line=False
        else:
            l=l.strip()
            words=l.split(" ")
            words_string = words[0]+ " "+ " ".join( [ g.split("/")[0] for g in words[1:] ] ) +'\n'
            fout.write(words_string.encode('utf-8'))

    input_file.close()
    fout.close()


# Forked workers. Only the 'fork' start method is used: the others pickle the target,
# and the tests load this script without registering it as a module. For the same
# reason, workers send only built-in types.

def fork_context():
    ''' The fork start method, or None where forking is unavailable or pointless. '''

    if not hasattr(os, 'fork') or (os.cpu_count() or 1) < 2:
        return None
    return multiprocessing.get_context('fork')


def fork_worker(ctx, target, inherited, *args):
    ''' Runs target(conn, *args) in a forked process; returns (process, connection that
    receives what target sends). The child closes the inherited connections. '''

    receiver, sender = ctx.Pipe(duplex=False)
    process = ctx.Process(target=worker_main,
                          args=(target, sender, [receiver] + list(inherited), args))
    process.daemon = True
    process.start()
    sender.close()
    return process, receiver


def worker_main(target, conn, inherited, args):
    for other in inherited:
        other.close()
    # a worker must not write to the parent's streams; what it prints is sent back
    sys.stdout = sys.stderr = OutputLog([], 'stderr')
    target(conn, *args)


def stop_workers(workers):
    for process, conn in workers:
        if process.is_alive():
            process.terminate()
        process.join()
        conn.close()


def receive(worker, what):
    try:
        return worker[1].recv()
    except EOFError:
        raise RuntimeError("the worker %s exited unexpectedly" % what)


class OutputLog(object):
    ''' Stands in for sys.stdout or sys.stderr, recording what is written to it. '''

    def __init__(self, log, stream):
        self.log = log
        self.stream = stream

    def write(self, text):
        self.log.append((self.stream, text))
        return len(text)

    def flush(self):
        pass


class WorkerTraceback(Exception):
    ''' The traceback of an exception raised in a forked worker. '''

    def __str__(self):
        return self.args[0]


def capture(func, *args):
    ''' Calls func with its output recorded; returns (log, exception, traceback text). '''

    log = []
    stdout, stderr = sys.stdout, sys.stderr
    sys.stdout, sys.stderr = OutputLog(log, 'stdout'), OutputLog(log, 'stderr')
    try:
        func(*args)
        return log, None, None
    except Exception as error:
        return log, error, traceback.format_exc()
    finally:
        sys.stdout, sys.stderr = stdout, stderr


def send_result(conn, log, error, tb, result=None):
    try:
        conn.send((log, error, tb, result))
    except Exception:
        if error is None:
            raise
        conn.send((log, RuntimeError("%s: %s" % (type(error).__name__, error)), tb, None))


def replay(log):
    for stream, text in log:
        getattr(sys, stream).write(text)


def remove_files(paths):
    for path in paths:
        try:
            os.remove(path)
        except OSError:
            pass


def contiguous_blocks(paths, n):
    ''' Cuts paths into at most n contiguous (start, end) blocks of about equal bytes. '''

    sizes = [os.path.getsize(path) for path in paths]
    total = float(sum(sizes))
    blocks, start, done = [], 0, 0
    for i, size in enumerate(sizes[:-1]):
        done += size
        if len(blocks) < n - 1 and done >= total * (len(blocks) + 1) / n:
            blocks.append((start, i + 1))
            start = i + 1
    blocks.append((start, len(paths)))
    return blocks


# Writing the counts of a pooled run

def sample_names(libs, options):
    ''' The sorted-file names and header names sort_junctions gives libs, or
    (None, None) when a name could be read back differently or two names collide on a
    case-insensitive file system. '''

    runName = "%s/%s"%(options.rundir, options.outprefix)
    names = ["%s/%s.%s.sorted.gz"%(options.rundir, lib.split('/')[-1], runName.split("/")[-1])
             for lib in libs]
    headers = [lib.split("/")[-1].split(".junc")[0] for lib in libs]
    for name in names + headers:
        if not name.isascii() or name != name.strip() or '\n' in name or '\r' in name:
            return None, None
    if len(set(name.lower() for name in names)) != len(names):
        return None, None
    return names, headers


class GzipWriter(threading.Thread):
    ''' Writes the texts given to put() to gzip files from its own thread. '''

    def __init__(self, handles):
        threading.Thread.__init__(self)
        self.daemon = True
        self.handles = handles
        self.queue = queue.Queue(4)
        self.error = None

    def run(self):
        while True:
            texts = self.queue.get()
            if texts is None:
                return
            if self.error is None:
                try:
                    for i, text in texts:
                        self.handles[i].write(text.encode('utf-8'))
                except Exception as error:
                    self.error = error

    def put(self, texts):
        if self.error is not None:
            raise self.error
        self.queue.put(texts)

    def finish(self):
        self.queue.put(None)
        self.join()
        if self.error is not None:
            raise self.error


class CountsWriter(object):
    ''' Writes a pooled run's sorted files and merged counts under temporary names.

    The sorted files and the numerators go to .part files next to their names, the
    perind file to tmpmerge.gz in a new temp directory, as merge_junctions writes it.
    install() puts them in place in the order sort_junctions, merge_junctions and
    get_numers write them. '''

    def __init__(self, options, libs, names, headers):
        self.options = options
        self.libs = libs
        self.names = names
        self.headers = headers
        self.tokens = [("chrom %s\n"%header).split()[1:] for header in headers]
        self.parts = ["%s.%d.%d.part" % (name, os.getpid(), k) for k, name in enumerate(names)]
        ext = "constcounts" if options.const else "counts"
        self.perind = "%s/%s_perind.%s.gz"%(options.rundir, options.outprefix, ext)
        self.numers = "%s/%s_perind_numers.%s.gz"%(options.rundir, options.outprefix, ext)
        self.numers_part = "%s.%d.part" % (self.numers, os.getpid())
        self.tmpdir = None
        self.handles = []
        self.writers = []
        self.nrows = 0

    def open(self):
        ''' Opens every file and starts the writer threads (after any fork). '''

        for part, header in zip(self.parts, self.headers):
            handle = gzip.open(part, 'wb', compresslevel=SORTED_LEVEL)
            self.handles.append(handle)
            handle.write(("chrom %s\n"%header).encode('utf-8'))
        sorted_handles = list(self.handles)

        tmpfile = tempfile.mktemp()
        os.mkdir(tmpfile)
        self.tmpdir = tmpfile
        perind = gzip.open(tmpfile+"/tmpmerge.gz", 'wb', compresslevel=OUTPUT_LEVEL)
        self.handles.append(perind)
        tokens = [t for header in self.tokens for t in header]
        perind.write(("chrom" + "".join(" " + t for t in tokens) + "\n").encode('utf-8'))
        numers = gzip.open(self.numers_part, 'wb', compresslevel=OUTPUT_LEVEL)
        self.handles.append(numers)
        numers.write((" ".join(tokens) + "\n").encode('utf-8'))

        for handles in (sorted_handles, [perind], [numers]):
            writer = GzipWriter(handles)
            writer.start()
            self.writers.append(writer)

    def write(self, names, counts, totals):
        ''' Writes rows: names[i] has counts[j][i] reads for sample j, whose cluster has
        totals[j][i]. '''

        nsamples = len(self.names)
        sorted_writer, perind_writer, numers_writer = self.writers
        for start in range(0, len(names), WRITE_ROWS):
            end = start + WRITE_ROWS
            block = names[start:end]
            n = len(block)
            columns = [(counts[j][start:end], totals[j][start:end]) for j in range(nsamples)]

            fmt = "%s %d/%d\n" * n
            args = [None] * (3 * n)
            args[0::3] = block
            texts = []
            for j, (c, t) in enumerate(columns):
                args[1::3] = c
                args[2::3] = t
                texts.append((j, fmt % tuple(args)))
            sorted_writer.put(texts)

            width = 1 + 2 * nsamples
            args = [None] * (width * n)
            args[0::width] = block
            for j, (c, t) in enumerate(columns):
                args[1 + 2 * j::width] = c
                args[2 + 2 * j::width] = t
            perind_writer.put([(0, (("%s" + " %d/%d" * nsamples + "\n") * n) % tuple(args))])

            width = 1 + nsamples
            args = [None] * (width * n)
            args[0::width] = block
            for j, (c, t) in enumerate(columns):
                args[1 + j::width] = c
            numers_writer.put([(0, (("%s" + " %d" * nsamples + "\n") * n) % tuple(args))])
            self.nrows += n

    def close(self):
        for writer in self.writers:
            writer.finish()
        self.writers = []
        for handle in self.handles:
            handle.close()

    def discard(self):
        for writer in self.writers:
            writer.queue.put(None)
        for writer in self.writers:
            writer.join()
        for handle in self.handles:
            try:
                handle.close()
            except Exception:
                pass
        remove_files(self.parts + [self.numers_part])
        if self.tmpdir is not None:
            shutil.rmtree(self.tmpdir, ignore_errors=True)

    def install(self):
        ''' What sort_junctions, merge_junctions and get_numers would write and print. '''

        options = self.options
        runName = "%s/%s"%(options.rundir, options.outprefix)
        done = 0
        try:
            fout_runlibs = open(runName+"_sortedlibs",'w')
            for k, lib in enumerate(self.libs):
                fout_runlibs.write(self.names[k]+'\n')
                sorting_messages(k + 1, len(self.libs), lib, [lib], options)
                os.replace(self.parts[k], self.names[k])
                done += 1
            fout_runlibs.close()
            merge_messages(options, self.tokens, self.nrows)
            landed = shutil.move(self.tmpdir+"/tmpmerge.gz", self.perind)
        except BaseException:
            remove_files(self.parts[done:] + [self.numers_part])
            raise
        if landed != self.perind:
            remove_files([self.numers_part])
            get_numers(options)
        try:
            os.replace(self.numers_part, self.numers)
        except OSError:
            remove_files([self.numers_part])
            get_numers(options)


def int_array(values):
    try:
        return array.array('q', values)
    except OverflowError:
        return values


def file_counts(ids, reads, slots, sizes, ends):
    ''' Reads per row from one file's counted lines and, per row, the total of its
    cluster. slots[id] is the row of intron id, len(rows) if none; sizes and ends are
    the clusters' row counts and cumulative row counts. '''

    nrows = ends[-1] if ends else 0
    acc = [0] * (nrows + 1)
    for i, r in zip(ids, reads):
        acc[slots[i]] += r
    del acc[nrows]
    prefix = [0]
    prefix.extend(itertools.accumulate(acc))
    cluster_totals = map(operator.sub, map(prefix.__getitem__, ends),
                         map(prefix.__getitem__, [0] + ends[:-1]))
    totals = itertools.chain.from_iterable(map(itertools.repeat, cluster_totals, sizes))
    return array.array('q', acc), array.array('q', totals)


def read_refined(path):
    ''' The clusters of a _refined file as (chrom, strand, introns in row order), or None
    when sort_junctions would fail on it or would not count it like pooling does. '''

    clusters = []
    try:
        with open(path) as handle:
            for ln in handle:
                fields = ln.split()
                key = fields[0].split(":")
                if len(key) != 2 or not fields[0].isascii():
                    return None
                ks = []
                for exon in fields[1:]:
                    A, B, count = exon.split(":")
                    int(A), int(B)
                    ks.append((A, B))
                ks.sort()
                clusters.append((key[0], key[1], [(int(A), int(B)) for A, B in ks]))
    except (OSError, ValueError, IndexError):
        return None
    return clusters


def sort_and_merge(libl, options, pooled):
    ''' sort_junctions, merge_junctions and get_numers for a pooled run, from the lines
    pool_junc_reads counted: forked workers sum each sample's reads per row and this
    process writes the outputs. Returns False, having written and printed nothing, when
    the serial code must run instead. '''

    ctx = fork_context()
    if ctx is None or len(libl) < 2 or len(pooled.files) != len(libl) or None in pooled.files:
        return False
    names, headers = sample_names(libl, options)
    clusters = read_refined("%s/%s_refined"%(options.rundir, options.outprefix))
    if names is None or clusters is None:
        return False

    rows, ids, sizes = [], [], []
    for cluN, (chromID, strand, introns) in enumerate(clusters, 1):
        sizes.append(len(introns))
        for A, B in introns:
            rows.append("%s:%d:%d:clu_%d_%s"%(chromID, A, B, cluN, strand))
            ids.append(pooled.introns.get((chromID, strand, A, B)))
    slots = [len(rows)] * len(pooled.tot)
    for row, i in enumerate(ids):
        if i is not None:
            if slots[i] != len(rows):
                return False
            slots[i] = row
    ends = list(itertools.accumulate(sizes))

    writer = CountsWriter(options, libl, names, headers)
    workers = []
    try:
        for w in range(COUNT_WORKERS):
            workers.append(fork_worker(ctx, count_worker, [x[1] for x in workers],
                                       pooled, slots, sizes, ends, w))
        writer.open()
        columns = [None] * len(libl)
        for worker in workers:
            for k, counts, totals in receive(worker, "counting"):
                columns[k] = (counts, totals)
        writer.write(rows, [c[0] for c in columns], [c[1] for c in columns])
        writer.close()
    except Exception:
        writer.discard()
        return False
    except BaseException:
        writer.discard()
        raise
    finally:
        stop_workers(workers)
    writer.install()
    return True


def count_worker(conn, pooled, slots, sizes, ends, first):
    block_slots = {}
    out = []
    for k in range(first, len(pooled.files), COUNT_WORKERS):
        block, ids, reads = pooled.files[k]
        if block not in block_slots:
            idmap = pooled.idmaps[block]
            block_slots[block] = slots if idmap is None else [slots[i] for i in idmap]
        out.append((k,) + file_counts(ids, reads, block_slots[block], sizes, ends))
    conn.send(out)


# pipeline_run: pooling, clustering and counting a chromosome at a time

class PipelineFallback(Exception):
    ''' The input needs the serial code. '''


def pipeline_run(libl, options):
    ''' Clusters a run without -c chromosome by chromosome: forked workers parse
    contiguous blocks of the junction files while this process pools, refines and
    writes each chromosome once every file has passed it. Returns False, having written
    and printed nothing, when the serial code must run instead: no fork, options or file
    names it would read differently, a line it reports, prints or fails on, files whose
    lines of a chromosome are not together or whose chromosome orders disagree, more
    than PIPE_MAX_CHROMS chromosomes, or any error before the outputs are installed. '''

    ctx = fork_context()
    if ctx is None or not pipeline_accepts(libl, options):
        return False
    names, headers = sample_names(libl, options)
    if names is None:
        return False
    run = Pipeline(libl, options, names, headers)
    try:
        run.start(ctx)
        run.collect()
    except Exception:
        run.discard()
        return False
    except BaseException:
        run.discard()
        raise
    run.install()
    return True


def pipeline_accepts(libl, options):
    try:
        int(options.maxintronlen), float(options.mincluratio), int(options.minclureads)
    except Exception:
        return False
    if len(libl) < 2 or len(set(libl)) != len(libl):
        return False
    for lib in libl:
        if lib != lib.strip() or not os.path.isfile(lib):
            return False
    # the workers decode lines as UTF-8, as text mode does with this encoding only
    return codecs.lookup(io.TextIOWrapper(io.BytesIO()).encoding).name == 'utf-8'


class FileOrder(object):
    ''' The sequence of chromosomes in each file, and where each chromosome comes. '''

    def __init__(self):
        self.where = {}
        self.seqs = {}

    def visit(self, f, chrom):
        seq = self.seqs.setdefault(f, [])
        self.where.setdefault(chrom, {})[f] = len(seq)
        seq.append(chrom)

    def before(self, a, b):
        ''' Whether some file holds a before b. '''

        wa, wb = self.where.get(a), self.where.get(b)
        if wa is None or wb is None:
            return False
        if len(wa) <= len(wb):
            return any(pa < wb.get(f, -1) for f, pa in wa.items())
        return any(wa.get(f, pb) < pb for f, pb in wb.items())


class PipeChrom(object):
    ''' A chromosome in the parent of pipeline_run. '''

    def __init__(self, name):
        self.name = name
        self.blocks = set()
        self.requested = False
        self.sums = {}
        self.counts = {}
        self.keys = []
        self.columns = None


class PipeKey(object):
    ''' The clusters of one (chrom, strand) and where it comes in the serial order:
    first is the first file that counts it, order sorts the keys that file meets first.

    pooled and refined hold its lines of _pooled and _refined, rows its refined
    clusters' introns in row order, start and end its rows within its chromosome. '''

    def __init__(self, chrom, strand, first, order):
        self.chrom = chrom
        self.strand = strand
        self.first = first
        self.order = order
        self.ngroups = 0
        self.pooled = self.refined = ''
        self.rows = []
        self.start = self.end = 0
        self.counts = self.totals = None

    def cluster(self, introns, minratio, minreads, const):
        ''' What pool_junc_reads and refine_clusters make of introns {(A, B): reads}. '''

        clusters = pool_clusters(introns)
        self.ngroups = len(clusters)
        pooled, refined = [], []
        for cl in clusters:
            if len(cl) > 1:
                clu = [(x, introns[x]) for x in cl]
                pooled.append(self.line(clu))
                refined.extend(refine_pooled(clu, minratio, minreads, const))
        self.pooled = "".join(pooled)
        self.refined = "".join(self.line(clu) for clu in refined)
        self.rows = [sorted((x for x, count in clu), key=lambda x: ("%d"%x[0], "%d"%x[1]))
                     for clu in refined]

    def line(self, clu):
        return "%s:%s "%(self.chrom, self.strand) + "".join(
            "%d:%d:%d "%(A, B, count) for (A, B), count in clu) + "\n"


class Pipeline(object):
    ''' The parent of pipeline_run's workers.

    Workers report where each of their files is (its head: the chromosome of its next
    counted line). Once no file can hold more lines of a chromosome, the parent asks the
    workers for its pooled reads, clusters it, and asks for each file's reads per row.
    Rows are written in the serial order: clusters of (chrom, strand) keys the first file
    counts as that file meets them, then the other keys in the order the later files
    first count them. '''

    def __init__(self, libs, options, names, headers):
        self.libs = libs
        self.options = options
        self.minratio = float(options.mincluratio)
        self.minreads = int(options.minclureads)
        self.blocks = contiguous_blocks(libs, PIPE_WORKERS)
        self.block_of = [k for k, (start, end) in enumerate(self.blocks) for _ in range(start, end)]
        self.heads = [None] * len(libs)
        self.head_files = {}
        self.order = FileOrder()
        self.chroms = {}
        self.waiting = collections.deque()
        self.reported = 0
        self.parsed = 0
        self.unfinished = 0
        self.written = 0
        self.keys = []
        self.late = []
        self.ncluster = 0
        self.events = queue.Queue()
        self.processes = []
        self.outboxes = []
        self.threads = []
        self.writer = CountsWriter(options, libs, names, headers)

    def start(self, ctx):
        pipes = [(ctx.Pipe(duplex=False), ctx.Pipe(duplex=False)) for _ in self.blocks]
        ends = [end for pair in pipes for pipe in pair for end in pipe]
        for (start, end), ((_, from_worker), (to_worker, _)) in zip(self.blocks, pipes):
            others = [x for x in ends if x is not from_worker and x is not to_worker]
            process = ctx.Process(target=worker_main, args=(
                pipe_worker, from_worker, others, (to_worker, self.libs[start:end], self.options)))
            process.daemon = True
            process.start()
            self.processes.append(process)
        for (_, from_worker), (to_worker, _) in pipes:
            from_worker.close()
            to_worker.close()
        self.writer.open()
        for k, ((receiver, _), (_, sender)) in enumerate(pipes):
            outbox = queue.Queue()
            self.outboxes.append(outbox)
            for target, args in ((pipe_reader, (k, receiver, self.events)),
                                 (pipe_sender, (sender, outbox))):
                thread = threading.Thread(target=target, args=args)
                thread.daemon = True
                thread.start()
                self.threads.append(thread)

    def collect(self):
        while self.parsed < len(self.blocks) or self.unfinished:
            k, message = self.events.get()
            getattr(self, 'on_' + message[0])(k, *message[1:])
        self.late.sort(key=lambda key: key.order)
        for key in self.late:
            self.write_key(key, key.counts, key.totals)
        for outbox in self.outboxes:
            outbox.put(('finish',))
            outbox.put(None)
        for process in self.processes:
            process.join()
        for thread in self.threads:
            thread.join()
        self.writer.close()

    def discard(self):
        for process in self.processes:
            if process.is_alive():
                process.terminate()
        for process in self.processes:
            process.join()
        for outbox in self.outboxes:
            outbox.put(None)
        for thread in self.threads:
            thread.join()
        self.writer.discard()

    def on_heads(self, k, heads):
        start = self.blocks[k][0]
        for i, head in enumerate(heads):
            self.move(start + i, head)
        self.reported += 1
        self.finalize()

    def on_step(self, k, chrom, moved):
        start = self.blocks[k][0]
        for i, head in moved:
            self.move(start + i, head)
        self.finalize()

    def on_parsed(self, k):
        self.parsed += 1

    def on_done(self, k):
        pass

    def on_error(self, k, text):
        raise PipelineFallback(text)

    def on_died(self, k):
        raise PipelineFallback("a worker exited unexpectedly")

    def move(self, g, head):
        old = self.heads[g]
        if old is not None:
            self.head_files[old] -= 1
            if not self.head_files[old]:
                del self.head_files[old]
        self.heads[g] = head
        if head is None:
            return
        chrom = self.chroms.get(head)
        if chrom is None:
            if len(self.chroms) == PIPE_MAX_CHROMS:
                raise PipelineFallback("more than %d chromosomes" % PIPE_MAX_CHROMS)
            chrom = self.chroms[head] = PipeChrom(head)
            self.waiting.append(chrom)
            self.unfinished += 1
        elif chrom.requested:
            raise PipelineFallback("chromosome orders disagree")
        chrom.blocks.add(self.block_of[g])
        self.order.visit(g, head)
        self.head_files[head] = self.head_files.get(head, 0) + 1

    def finalize(self):
        ''' Asks for the reads of the chromosomes no file can still hold, in the order
        the files first reach them, up to the first one a file may still hold. '''

        if self.reported < len(self.blocks):
            return
        waiting = self.waiting
        while waiting:
            chrom = waiting[0]
            if not all(head != chrom.name and self.order.before(chrom.name, head)
                       for head in self.head_files):
                return
            waiting.popleft()
            chrom.requested = True
            for k in chrom.blocks:
                self.outboxes[k].put(('sums', chrom.name))

    def on_sums(self, k, name, strands, firsts, groups):
        chrom = self.chroms[name]
        chrom.sums[k] = (strands, firsts, groups)
        if len(chrom.sums) == len(chrom.blocks):
            self.pipe_chrom(chrom)

    def pipe_chrom(self, chrom):
        ''' Pools and refines a chromosome's keys and asks for each file's reads per row. '''

        name = chrom.name.decode('ascii')
        strands, firsts = [], {}
        by_strand = {}
        block_keys = []
        for k in sorted(chrom.sums):
            bstrands, bfirsts, groups = chrom.sums[k]
            keys = []
            for s, f, (A, B, tot) in zip(bstrands, bfirsts, groups):
                if s not in firsts:
                    firsts[s] = self.blocks[k][0] + f
                    strands.append(s)
                keys.append(list(zip(A, B)))
                introns = by_strand.get(s)
                if introns is None:
                    by_strand[s] = dict(zip(keys[-1], tot))
                else:
                    get = introns.get
                    for key, t in zip(keys[-1], tot):
                        introns[key] = get(key, 0) + t
            block_keys.append((k, bstrands, keys))
        chrom.sums = None

        rows = dict((s, {}) for s in strands)
        nrows = 0
        sizes = []
        where = self.order.where[chrom.name]
        for index, s in enumerate(strands):
            first = firsts[s]
            key = PipeKey(name, s.decode('ascii'), first, (first, where[first], index))
            key.cluster(by_strand[s], self.minratio, self.minreads, self.options.const)
            if key.rows and (':' in key.chrom or ':' in key.strand):
                raise PipelineFallback("a cluster name with more than one ':'")
            key.start = nrows
            for introns in key.rows:
                sizes.append(len(introns))
                for intron in introns:
                    rows[s][intron] = nrows
                    nrows += 1
            key.end = nrows
            chrom.keys.append(key)
        chrom.nrows = nrows
        ends = list(itertools.accumulate(sizes))
        for k, bstrands, keys in block_keys:
            slots = [list(map(rows[s].get, ks, itertools.repeat(nrows)))
                     for s, ks in zip(bstrands, keys)]
            self.outboxes[k].put(('map', chrom.name, slots, sizes, ends))

    def on_counts(self, k, name, files):
        chrom = self.chroms[name]
        chrom.counts[k] = files
        if len(chrom.counts) < len(chrom.blocks):
            return
        zeros = array.array('q', bytes(8 * chrom.nrows))
        columns = [(zeros, zeros)] * len(self.libs)
        for block, files in chrom.counts.items():
            start = self.blocks[block][0]
            for i, data in enumerate(files):
                if data is not None:
                    columns[start + i] = data
        chrom.columns = columns
        chrom.counts = None
        self.unfinished -= 1
        if 0 not in self.order.where[chrom.name]:
            self.pipe_write(chrom)
        seq = self.order.seqs.get(0, [])
        while self.written < len(seq) and self.chroms[seq[self.written]].columns is not None:
            self.pipe_write(self.chroms[seq[self.written]])
            self.written += 1

    def pipe_write(self, chrom):
        ''' Writes the keys the first file counts and keeps the others for the end. '''

        for key in chrom.keys:
            counts = [c[key.start:key.end] for c, t in chrom.columns]
            totals = [t[key.start:key.end] for c, t in chrom.columns]
            if key.first == 0:
                self.write_key(key, counts, totals)
            else:
                key.counts, key.totals = counts, totals
                self.late.append(key)
        chrom.columns = ()

    def write_key(self, key, counts, totals):
        names = []
        for introns in key.rows:
            self.ncluster += 1
            for A, B in introns:
                names.append("%s:%d:%d:clu_%d_%s"%(key.chrom, A, B, self.ncluster, key.strand))
        self.writer.write(names, counts, totals)
        key.counts = key.totals = None
        self.keys.append(key)

    def install(self):
        ''' Writes what pool_junc_reads and refine_clusters write, then installs. '''

        options = self.options
        outPrefix, rundir = options.outprefix, options.rundir
        try:
            if options.verbose:
                for k, lib in enumerate(self.libs, 1):
                    sys.stderr.write("scanning %d/%d: %s ...\n"%(k, len(self.libs), lib))
            fout = open("%s/%s_pooled"%(rundir,outPrefix), 'w')
            Ncluster = 0
            sys.stderr.write("Parsing...\n")
            for key in self.keys:
                sys.stderr.write("%s:%s.."%(key.chrom, key.strand))
                fout.write(key.pooled)
                Ncluster += key.ngroups
            sys.stderr.write("\nWrote %d clusters...\n"%Ncluster)
            fout.close()

            fout = open("%s/%s_refined"%(rundir,outPrefix),'w')
            Ncl = 0
            for key in self.keys:
                fout.write(key.refined)
                Ncl += len(key.rows)
            sys.stderr.write("Split into %s clusters...\n"%Ncl)
            fout.close()
        except BaseException:
            self.writer.discard()
            raise
        self.writer.install()


def pipe_reader(k, conn, events):
    try:
        while True:
            message = conn.recv()
            events.put((k, message))
            if message[0] in ('done', 'error'):
                return
    except Exception:
        events.put((k, ('died',)))
    finally:
        conn.close()


def pipe_sender(conn, outbox):
    try:
        while True:
            message = outbox.get()
            if message is None:
                return
            conn.send(message)
    except Exception:
        pass
    finally:
        conn.close()


def pipe_worker(out_conn, in_conn, libs, options):
    # a thread sends, so that parsing goes on while the parent is slow to read
    outbox = queue.Queue()
    sender = threading.Thread(target=pipe_sender, args=(out_conn, outbox))
    sender.daemon = True
    sender.start()
    try:
        PipeWorker(outbox, in_conn, libs, options).run()
    except Exception as error:
        outbox.put(('error', '%s: %s' % (type(error).__name__, error)))
    outbox.put(None)
    sender.join()


class FileLines(object):
    ''' One file's counted lines of a chromosome: intron ids, reads and strands met. '''

    def __init__(self):
        self.ids = array.array('i')
        self.reads = array.array('i')
        self.strands = []


class ChromCounts(object):
    ''' A chromosome's counted lines in a worker's block of files.

    introns maps (strand, A, B) to an id below nids; ids are not contiguous: a line
    whose intron is new gets the id of its position among all lines added. '''

    def __init__(self):
        self.introns = {}
        self.nids = 0
        self.files = {}
        self.groups = []

    def add(self, index, strands, starts, ends, reads):
        lines = self.files.get(index)
        if lines is None:
            lines = self.files[index] = FileLines()
        lines.ids.extend(map(self.introns.setdefault, zip(strands, starts, ends),
                             itertools.count(self.nids)))
        self.nids += len(reads)
        lines.reads.extend(reads)
        for s in dict.fromkeys(strands):
            if s not in lines.strands:
                lines.strands.append(s)

    def sums(self):
        ''' The strands in the order the block's files first count them, the first file
        that counts each, and per strand the introns' A, B and pooled reads. '''

        strands, firsts = [], []
        for index in sorted(self.files):
            for s in self.files[index].strands:
                if s not in strands:
                    strands.append(s)
                    firsts.append(index)
        tot = [0] * self.nids
        for lines in self.files.values():
            for i, r in zip(lines.ids, lines.reads):
                tot[i] += r
        groups = dict((s, ([], [], [])) for s in strands)
        for (s, A, B), i in self.introns.items():
            group = groups[s]
            group[0].append(A)
            group[1].append(B)
            group[2].append(i)
        self.groups = [ids for A, B, ids in groups.values()]
        return strands, firsts, [(int_array(A), int_array(B),
                                  array.array('q', map(tot.__getitem__, ids)))
                                 for A, B, ids in groups.values()]

    def slots(self, rows, spare):
        ''' The row of each id, given rows per strand in the order of sums(). '''

        slots = [spare] * self.nids
        for ids, group_rows in zip(self.groups, rows):
            for i, row in zip(ids, group_rows):
                slots[i] = row
        return slots


class PipeWorker(object):
    ''' A forked worker of pipeline_run, reading a block of junction files a chromosome
    at a time and answering the parent's requests between files. '''

    def __init__(self, outbox, in_conn, libs, options):
        self.out = outbox
        self.inbox = in_conn
        self.nochromcheck = options.nochromcheck
        self.maxlen = int(options.maxintronlen)
        self.readers = [JuncReader(lib, self) for lib in libs]
        self.order = FileOrder()
        self.chroms = {}
        self.released = set()
        self.finished = False

    def counted(self, chrom):
        return bool(self.nochromcheck) or chrom in CHROMOSOME_BYTES

    def run(self):
        for i, reader in enumerate(self.readers):
            reader.advance(None, None, i)
            self.visit(i, reader)
        self.out.put(('heads', [reader.head for reader in self.readers]))
        chrom = self.next_chrom()
        while chrom is not None:
            if chrom in self.released:
                raise PipelineFallback("chromosome orders disagree")
            counts = self.chroms.get(chrom)
            if counts is None:
                counts = self.chroms[chrom] = ChromCounts()
            moved = []
            for i, reader in enumerate(self.readers):
                if reader.head == chrom:
                    reader.passed.add(chrom)
                    reader.advance(chrom, counts, i)
                    self.visit(i, reader)
                    moved.append((i, reader.head))
                    self.serve(False)
            self.out.put(('step', chrom, moved))
            chrom = self.next_chrom()
        self.out.put(('parsed',))
        while not self.finished:
            self.serve(True)
        self.out.put(('done',))

    def visit(self, i, reader):
        head = reader.head
        if head is not None:
            if head in reader.passed or head in self.released:
                raise PipelineFallback("a file leaves a chromosome and returns to it")
            self.order.visit(i, head)

    def next_chrom(self):
        ''' The head to read next: one no other head is known to precede, at the files
        that have read the fewest chromosomes, then held by the most files. '''

        files = {}
        for i, reader in enumerate(self.readers):
            if reader.head is not None:
                files.setdefault(reader.head, []).append(i)
        if not files:
            return None
        heads = list(files)
        free = [h for h in heads if not any(self.order.before(o, h) for o in heads if o != h)]
        return min(free or heads, key=lambda h: (min(len(self.readers[i].passed) for i in files[h]),
                                                 -len(files[h]), files[h][0]))

    def serve(self, wait):
        while wait or self.inbox.poll():
            wait = False
            message = self.inbox.recv()
            if message[0] == 'sums':
                self.send_sums(message[1])
            elif message[0] == 'map':
                self.send_counts(*message[1:])
            else:
                self.finished = True
                return

    def send_sums(self, chrom):
        if any(reader.head == chrom for reader in self.readers):
            raise PipelineFallback("chromosome orders disagree")
        self.released.add(chrom)
        self.out.put(('sums', chrom) + self.chroms[chrom].sums())

    def send_counts(self, chrom, rows, sizes, ends):
        counts = self.chroms.pop(chrom)
        slots = counts.slots(rows, ends[-1] if ends else 0)
        files = []
        for i in range(len(self.readers)):
            lines = counts.files.get(i)
            if lines is None:
                files.append(None)
            else:
                files.append(file_counts(lines.ids, lines.reads, slots, sizes, ends))
        self.out.put(('counts', chrom, files))

    def text_record(self, text):
        ''' (chrom, strand, A, B, reads) of a counted line read as text, None for a line
        neither pool_junc_reads nor sort_junctions counts; raises for a line either one
        reports, prints or fails on. '''

        fields = text.split()
        if len(fields) != 12:
            raise PipelineFallback("a line without 12 fields")
        chrom, A, B, dot, counts, strand, rA, rb, rgb, blockCount, blockSize, blockStarts = fields
        if int(blockCount) > 2:
            raise PipelineFallback("a line with more than 2 blocks")
        if not self.nochromcheck and (chrom not in CHROMOSOMES):
            return None
        Aoff, Boff = blockSize.split(",")
        A, B = int(A)+int(Aoff), int(B)-int(Boff)+1
        reads = int(counts)
        if strand == "?" or B-A > self.maxlen:
            return None
        if not -2**31 <= reads < 2**31:
            raise PipelineFallback("a count out of range")
        return chrom.encode('ascii'), strand.encode('ascii'), A, B, reads


class JuncReader(object):
    ''' One junction file of a PipeWorker, read in binary chunks.

    A run of lines that start with one chromosome is checked and parsed in one batch
    (pipe_batch); any other line is decoded and parsed as text exactly as
    pool_junc_reads and sort_junctions parse it (pipe_text). '''

    def __init__(self, path, worker):
        self.worker = worker
        self.handle = gzip.open(path, 'rb') if path.endswith('.gz') else open(path, 'rb')
        self.lines = []
        self.pos = 0
        self.tail = b''
        self.texts = []
        self.head = None
        self.passed = set()

    def fill(self):
        data = self.handle.read(PIPE_CHUNK)
        if not data:
            if not self.tail:
                return False
            self.lines, self.tail = [self.tail], b''
        else:
            self.lines = (self.tail + data).split(b'\n')
            self.tail = self.lines.pop()
        self.pos = 0
        return True

    def advance(self, chrom, counts, index):
        ''' Counts lines of chrom until a counted line of another chromosome, which
        becomes the head and stays unread, or the end of the file (head None). '''

        worker = self.worker
        prefix = None if chrom is None else chrom + b'\t'
        while True:
            if self.texts:
                record = worker.text_record(self.texts[0])
                if record is not None and record[0] != chrom:
                    self.head = record[0]
                    return
                del self.texts[0]
                if record is not None:
                    counts.add(index, [record[1]], [record[2]], [record[3]], [record[4]])
                continue
            if self.pos == len(self.lines):
                if not self.fill():
                    self.handle.close()
                    self.head = None
                    return
                continue
            line = self.lines[self.pos]
            if prefix is not None and line.startswith(prefix):
                self.take_run(chrom, counts, index)
                continue
            tab = line.find(b'\t')
            if tab > 0 and not worker.counted(line[:tab]):
                self.take_run(line[:tab], None, index)
                continue
            self.texts = pipe_text(line)
            self.pos += 1

    def take_run(self, chrom, counts, index):
        lines, pos = self.lines, self.pos
        prefix = chrom + b'\t'
        end = run_end(lines, pos, prefix)
        batch = pipe_batch(lines[pos:end], chrom, counts is not None, self.worker.maxlen)
        if batch is False:
            end = pos + 1
            while end < len(lines) and lines[end].startswith(prefix):
                end += 1
            batch = pipe_batch(lines[pos:end], chrom, counts is not None, self.worker.maxlen)
        if batch is None or batch is False:
            self.texts = [text for line in lines[pos:end] for text in pipe_text(line)]
        elif counts is not None:
            counts.add(index, *batch)
        self.pos = end


def run_end(lines, pos, prefix):
    ''' The end of the run of lines starting with prefix at pos, found by bisection;
    pipe_batch catches a run that does not hold together. '''

    hi = len(lines)
    if lines[hi - 1].startswith(prefix):
        return hi
    lo = pos + 1
    while lo < hi:
        mid = (lo + hi) // 2
        if lines[mid].startswith(prefix):
            lo = mid + 1
        else:
            hi = mid
    return lo


def pipe_batch(lines, chrom, count, maxlen):
    ''' Checks and parses lines that start with chrom and a tab in one go.

    Returns None when a line must be read as text (str.split() could see other fields
    than a split on tabs, or it is not a 2-block junction), False when a line does not
    start with chrom, and otherwise, if count is set, the (strands, A, B, reads) of the
    lines pool_junc_reads counts. '''

    n = len(lines)
    blob = b'\t\n\t'.join(lines)
    if (not blob.isascii() or any(byte in blob for byte in SPLIT_BYTES) or b'\t\t' in blob
            or blob.endswith(b'\t')):
        return None
    # a b'\n' field between lines: every line has exactly 12 fields
    fields = blob.split(b'\t')
    if len(fields) != 13 * n - 1 or fields[12::13].count(b'\n') != n - 1:
        return None
    if fields[0::13].count(chrom) != n:
        return False
    if fields[9::13].count(b'2') != n:
        return None
    if not count:
        return ()

    sizes = fields[10::13]
    aoff, boff = {}, {}
    for size in set(sizes):
        parts = size.split(b',')
        if len(parts) != 2:
            return None
        aoff[size] = int(parts[0])
        boff[size] = int(parts[1]) - 1
    if len(aoff) == 1:
        aoffs, boffs = itertools.repeat(aoff[sizes[0]]), itertools.repeat(boff[sizes[0]])
    else:
        aoffs, boffs = map(aoff.__getitem__, sizes), map(boff.__getitem__, sizes)
    starts = list(map(operator.add, map(int, fields[1::13]), aoffs))
    ends = list(map(operator.sub, map(int, fields[2::13]), boffs))
    reads = list(map(int, fields[4::13]))
    strands = fields[5::13]
    if b'?' in strands or max(map(operator.sub, ends, starts)) > maxlen:
        keep = [s != b'?' and B - A <= maxlen for s, A, B in zip(strands, starts, ends)]
        strands = list(itertools.compress(strands, keep))
        starts = list(itertools.compress(starts, keep))
        ends = list(itertools.compress(ends, keep))
        reads = list(itertools.compress(reads, keep))
    return strands, starts, ends, array.array('i', reads)


def pipe_text(line):
    ''' The lines text mode reads from one line of bytes (its newline removed). '''

    if b'\r' in line:
        return io.TextIOWrapper(io.BytesIO(line + b'\n'), encoding='utf-8').readlines()
    return [line.decode('utf-8') + '\n']


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

    parser.add_option("-c", "--cluster", dest="cluster", default = None,
                  help="refined cluster file when clusters are already made")

    parser.add_option("-k", "--nochromcheck", dest="nochromcheck", default = False,
                  help="Don't check that the chromosomes are well formated e.g. chr1, chr2, ..., or 1, 2, ...")

    parser.add_option("-C", "--includeconst", dest="const", action="store_true", default = False,
                  help="also include constitutive introns")

    (options, args) = parser.parse_args()

    if options.juncfiles is None:
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
