
import subprocess

def run(cmd, max_minutes = 6000):
    import time
    sys.stderr.write("Running cmd: %s\n"%cmd)
    p = subprocess.Popen(cmd ,
                         shell=True,
                         stdin=subprocess.PIPE,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE,
                         close_fds=True)


    (file_stdin, file_stdout, file_stderr) = (p.stdin, p.stdout, p.stderr)

    t = 0
    r = ''
    e = ''
    while t < 60*max_minutes and p.poll() is None:
        time.sleep(1)  # (comment 1)
        t += 1
        r += file_stdout.read()
        e += file_stderr.read()
    r += file_stdout.read()
    e += file_stderr.read()

    
    file_stdin.close()
    #lines = file_stdout.read()
    lines_stderr = file_stderr.read()
    exit_code = file_stdout.close()
    
    file_stdout.close()
    file_stderr.close()
    return (r, e, exit_code)













if __name__ == "__main__":
    import sys, os
    from optparse import OptionParser

    parser = OptionParser()

    parser.add_option("-o", "--outprefix", dest="outprefix", default = 'leafcutter',
                      help="output prefix (default leafcutter)")

    parser.add_option("-t", "--tempdir", dest="tmpdir", default='./tmp/',
                      help="write to directory (default ./)")

    parser.add_option("-d", "--leafdir", dest="leafd", default='./',
                      help="LeafCutter directory")
    
    parser.add_option("-l", "--maxintronlen", dest="maxintronlen", default = 100000,
                  help="maximum intron length in bp (default 100,000bp)")

    parser.add_option("-m", "--minclureads", dest="minclureads", default = 30,
                  help="minimum reads in a cluster (default 30 reads)")

    parser.add_option("-p", "--mincluratio", dest="mincluratio", default = 0.001,
                  help="minimum fraction of reads in a cluster that support a junction (default 0.001)")

    parser.add_option("-a", "--annot", dest="annotation", default = "gencode.v19.annotation.gtf.gz",
                  help="path of annotation file e.g. ~/tools/leafcutter/clustering/gencode.v19.annotation.gtf.gz")

    parser.add_option("-A", "--A", dest="bamA",
                  help="bam files from condition A.")

    parser.add_option("-B", "--B", dest="bamB",
                  help="bam files from condition B.")

    

    (options, args) = parser.parse_args()


    try: open(options.leafd+"/leafcutter/R/differential_splicing.R")
    except:
        sys.stderr.write("Please specify correct LeafCutter directory e.g. -d tools/leafcutter/.\n")
        exit(0)

    if options.bamA == None or options.bamB == None:
        sys.stderr.write("Error: no bam file provided...\n")
        exit(0)

    bamA, bamB = options.bamA.split(','), options.bamB.split(',')

    # create tmp file
    try: os.mkdir(options.tmpdir)
    except: pass

    # check samtools

    sys.stderr.write("processing bam files...\n")
    fout = file("%s/junction_files.txt"%options.tmpdir,'w')
    for bam in bamA+bamB:
        bedfile = "%s/%s.bed"%(options.tmpdir,bam.split('/')[-1])
        juncfile = "%s/%s.junc"%(options.tmpdir,bam.split('/')[-1])
        fout.write(juncfile+'\n')
        try: open(juncfile)
        except: pass
        else:
            sys.stderr.write("%s exists..skipping\n"%juncfile)
            continue
        print run("samtools view %s | python %s/scripts/filter_cs.py | %s/scripts/sam2bed.pl --use-RNA-strand - %s"%(bam, options.leafd, options.leafd,bedfile))[1]
        print run("%s/scripts/bed2junc.pl %s %s; rm %s"%(options.leafd,bedfile,juncfile, bedfile))[1]
        
    fout.close()
    
    print run("python %s/clustering/leafcutter_cluster.py -j %s/junction_files.txt -m %s -o %s -l %s -r %s -p %s"%(options.leafd,options.tmpdir,options.minclureads, options.outprefix,str(options.maxintronlen), options.tmpdir,str(options.mincluratio)))[1]

    if options.annotation != None:
        print run("python %s/clustering/get_cluster_gene.py %s %s/%s_perind.counts.gz"%(options.leafd,options.annotation, options.tmpdir,options.outprefix))[1]
        
    
    
    
    fout = file("%s/ds_test"%(options.tmpdir),'w')
    for bam in bamA:
        fout.write("%s group1\n"%bam.split("/")[-1])
    for bam in bamB:
        fout.write("%s group2\n"%bam.split("/")[-1])
    fout.close()

    
    print run("Rscript %s/scripts/leafcutter_ds.R --num_threads 1 -i 3 %s/%s_perind_numers.counts.gz %s/ds_test"%(options.leafd,options.tmpdir,options.outprefix,options.tmpdir))[1]

