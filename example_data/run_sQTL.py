
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
                      help="Output prefix (default leafcutter)")

    parser.add_option("-t", "--tempdir", dest="tmpdir", default='./tmp/',
                      help="Where to output files (default ./)")

    parser.add_option("-d", "--leafdir", dest="leafd", default='./',
                      help="Top-level LeafCutter directory")
    
    parser.add_option("-l", "--maxintronlen", dest="maxintronlen", default = 100000,
                  help="Maximum intron length in bp (default 100,000bp)")

    parser.add_option("-m", "--minclureads", dest="minclureads", default = 30,
                  help="Minimum reads in a cluster (default 30 reads)")

    parser.add_option("-p", "--mincluratio", dest="mincluratio", default = 0.001,
                  help="Minimum fraction of reads in a cluster that support a junction (default 0.001)")

    parser.add_option("-a", "--annot", dest="annotation", default = None,
                  help="[optional] Path of annotation GTF file e.g. ~/tools/leafcutter/clustering/gencode.v19.annotation.gtf.gz")

    parser.add_option("-b", "--bams", dest="bam",
                  help="Text file listing bam files to quantify")
    
    (options, args) = parser.parse_args()

    try: open(options.leafd+"/clustering/leafcutter_cluster.py")
    except:
        sys.stderr.write("Please specify correct LeafCutter directory e.g. -d tools/leafcutter/.\n")
        exit(0)

    if options.bam == None:
        sys.stderr.write("Error: no bam file provided...\n")
        exit(0)

    bams = open(options.bam).readlines()

    # create tmp file directory
    try: os.mkdir(options.tmpdir)
    except: pass

    # (should check if samtools are installed)

    sys.stderr.write("processing bam files...\n")
    fout = file("%s/junction_files.txt"%options.tmpdir,'w')
    for bam in bams:
        bam = bam.strip()
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
        pass
    

    print run("python %s/scripts/prepare_phenotype_table.py %s/%s_perind.counts.gz"%(options.leafd,options.tmpdir,options.outprefix))

    sys.stdout.write("\n*******fastQTL instructions (also see http://fastqtl.sourceforge.net/) *******\n")
    sys.stdout.write("\n(1) Prepare phenotypes: Use `sh %s/%s_perind.counts.gz_prepare.sh' to create index for fastQTL (requires tabix and bgzip).\n"%(options.tmpdir,options.outprefix))
    sys.stdout.write("(2) Prepare covariates: To take the top 5 PCs, use `head -6 %s/%s_perind.counts.gz.PCs > %s/%s_perind.counts.gz.PCs.PC5\n"%(options.tmpdir,options.outprefix,options.tmpdir,options.outprefix))
    sys.stdout.write("(3) Prepare genotypes: bgzip + tabix your genotype (VCF) file > SNPs.MAF05.txt.gen.gz (not included)\n")
    sys.stdout.write("(4) Run fastQTL: Use e.g. `FastQTL/bin/fastQTL.static --vcf SNPs.MAF05.txt.gen.gz --bed %s/%s_perind.counts.gz.qqnorm_chr1.gz -out %s/%s_output_chr1 --chunk 1 1 --window 1e5 --cov %s/%s_perind.counts.gz.PCs.PC5\n\n\n"%(options.tmpdir,options.outprefix,options.tmpdir,options.outprefix,options.tmpdir,options.outprefix))
