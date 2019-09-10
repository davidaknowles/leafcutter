#!/bin/bash

# wget bam files (4Gb!)
wget https://www.dropbox.com/s/pni1zq5y6cr4tx5/example_geuvadis.tar.gz?dl=0 -O example_geuvadis.tar.gz
tar -xvf example_geuvadis.tar.gz

# download gencode annotation
wget https://web.stanford.edu/~dak33/leafcutter/gencode.v19.annotation.gtf.gz

## Differential splicing analysis using wrapper script
python run_ds.py -A run/geuvadis/RNA.NA06984_CEU.chr1.bam,run/geuvadis/RNA.NA06985_CEU.chr1.bam,run/geuvadis/RNA.NA06986_CEU.chr1.bam,run/geuvadis/RNA.NA06989_CEU.chr1.bam -B run/geuvadis/RNA.NA18486_YRI.chr1.bam,run/geuvadis/RNA.NA18487_YRI.chr1.bam,run/geuvadis/RNA.NA18488_YRI.chr1.bam,run/geuvadis/RNA.NA18489_YRI.chr1.bam -d .. -a gencode.v19.annotation.gtf.gz -t tmpdir

## Splicing QTL analysis using wrapper script
ls run/geuvadis/*.bam > bamfiles.txt
python run_sQTL.py -b bamfiles.txt -d .. -a gencode.v19.annotation.gtf.gz -t sqtl_tmpdir
