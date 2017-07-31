#!/bin/bash 

# wget bam files (4Gb!)
wget https://www.dropbox.com/s/pni1zq5y6cr4tx5/example_geuvadis.tar.gz?dl=0 -O example_geuvadis.tar.gz
tar -xvf example_geuvadis.tar.gz

## Differential splicing analysis using wrapper script
# load samtools/python

python run_ds.py -A run/geuvadis/RNA.NA06984_CEU.chr1.bam,run/geuvadis/RNA.NA06985_CEU.chr1.bam,run/geuvadis/RNA.NA06986_CEU.chr1.bam,run/geuvadis/RNA.NA06989_CEU.chr1.bam -B run/geuvadis/RNA.NA18486_YRI.chr1.bam,run/geuvadis/RNA.NA18487_YRI.chr1.bam,run/geuvadis/RNA.NA18488_YRI.chr1.bam,run/geuvadis/RNA.NA18489_YRI.chr1.bam -d ~/tools/leafcutter/ -a ../clustering/gencode.v19.annotation.gtf.gz -t tmpdir


## Splicing QTL analysis using wraper script
#load samtools/python

ls run/geuvadis/*.bam > bamfiles.txt
python run_sQTL.py -b bamfiles.txt -d ~/tools/leafcutter/ -a ~/tools/leafcutter/clustering/gencode.v19.annotation.gtf.gz -t sqtl_tmpdir



## Manually working through the LeafCutter example (helps you understand the steps and modify code/usage)...

# Convert bam to junction files                                                                                                                                                                                                                                                                                                                                           
if [ -e test_juncfiles.txt ]
then
    rm test_juncfiles.txt
fi
for bamfile in `ls run/geuvadis/*chr1.bam`
do
    echo Converting $bamfile to $bamfile.junc
    sh ../scripts/bam2junc.sh $bamfile $bamfile.junc
    echo $bamfile.junc >> test_juncfiles.txt
done

python ../clustering/leafcutter_cluster.py -j test_juncfiles.txt -m 50 -o testYRIvsEU -l 500000

if [ -e test_diff_intron.txt ]
then 
    rm test_diff_intron.txt 
fi


for lib in `cut -f 3 -d'/' ../example_data/test_juncfiles.txt | cut -f 1-4 -d'.' | grep YRI`
do
    echo $lib YRI >> test_diff_intron.txt
done

for lib in `cut -f 3 -d'/' ../example_data/test_juncfiles.txt | cut -f 1-4 -d'.' | grep CEU`
do
    echo $lib CEU >> test_diff_intron.txt
done

../scripts/leafcutter_ds.R --num_threads 4 ../example_data/testYRIvsEU_perind_numers.counts.gz ../example_data/test_diff_intron.txt

#../scripts/ds_plots.R -e ../leafcutter/data/gencode19_exons.txt.gz ../example_data/testYRIvsEU_perind_numers.counts.gz ../example_data/test_diff_intron.txt leafcutter_ds_cluster_significance.txt -f 0.05
