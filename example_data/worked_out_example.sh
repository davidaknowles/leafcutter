# Convert bam to junction files
if [ -e test_juncfiles.txt ]
then
    rm test_juncfiles.txt
fi

# wget bam files

for bamfile in `ls ../../run/geuvadis/*chr1.bam`
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


for lib in `cut -f 5 -d'/' ../example_data/test_juncfiles.txt | cut -f 1-4 -d'.' | grep YRI`
do
    echo $lib YRI >> test_diff_intron.txt
done

for lib in `cut -f 5 -d'/' ../example_data/test_juncfiles.txt | cut -f 1-4 -d'.' | grep CEU`
do
    echo $lib CEU >> test_diff_intron.txt
done

../differential_splicing/leafcutter_ds.R ../example_data/testYRIvsEU_perind_numers.counts.gz ../example_data/test_diff_intron.txt
