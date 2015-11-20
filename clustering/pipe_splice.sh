#juncFiles="../run/intron_files/*.junc"
#juncFiles="../run/gtex/processed/*.junc ../run/evol/*recomp"        
#juncFiles="../run/evol/*recomp" 
juncFiles="../run/evol/*brain*recomp ../run/gtex/processed/*Brain*.junc"
#juncFiles="../run/gtex/processed/male_Lung_SRR613306_2.junc ../run/gtex/processed/male_Lung_SRR613306_1.junc"
outPrefix="../run/test"

maxIntronLen=100000
minReadsClu=50

echo "Checking junction files..."
#python check_files.py $outPrefix\_chrom $juncFiles

echo "Sorting junction files..."
#python sort_junc.py $outPrefix $juncFiles

echo "Pooling individuals to identify intron clusters..."
#python get_intron_splicing_pooled.py $outPrefix $maxIntronLen $juncFiles

echo "Spliting intron clusters..."
#python refine_intron_splicing.py $outPrefix $minReadsClu

echo "Getting per individuals counts..."
python get_intron_splicing_perind.py $outPrefix

echo "Assigning genes to clusters..."
#python get_cluster_gene.py $outPrefix\_perind.counts.gz
