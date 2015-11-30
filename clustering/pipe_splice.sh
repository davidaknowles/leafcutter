#juncFiles="../run/intron_files/*.junc"
#juncFiles="../run/gtex/processed/*.junc ../run/evol/*recomp"        
#juncFiles="../run/evol/*recomp" 
#juncFiles="../run/evol/*brain*recomp ../run/gtex/processed/*Brain*.junc"
juncFiles="~/scailscratch/splicing/processed/*.junc"
#juncFiles="../run/gtex/processed/male_Lung_SRR613306_2.junc ../run/gtex/processed/male_Lung_SRR613306_1.junc"
outPrefix="~/scailscratch/splicing/gtex_final"

mkdir outPrefix

maxIntronLen=100000
minReadsClu=50

echo "Checking junction files..."
#eval "python check_files.py $outPrefix\_chrom $juncFiles"

echo "Sorting junction files..."
#eval "python sort_junc.py $outPrefix $juncFiles"

echo "Pooling individuals to identify intron clusters..."
eval "python get_intron_splicing_pooled.py $outPrefix $maxIntronLen $juncFiles"

echo "Spliting intron clusters..."
#python refine_intron_splicing.py $outPrefix $minReadsClu

echo "Getting per individuals counts..."
#python get_intron_splicing_perind.py $outPrefix

echo "Assigning genes to clusters..."
#python get_cluster_gene.py $outPrefix\_perind.counts.gz
