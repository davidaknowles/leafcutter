# leafcutter
LeafCutter: Annotation-free quantification of RNA splicing. 

Yang I Li, David A Knowles, Jonathan K Pritchard. 

For details please see
http://biorxiv.org/content/early/2016/03/16/044107

## Installation

To compile the R package to perform differential splicing analysis and make junction plots you can either...

### Option 1. Install from source
```
cd leafcutter
R CMD INSTALL --build .
```
You'll need the following R packages: `Rcpp, rstan, foreach, ggplot2, R.utils, gridExtra, reshape2, Hmisc, dplyr, doMC, optparse`. 

### Option 2. Install using devtools

This has the advantage of installing the required package dependencies for you. 
```
library(devtools)
install_github("davidaknowles/leafcutter/leafcutter")
```

## Usage

For a (hopefully) complete example of the complete pipeline, take a look at
```
example_data/worked_out_example.sh
```

LeafCutter has two main components: 

1. Python code to 
   - generate intron excision counts from `junc` files (which can be obtained easily from `.bam` files)
   - group introns into clusters
2. `R` code to 
   * perform differential splicing (here, differential intron excision) analysis
   * plot (differentially spliced) clusters
   
### Step 1. Converting `bam`s to `junc`s

I'm skipping Step 0 which would be mapping `fastq` files, e.g. using STAR, to obtain `.bam` files. 

We provide a helper script `scripts/bam2junc.sh` to (you guessed it) convert `bam` files to `junc` files. This step uses the CIGAR strings in the `bam` to quantify the usage of each intron. 

`example_data/worked_out_example.sh` gives you an example of how to do this in batch:
```
for bamfile in `ls my_awesome_data/*.bam`
do
    echo Converting $bamfile to $bamfile.junc
    sh ../scripts/bam2junc.sh $bamfile $bamfile.junc
    echo $bamfile.junc >> awesome_juncfiles.txt
done
```

This step is pretty fast but if you have samples numbering in the 100s you might want to do this on a cluster. Note that we also make a list of the generated `junc` files in `awesome_juncfiles.txt`. 

### Step 2. 
