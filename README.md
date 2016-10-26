# leafcutter
LeafCutter: Annotation-free quantification of RNA splicing. 

Yang I Li, David A Knowles, Jonathan K Pritchard. 

For details please see
http://biorxiv.org/content/early/2016/03/16/044107

For a (hopefully) complete example of the complete pipeline, take a look at
```
example_data/worked_out_example.sh
```

To compile the R package to perform differential splicing analysis and make junction plots you can either...

## Option 1. Install from source
```
cd leafcutter
R CMD INSTALL --build .
```
Note you'll need the following R packages: `Rcpp, rstan, foreach, ggplot2, R.utils, gridExtra, reshape2, Hmisc, dplyr, doMC, optparse`. 

## Option 2. Install using devtools

This has the advantage of installing the required package dependencies for you. 
```
library(devtools)
install_github("davidaknowles/leafcutter/leafcutter")
```
