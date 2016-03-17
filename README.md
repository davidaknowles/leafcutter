# leafcutter
LeafCutter: Annotation-free quantification of RNA splicing. 

Yang I Li, David A Knowles, Jonathan K Pritchard. 

For details please see
http://biorxiv.org/content/early/2016/03/16/044107

To compile the R package to perform differential splicing analysis and make junction plots:
```
cd leafcutter
R CMD INSTALL --build .
```

Note you'll need the following R packages: `rstan, foreach, ggplot2, R.utils, gridExtra, reshape2, Hmisc`. 

For a (hopefully) complete example of the complete pipeline, take a look at
```
example_data/worked_out_example.sh
```
