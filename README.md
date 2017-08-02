# leafcutter

<img src="./docs/logo.png" width="200"> **Annotation-free quantification of RNA splicing.**

*Yang I Li, David A Knowles, Jonathan K Pritchard.*

Leafcutter quantifies RNA splicing variation using short-read RNA-seq data. The core idea is to leverage spliced reads (reads that span an intron) to quantify (differential) intron usage across samples. The advantages of this approach include
* easy detection of novel introns
* modeling of more complex splicing events than exonic PSI
* avoiding the challenge of isoform abundance estimation
* simple, computationally efficient algorithms scaling to 100s or even 1000s of samples

For details please see our [bioRxiv preprint](http://biorxiv.org/content/early/2016/03/16/044107)

Full documentation is available at <http://davidaknowles.github.io/leafcutter/>

We've developed a leafcutter [shiny](https://shiny.rstudio.com/) app for visualizing leafcutter results: you can view an example [here](https://leafvis.shinyapps.io/leafvis/). This shows leafcutter differential splicing results for a comparison of 10 brain vs. 10 heart samples (5 male, 5 female in each group) from [GTEx](https://www.gtexportal.org/home/). 
