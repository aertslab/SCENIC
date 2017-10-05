# SCENIC

SCENIC is an R package to infer Gene Regulatory Networks and cell types from single-cell RNA-seq data. 

## Installation
To use SCENIC you will also need to install GENIE3, RcisTarget and AUCell. The current version of SCENIC (v0.1.6) requires AUCell 0.99.5 and RcisTarget 0.99.0. Newer versions are available in the corresponding Github repositories and will soon be available in Bioconductor. We will update SCENIC to work with these updated versions once they are available in Bioconductor (begining of November 2017).

To install the versions compatible with this version of SCENIC, you can use:
```
# You might need to install some of these dependencies first:
c("R.utils", "utils", "graphics", "stats", "data.table", "mixtools", "GSEABase", 
"SummarizedExperiment", "methods", "Biobase", "zoo", "DT", "NMF", "plotly", 
"BiocStyle", "rmarkdown", "doMC", " doRNG", " zoo", "doParallel", "foreach")

# GENIE3:
install.packages("https://bioconductor.org/packages/devel/bioc/src/contrib/GENIE3_0.99.7.tar.gz", repos=NULL)
# AUCell:
install.packages("http://scenic.aertslab.org/downloads/Rpackages/AUCell_0.99.5.tar.gz", repos=NULL)
# RcisTarget:
install.packages("http://scenic.aertslab.org/downloads/Rpackages/RcisTarget_0.99.0.tar.gz", repos=NULL)
```

RcisTarget databases (choose the appropiate organism):
```
# Human: install.packages("http://scenic.aertslab.org/downloads/databases/RcisTarget.hg19.motifDatabases.20k_0.1.1.tar.gz", repos=NULL)
# Mouse: install.packages("http://scenic.aertslab.org/downloads/databases/RcisTarget.mm9.motifDatabases.20k_0.1.1.tar.gz", repos=NULL)
```

Finally, to install SCENIC:
```
# SCENIC:
devtools::install_github("aertslab/SCENIC")
```

## Tutorials

The package includes tutorials describing the functionality and usage of SCENIC:
- [Introduction to SCENIC](https://htmlpreview.github.io/?https://github.com/aertslab/SCENIC/blob/master/inst/doc/SCENIC_Intro.html)
- [Running SCENIC](https://htmlpreview.github.io/?https://github.com/aertslab/SCENIC/blob/master/inst/doc/SCENIC_runningStep1nWrapper.html)

Advanced tutorials:
- [Step 1 - Part 2: Co-expression modules](https://htmlpreview.github.io/?https://github.com/aertslab/SCENIC/blob/master/inst/doc/Step1.2_CoexpressionModules.html)
- [Step 2: Regulons (DNA motif analysis with RcisTarget)](https://htmlpreview.github.io/?https://github.com/aertslab/SCENIC/blob/master/inst/doc/Step2_Regulons.html)
- [Step 3 - Part 1: Network activity in each cell (AUCell)](https://htmlpreview.github.io/?https://github.com/aertslab/SCENIC/blob/master/inst/doc/Step3.1_NwActivity.html)
- [Step 3 - Part 2: Binary network activity](https://htmlpreview.github.io/?https://github.com/aertslab/SCENIC/blob/master/inst/doc/Step3.2_BinaryNwActivity.html)
- [Step 4: GRN-based cell states (Clustering)](https://htmlpreview.github.io/?https://github.com/aertslab/SCENIC/blob/master/inst/doc/Step4_Clustering.html)

(examples for follow-up coming soon)


## Website

For more information, please visit http://scenic.aertslab.org
