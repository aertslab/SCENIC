
# SCENIC

SCENIC is an R package to infer Gene Regulatory Networks and cell types from single-cell RNA-seq data. 

For more details and installation instructions see the tutorials: [Introduction and setup](https://rawcdn.githack.com/aertslab/SCENIC/701cc7cc4ac762b91479b3bd2eaf5ad5661dd8c2/inst/doc/SCENIC_Setup.html), and [running SCENIC](https://rawcdn.githack.com/aertslab/SCENIC/0a4c96ed8d930edd8868f07428090f9dae264705/inst/doc/SCENIC_Running.html).

The output from the examples is available at: http://scenic.aertslab.org/examples/ . 

Frequently asked questions: [FAQ](https://github.com/aertslab/SCENIC/blob/master/vignettes/FAQ.md)

### News

24/01/2019
- [Tutorial](https://rawcdn.githack.com/aertslab/SCENIC/0a4c96ed8d930edd8868f07428090f9dae264705/inst/doc/importing_pySCENIC.html) for importing results from [loom](http://scope.aertslab.org/) files/[pySCENIC](http://pyscenic.readthedocs.io)

20/06/2018
- Added function `export2scope()` (see http://scope.aertslab.org/). Version bump to 1.0.

01/06/2018
- Updated SCENIC pipeline to support the new version of RcisTarget and AUCell.

01/05/2018
- [RcisTarget](https://bioconductor.org/packages/RcisTarget) is now available in Bioconductor. The new databases can be downloaded from [https://resources.aertslab.org/cistarget/]. 

30/03/2018 New releases:
- [pySCENIC](http://pyscenic.readthedocs.io): lightning-fast python implementation of the SCENIC pipeline.
- [Arboreto](https://arboreto.readthedocs.io/) package including *GRNBoost2* and scalable GENIE3: Easy to install Python library that supports distributed computing. It allows fast co-expression module inference (Step1) on large datasets, compatible with both, the R and python implementations of SCENIC.
- [Drosophila databases](https://resources.aertslab.org/cistarget/) for RcisTarget
