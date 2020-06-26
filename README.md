
# SCENIC

**SCENIC (Single-Cell rEgulatory Network Inference and Clustering)** is an R package to infer Gene Regulatory Networks and cell types from single-cell RNA-seq data. 

For more details and installation instructions see the tutorials:
  - [Introduction and setup](https://rawcdn.githack.com/aertslab/SCENIC/701cc7cc4ac762b91479b3bd2eaf5ad5661dd8c2/inst/doc/SCENIC_Setup.html)
  - [running SCENIC](https://rawcdn.githack.com/aertslab/SCENIC/0a4c96ed8d930edd8868f07428090f9dae264705/inst/doc/SCENIC_Running.html)

The output from the examples is available at: [http://scenic.aertslab.org/examples/](http://scenic.aertslab.org/examples/)
 

Frequently asked questions: [FAQ](https://github.com/aertslab/SCENIC/blob/master/vignettes/FAQ.md)

Use [pySCENIC](http://pyscenic.readthedocs.io) for a lightning-fast python implementation of SCENIC.


### News

2019/01/24:
  - [Tutorial](https://rawcdn.githack.com/aertslab/SCENIC/0a4c96ed8d930edd8868f07428090f9dae264705/inst/doc/importing_pySCENIC.html)
    for importing [pySCENIC](http://pyscenic.readthedocs.io) results in SCENIC by using [loom](http://scope.aertslab.org/) files.

2018/06/20:
  - Added function `export2scope()` (see http://scope.aertslab.org/).
  - Version bump to 1.0.

2018/06/01:
  - Updated SCENIC pipeline to support the new version of RcisTarget and AUCell.

2018/05/01:
  - [RcisTarget](https://bioconductor.org/packages/RcisTarget) is now available in Bioconductor.
  - The new databases can be downloaded from [https://resources.aertslab.org/cistarget/](https://resources.aertslab.org/cistarget/). 

2018/03/30: New releases
  - [pySCENIC](https://pyscenic.readthedocs.io): lightning-fast python implementation of the SCENIC pipeline.
  - [Arboreto](https://arboreto.readthedocs.io) package including **GRNBoost2** and scalable **GENIE3**:
      - Easy to install Python library that supports distributed computing.
      - It allows fast co-expression module inference (Step1) on large datasets, compatible with both, the R and python implementations of SCENIC.
  - [Drosophila databases](https://resources.aertslab.org/cistarget/) for RcisTarget.
