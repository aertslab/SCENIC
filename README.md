
# SCENIC

**SCENIC (Single-Cell rEgulatory Network Inference and Clustering)** is a computational method to infer Gene Regulatory Networks and cell types from single-cell RNA-seq data. 

The description of the method and some usage examples are available in [Nature Methods (2017)](https://www.nature.com/articles/nmeth.4463).

There are currently **implementations** of SCENIC in R (this repository), and in Python. If you don't have a strong prefference for using R, we would recommend to check out the [SCENIC protocol repository](https://github.com/aertslab/SCENICprotocol/), which contains the *Nextflow workflow*, and *Python/Jupyter notebooks* to easily run SCENIC (highly recommended for running it in batch or bigger datasets). The **output** from any of the implementations can then be explored either in R, Python or [SCope](http://scope.aerslab.org) (a web interface).

For more details and installation instructions on running SCENIC in `R` see the **tutorials**:
  - [Introduction and setup](http://htmlpreview.github.io/?https://github.com/aertslab/SCENIC/blob/master/inst/doc/SCENIC_Setup.html)
  - [Running SCENIC](http://htmlpreview.github.io/?https://github.com/aertslab/SCENIC/blob/master/inst/doc/SCENIC_Running.html)
  - The output from these examples is available at: [https://scenic.aertslab.org/scenic_paper/examples/](https://scenic.aertslab.org/scenic_paper/examples/)

Frequently asked questions: [FAQ](https://github.com/aertslab/SCENIC/blob/master/vignettes/FAQ.md)

### News

2020/06/26:
- The **SCENICprotocol** including the Nextflow workflow, and `pySCENIC` notebooks are now officially released. For details see the [Github repository](https://github.com/aertslab/SCENICprotocol/), and the associated publication in [Nature Protocols](https://doi.org/10.1038/s41596-020-0336-2).

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
