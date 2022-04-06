
# SCENIC

**SCENIC (Single-Cell rEgulatory Network Inference and Clustering)** is a computational method to infer Gene Regulatory Networks and cell types from single-cell RNA-seq data. 

The description of the method and some usage examples are available in [Nature Methods (2017)](https://www.nature.com/articles/nmeth.4463).

There are currently **implementations** of SCENIC in R (this repository), in Python ([pySCENIC](https://github.com/aertslab/pySCENIC)), as well as wrappers to automate analyses with Nextflow ([VSN-pipelines](https://vsn-pipelines.readthedocs.io/en/latest/)).

The **output** from any of the implementations can be explored either in R, Python or [SCope](https://scope.aertslab.org) (a web interface).

### Tutorials

If you have access to Nextflow and a container system (e.g. Docker or Singularity), we **recommend** to run SCENIC through the VSN-pipeline. 
> This option is specially useful for running SCENIC on large datasets, or in batch on multiple samples. 

- [1. Run SCENIC from VSN](http://htmlpreview.github.io/?https://github.com/aertslab/SCENIC/blob/master/Tutorials_JupyterNotebooks/SCENIC_tutorial_1-RunningVSN.html) 
- [2. Explore SCENIC output (with SCope and R)](http://htmlpreview.github.io/?https://github.com/aertslab/SCENIC/blob/master/Tutorials_JupyterNotebooks/SCENIC_tutorial_2-ExploringOutput.html)

If you prefer to use **R** for the whole analysis, these are the main tutorials:
> The tutorials in R include a more detailed explanation of the workflow and source code. 
  - [Introduction and setup](http://htmlpreview.github.io/?https://github.com/aertslab/SCENIC/blob/master/inst/doc/SCENIC_Setup.html)
  - [Running SCENIC](http://htmlpreview.github.io/?https://github.com/aertslab/SCENIC/blob/master/inst/doc/SCENIC_Running.html)
  - The output from these examples is available at: [https://scenic.aertslab.org/scenic_paper/examples/](https://scenic.aertslab.org/scenic_paper/examples/)

**Python/Jupyter notebooks** with examples running SCENIC in different settings are available in the [SCENIC protocol repository](https://github.com/aertslab/SCENICprotocol/).

Frequently asked questions: [FAQ](https://github.com/aertslab/SCENIC/blob/master/vignettes/FAQ.md)

---

### News

2021/03/26:

- New tutorials to [run SCENIC from VSN](http://htmlpreview.github.io/?https://github.com/aertslab/SCENIC/blob/master/Tutorials_JupyterNotebooks/SCENIC_tutorial_1-RunningVSN.html) 
and [explore its output (with SCope and R)](http://htmlpreview.github.io/?https://github.com/aertslab/SCENIC/blob/master/Tutorials_JupyterNotebooks/SCENIC_tutorial_2-ExploringOutput.html)

- Tutorial to [create new databases](https://github.com/aertslab/create_cisTarget_databases)


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
