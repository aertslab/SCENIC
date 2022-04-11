
*Frequently asked questions (FAQ):*

  *Interpretation of results:*
   - Missing relevant TFs / Known motifs are not annotated to the TF
   - Is it possible to use pseudotime/trajectory inference methods on the output of SCENIC? (e.g. to explore how TF activity changes over time)
   - [R] What is the meaning of the sufix '_extended' and the number in parenthesis in the regulon name?
   
  *Modifying the pipeline:*
   - Is it possible to apply SCENIC to other organisms?
   - Can I use my own list of TFs and targets (e.g. from an external analysis or resource) with AUCell and skip the network inference steps?
   - Is it possible to use another co-expression network tool for the initial TF-target links? 
   
   *Technical issues:*
   - Cannot download the databases (webserver is down)
   - R vs Python version (and VSN-pipelines)
   - Previous versions of SCENIC & SCENIC for previous versions of R/Bioconductor
   - Errors reading the .feather databases
   - Other issues or questions

## Interpretation of results

### Missing relevant TFs / Known motifs are not annotated to the TF

This thread might be related to your question: https://github.com/aertslab/SCENIC/issues/14#issuecomment-357009514

### Is it possible to use pseudotime/trajectory inference methods on the output of SCENIC? (e.g. to explore how TF activity changes over time)
Yes, the regulon activity matrix can be used as input for other methods, such as dimensionality reduction (e.g. t-SNE/UMAP, difussion maps) or pseudotime/trajectory analysis. 

However, you should check the assumptions/requirements of the input data of the specific tool (e.g. branches, expectations in regards to distribution and continuous/discrete values, etc...).

### [R] What is the meaning of the sufix '_extended' and the number in parenthesis in the regulon name?

The numbers in parenthesis in the regulon names indicate the number of genes in the regulon (i.e. "Dlx (35g)").

The "extended" regulons include motifs that have been linked to the TF by lower confidence annotations (e.g. inferred by motif similarity). 

The full explanation is available in the vignette *[detailedStep_2_createRegulons.Rmd](https://rawcdn.githack.com/aertslab/SCENIC/a0a00644b2f3589a3e2bc65486fc5f6cc00f48e1/inst/doc/detailedStep_2_createRegulons.html)* : 

> 1.2 Annotate motifs to TFs
...
The annotations provided by the cisTarget databases can be divided into high-confidence or low-confidence, depending on the annotation source (annotated in the original database, inferred by orthology, or inferred by motif similarity). The main regulons only use the “high confidence” annotations, which by default are “direct annotation” and “inferred by orthology”. **The sufix _extended** in the regulon name indicates **lower confidence annotations** (by default “inferred by motif similarity”) are also used.

## Modifying the pipeline

### Is it possible to apply SCENIC to other organisms?
At the moment we only provide the motif databases (required for the standard pipeline) for human, mouse and fruit fly.

The second step in SCENIC's workflow is to perform motif enrichment analysis on the co-expression modules to prune them into regulons. The standard SCENIC pipeline uses the cisTarget framework to perform this analysis, which requires species-specific databases.
Therefore, to apply SCENIC to a new organism (such as zebrafrish) would require some modifications on this step. Either [creating the cisTarget motif databases](https://github.com/aertslab/create_cisTarget_databases) for the new organism (which implies defining regulatory regions and scoring thousands of motifs PWM), or using an alternative tool (i.e. instead of RcisTarget) for the motif analysis.

### Can I use my own list of TFs and targets (e.g. from an external analysis or resource) with AUCell and skip the network inference steps?

Yes, AUCell accepts any type of gene set as input (e.g. TFs and potential target genes, pathway... etc). Examples are available in the [AUCell vignette (R)](https://bioconductor.org/packages/release/bioc/vignettes/AUCell/inst/doc/AUCell.html). 

Once you have the AUC matrix, you can decide whether to continue SCENIC's workflow (e.g. binarization, heatmap... etc) or do an independent analysis.
Note that you will just need to adjust the interpretation of the results according to the input gene-sets... 

### Is it possible to use another co-expression network tool for the initial TF-target links? 

We chose GENIE3/GRNBoost for building the initial co-expression networks because it was the best peformer in previous benchmarks. 
However, other tools can also be used instead. Of course, the final network results will be different... but if you prefer to use another tool for any reason, you can certainly give it a try. 

Once you have a co-expression network (sets of co-expressed genes, and/or links between TF-targets), in R you can start SCENIC on the "create regulons step" to perform the motif enrichment analysis and prune it into a GRN (regulons). Either with the function `runSCENIC_2_createRegulons()` or running it manually ([detailedStep_2_createRegulons.Rmd](https://github.com/aertslab/SCENIC/blob/master/vignettes/detailedStep_2_createRegulons.Rmd)).

An example of the input format is available at: 
http://scenic.aertslab.org/examples/SCENIC_MouseBrain/int/2.1_tfModules_forMotifEnrichmet.Rds 
(e.g. co-expression gene sets as a list, including the candidate TF/regulator in the name).


## Technical issues

### Cannot download the databases (webserver is down)

An alternative mirror for downloading the RcisTarget databases and annotation files is available at: https://resources-mirror.aertslab.org/cistarget/

### R vs Python version

There are **implementations** of SCENIC in R (this repository), in Python ([pySCENIC](https://github.com/aertslab/pySCENIC)), as well as wrappers to automate analyses with Nextflow ([VSN-pipelines](https://vsn-pipelines.readthedocs.io/en/latest/)).

The results are equivalent across versions and provide output `.loom` files that can be explored in [SCope](http://scope.aertslab.org) or used as interface between R and Python.

If you have access to Nextflow and a container system (e.g. Docker or Singularity), we **recommend** to run SCENIC through the [VSN-pipeline](https://vsn-pipelines.readthedocs.io/en/latest/). This option is specially useful for running SCENIC on large datasets, in batch on multiple samples, and for users with limited experience in Python or R.

> This ([tutorial](http://htmlpreview.github.io/?https://github.com/aertslab/SCENIC/blob/master/Tutorials_JupyterNotebooks/SCENIC_tutorial_1-RunningVSN.html)) explains how to run VSN, which uses containers (in [Docker](https://cloud.docker.com/u/aertslab/repository/docker/aertslab/pyscenic) or [Singularity](https://www.singularity-hub.org/collections/2033)) including all the required dependencies, in a Nextflow pipeline ([VSN](https://github.com/aertslab/scenic-nf).

### [R] Previous versions of SCENIC & SCENIC for previous versions of R/Bioconductor

You can find previous versions of SCENIC (and which R/Bioconductor version they correspond to) in "releases": https://github.com/aertslab/SCENIC/releases

To install them, just add the tag: `devtools::install_github("aertslab/SCENIC@v1.1.1")`

### [R] Errors reading the .feather databases

Errors reading the databases usually happens when the files are incomplete/corrupt (e.g. by a failed download).

We recommended to download the databases using zsync_curl (https://resources.aertslab.org/cistarget/help.html). Once you have the files, make sure the sha256sum match the reported ones (https://resources.aertslab.org/cistarget/databases/sha256sum.txt). 

You can also check whether you are using the latest R feather package (versions older than 0.3.1 are more likely to crash due to corrupt downloads).

### Other issues or questions

You can check whether someone has already had a similar question in the [Discussions](https://github.com/aertslab/SCENIC/discussions/).

If you find a **bug** or you have a feature request, you can post it in the appropriate Github repo: [R](https://github.com/aertslab/SCENIC/issues?q=+is%3Aclosed) or [Python]([R](https://github.com/aertslab/pySCENIC/issues?q=+is%3Aclosed)) (please, before posting, check that it has not been already posted/solved).

