# HumanTissuesDEU

This repository contains the code and pre-computed data objects to reproduce each figure, statistic 
and tables from the paper:

* A Reyes and W Huber. Alternative start and termination sites of transcription drive most transcript isoform 
differences across human tissues. Nucleic Acids Research, 2017. 
doi: [10.1093/nar/gkx1165](https://www.doi.org/10.1093/nar/gkx1165)

The code is organized as an **R** package containing 6 different vignettes that reproduce 
different aspects of the manuscript:

* Vignette `01_precomputedObjects.Rmd` contains the code to generated the pre-computed objects from this package. 
* Vignette `02_figure1Analysis.Rmd` contains analyses and code related to figure 1. 
* Vignette `03_figure2Analysis.Rmd` contains analyses and code related to figure 2.
* Vignette`04_figure3Analysis.Rmd` contains analyses and code related to figure 3.
* Vignette `05_figure4Analysis.Rmd` contains analyses and code related to figure 4.
* Vignette `06_figure5Analysis.Rmd` contains analyses and code related to figure 5.
* Vignette `99_revisions.Rmd` contains plots that were generated to convince reviewers to accept our paper.

Since Github repositories have size limits, the **R** data objects of this package were deposited to
zenodo. In order to reproduce the results of the manuscript, you will have to build this package by 
following the next steps:

1. Clone this repository:

2. Download **R** objects from zenodo:

3. Copy the downloaded objects to the corresponding data directories of the package

4. Build and install the package

Now you should be able to open each vignette and reproduce each figure. If you run into any problem, please send me an e-mail and I will try to respond shortly.  
