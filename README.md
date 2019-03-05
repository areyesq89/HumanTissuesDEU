# HumanTissuesDEU

This repository contains the code and pre-computed data objects to reproduce each figure, statistic 
and tables from the paper:

* A Reyes and W Huber. Alternative start and termination sites of transcription drive most transcript isoform 
differences across human tissues. Nucleic Acids Research, 2017. 
doi: [10.1093/nar/gkx1165](https://www.doi.org/10.1093/nar/gkx1165)

The code is organized as an **R** package containing 6 different vignettes that reproduce 
different aspects of the manuscript:

* Vignette `01_precomputedObjects.Rmd` contains the code to generated the pre-computed objects from this package. 
* Vignette `02_figure1Analysis.Rmd` contains analyses and code related to Figure 1. 
* Vignette `03_figure2Analysis.Rmd` contains analyses and code related to Figure 2.
* Vignette `04_figure3Analysis.Rmd` contains analyses and code related to Figure 3.
* Vignette `05_figure4Analysis.Rmd` contains analyses and code related to Figure 4.
* Vignette `06_figure5Analysis.Rmd` contains analyses and code related to Figure 5.
* Vignette `99_revisions.Rmd` contains plots that were generated to convince reviewers to accept our paper.

Since Github repositories have size limits, the **R** data objects of this package were deposited to
zenodo. In order to reproduce the results of the manuscript, you will have to build this package by 
following the next steps in a terminal:

1. Clone this repository:

```

```

2. Download **R** objects from zenodo and copy them to the corresponding data directories of the package:

```

```

3. Build and install the package:

```

```

Now you should be able to go through each vignette and reproduce the paper. If you run into any problem, please send me an e-mail and I will try to respond shortly.

**Note:** The code to plot the sashimi plots presented in the manuscript is available though this Github repository. These functions input bam files and, unfortunately, we are unable to provide these bam files since these contain potentially identifiable data. 
**Note 2:** You will need a big computer (>20Gb of RAM) to reproduce the analysis.
