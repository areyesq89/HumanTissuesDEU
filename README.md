# HumanDEU

This repository contains the code and pre-computed data objects to reproduce each figure, statistic 
and tables from the paper:

* A Reyes and W Huber. Alternative start and termination sites of transcription drive most transcript isoform 
differences across human tissues. Nucleic Acids Research, 2017. 
doi: [10.1093/nar/gkx1165](https://www.doi.org/10.1093/nar/gkx1165)

The code is organized as an **R** package containing 6 different vignettes that reproduce 
different aspects of the manuscript:

* Vignette `01_precomputedObjects.Rmd` contains the code used to generate the pre-computed objects provided in this package. 
* Vignette `02_figure1Analysis.Rmd` contains analyses and code related to Figure 1. 
* Vignette `03_figure2Analysis.Rmd` contains analyses and code related to Figure 2.
* Vignette `04_figure3Analysis.Rmd` contains analyses and code related to Figure 3.
* Vignette `05_figure4Analysis.Rmd` contains analyses and code related to Figure 4.
* Vignette `06_figure5Analysis.Rmd` contains analyses and code related to Figure 5.
* Vignette `99_revisions.Rmd` contains plots that were generated to convince reviewers to accept our paper.

Since Github repositories have size limits, the **R** data objects of this package were deposited in
zenodo. In order to reproduce the results of the manuscript, you will have to build the **R** package by 
following the next steps in a terminal:

**1. Clone this repository:**

```
git clone https://github.com/areyesq89/HumanTissuesDEU.git
```

**2. Download **R** objects from zenodo placing them in their corresponding data directories inside the package:**

```
for i in `cat HumanTissuesDEU/inst/files/Robjects.txt`
do
curl -o HumanTissuesDEU/data/$i https://zenodo.org/record/2583270/files/$i
done

for i in `cat HumanTissuesDEU/inst/files/extdata.txt`
do
curl -o HumanTissuesDEU/inst/extdata/$i https://zenodo.org/record/2583270/files/$i
done

```

**3. Make sure that the files were downloaded correctly by comparing the *md5sum* results. You could do this in **R** by doing:**

```
uploaded <- read.table("HumanTissuesDEU/inst/files/md5sum.check")
colnames(uploaded) <- c("md5up", "file")
files <- list.files("HumanTissuesDEU", pattern="(RData|sqlite)$", full.names=TRUE, recursive=TRUE)
downloaded <- tools::md5sum(files)
names(downloaded) <- basename(names( downloaded ))
downloaded <- data.frame( md5down=downloaded, file=names(downloaded) )
isOK <- all( with( merge( downloaded, uploaded ), md5down == md5up ) )
if( !isOK ){ 
   stop("Files were not downloaded correctly, please try again!") 
} else { 
   message("Files downloaded succesfully")
}
```

**4. Build and install the package:**

```
R CMD build HumanTissuesDEU
R CMD INSTALL HumanDEU_0.0.99.tar.gz
```

If the build of the package was successful, the vignettes of the paper were compiled and the paper was reproduced! Now, you should also be able to go through each vignette and reproduce the paper. If you run into problems, please send me an e-mail and I will try to respond shortly.

**Note:** The code to plot the sashimi plots of this manuscript manuscript is available though this Github repository. These functions input bam files and, unfortunately, we are unable to provide these bam files since they contain potentially identifiable data. If you have access to the **GTEX** alignment data, the functions are expecting the bam files to have the format `/path/to/directory/alignments/<SRRID>/<SRRID>_Aligned.sortedByCoord.out.bam`. If you define an environent variable `export gtex=/path/to/directory`, the vignette will try to find the files and reproduce the sashimi plots. 

**Note 2:** You will need a relatively big machine (>20Gb of RAM) to reproduce the analysis.
