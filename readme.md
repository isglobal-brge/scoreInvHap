# scoreInvHap

## Introduction 

`r Rpackage("scoreInvHap")` genotypes inversions using SNP data. This method computes a similarity score between the available SNPs of an individual and experimental references for the three inversion-genotypes; NN: non-inverted homozygous, NI: inverted heterozygous and II: inverted homozygous. Individuals are classified into the inversion-genotypes with the highest score. There are other approaches to genotype inversions from SNP data: [inveRsion](https://bioconductor.org/packages/release/bioc/html/inveRsion.html), [invClust](http://nar.oxfordjournals.org/content/early/2015/02/05/nar.gkv073.full) or [PFIDO](https://github.com/MaxSalm/pfido). However, these approaches have limitations including:

* __they are based in population sample inferences__. When a new individual is included in the study the whole population sample needs to be reanalyzed. As they depend on dimensionality reduction methods, they are slow in the analysis of several individuals.  

* __they can only handle limited number of samples__. They need large samples sizes for accurate inferences. 
* __They are highly sensitive to SNP array densities__.  Their accuracy depends on how clear the population sample can be partitioned, which depends on the amount and quality of informative SNPs. 
 
* __They need external information to label the inversion-genotype clusters__. External validation of the clustering is performed at the end of the inference. If there are more than three clusters, it is not clear how which inversion genotype should be assigned to the clusters.  

`r Rpackage("scoreInvHap")` overcomes these difficulties by using a set of reference genotypes from the inversion of interest. The methods is a supervised classifier  that genotypes each individual according to their SNP similarities to the reference genotypes across the inverted segment. The classifier in particular uses the linkage desequillibrium (R^2^) between the SNPs and the inversion genoptypes, and the SNPs of each reference inversion-genotypes. 

The methods described in the paper **"scoreInvHap: Inversion genotyping for genome-wide association studies" by Carlos Ruiz-Arenas et al. (2019)**[[1]](#1)  publicly available at  **[PLOS Genetics](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1008203)** can be found in release [v 1.0.0](https://github.com/isglobal-brge/scoreInvHap/releases/tag/1.0.0).


## Data and Vignette

The package uses the dataset `snpsVCF` available at our BRGE group package called [`brgedata`](https://github.com/isglobal-brge/brgedata). Vignette can be found [here](https://github.com/isglobal-brge/snpfier/blob/master/vignettes/scoreInvHap.html) 

### Installation

scoreInvHap stable version can be installed from Bioconductor

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("scoreInvHap")

```

or from github repository with the latest changes using devtools package:

```r
# Install devtools
install.packages("devtools")

# Install scoreInvHap package
devtools::install_github("isglobal-brge/scoreInvHap")
```

### Usage

In the folder vignettes, scoreInvHap.Rmd contains an example on how to use scoreInvHap. 


### References
<a id="1">[1]</a> 
Ruiz-Arenas, C., Cáceres, A., López-Sánchez, M., Tolosana, I., Pérez-Jurado, L., & González, J. R. (2019). scoreInvHap: Inversion genotyping for genome-wide association studies. PLoS Genetics, 15(7), [e1008203]. [https://doi.org/10.1371/journal.pgen.1008203](https://doi.org/10.1371/journal.pgen.1008203)
