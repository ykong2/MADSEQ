# MADSEQ
MADSEQ-R package for mosaic aneuploidy detection  

[GitHub page: http://ykong2.github.io/MADSEQ/](http://ykong2.github.io/MADSEQ/)

# Quick Start
## Dependencies
* R (latest version recommended)
* [JAGS](http://mcmc-jags.sourceforge.net/)
* R package: [rjags](https://cran.r-project.org/web/packages/rjags/index.html), or visit their [GitHub](https://github.com/cran/rjags)

## Install 
```{r}
source("https://bioconductor.org/biocLite.R")
biocLite("MADSEQ")
```

###Install from GitHub
```{r}
library(devtools)
install_github("ykong2/MADSEQ", build_vignettes=TRUE)
```
###Install from from [tarball](http://ykong2.github.io/MADSEQ/)
```{r}
install.packages("directory_to_the_downloaded_package/file_name_MADSEQ.tar.gz")
```

## User Guide
Please follow the instruction of the [documentation](http://ykong2.github.io/MADSEQ/documentation.html)

## Other 
For more information, please check our page for [MADSEQ](http://ykong2.github.io/MADSEQ/)
