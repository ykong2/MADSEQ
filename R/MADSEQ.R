#-----------------------DOCUMENTATION PACKAGE--------------------
#' Mosaic Aneuploidy Detection using Massive Parallel Sequencing Data (MADSEQ)
#' 
#' The MADSEQ package provides a group of hierarchical Bayesian models for the
#' detection and quantification of mosaic aneuploidy using massive parallele 
#' sequencing data.  
#' 
#' MADSEQ is a group of hierarchical Bayesian models used for the detection 
#' and quantification of mosaic aneuploidy. The package takes bam file
#' and vcf file as input. There are functions for the calculation of the 
#' coverage for the sequencing data; the normalization of the coverage to 
#' correct GC bias; the detection and quantification of mosaic aneuploidy and
#' the inference of the type of aneuploidy (monosomy, mitotic trisomy, 
#' meiotic trisomy, loss of heterozygosity). The package also includes 
#' function to visualize the estimated distribution for detected mosaic 
#' aneuploidy. To fully understand how to use the MADSEQ package, please check
#' the documentation. The manual explains what data do you need, and how to 
#' process the data to be ready for the model, what steps to follow and how to
#' interpret the output from our model.
#' 
#' @author Yu Kong
#' @docType package
#' @name MADSEQ-package
#' @references Martyn Plummer (2016). rjags: Bayesian Graphical Models using 
#' MCMC. R package version 4-6.
#' \url{http://CRAN.R-project.org/package=rjags} \cr
#' C. Alkan, J. Kidd, T. Marques-Bonet et al (2009).
#' Personalized copy number and segmental duplication maps using 
#' next-generation sequencing. Nature Genetics, 41(10):1061-7.
NULL
#>NULL