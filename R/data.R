#----------------------Document the Example Data--------------------

#----------------------Aneuploidy Data----------------------------
#' Sequencing data for one chromosome with mosaic aneuploidy
#' 
#' Sequencing data for 200 heterozygous sites of one chormosome with mosaic aneuploidy, and sequencing depth for
#' all the sites from this chromosome, and average depth for the whole genome
#' 
#' @docType data
#' @usage data(aneuploidyChrom)
#' @format A data frame with 200 rows and 4 variables:
#' \describe{
#' \item{chr}{chromosome}
#' \item{position}{locus of the site, in hg19}
#' \item{Ref_D}{read depth for reference allele}
#' \item{Alt_D}{read depth for alternative allele}
#' A numeric vector with 300 values indicating sequencing read depth for all 300 sites.
#' An integer indicates the average sequencing depth for the whole genome.
#' }  
#' @keywords datasets
"aneuploidyChrom"

#' @rdname aneuploidyChrom
"aneuploidy_coverage"

#' @rdname aneuploidyChrom
"genome_coverage"


#----------------------Normal Data----------------------------
#' Sequencing data for a normal chromosome
#' 
#' Sequencing data for 500 heterozygous sites of one normal chormosome, and sequencing depth for all the sites from this chromosome
#' 
#' @docType data
#' @usage data(normalChrom)
#' @format A data frame with 500 rows and 4 variables:
#' \describe{
#' \item{chr}{chromosome}
#' \item{position}{locus of the site, in hg19}
#' \item{Ref_D}{read depth for reference allele}
#' \item{Alt_D}{read depth for alternative allele}
#' A numeric vector with 1060 values indicating sequencing read depth for all 1060 sites.
#' An integer indicates the average sequencing depth for the whole genome.
#' }
#' @keywords datasets
"normalChrom"

#' @rdname normalChrom
"normal_coverage"


#------------------------ Invalid Data----------------------------
#' A dataset does not pass the requirement of MADSEQ package
#' 
#' An example of data that does not meet the requirement of the package.
#' 
#' @docType data
#' @usage data(invalid)
#' @format A data frame with 100 rows and 4 variables:
#' \describe{
#' \item{chr}{chromosome}
#' \item{position}{locus of the site, in hg19}
#' \item{Ref_D}{read depth for reference allele}
#' \item{Alt_D}{read depth for alternative allele}
#' }
#' @keywords datasets
"invalid"


#------------------------ Invalid Data----------------------------
#' A MadSeqOutput from meiotic trisomy model
#' 
#' An example of MadSeqOutput from one of the MADSEQ models
#' 
#' @docType data
#' @usage data(meiotic_trisomy)
#' @format An object of MadSeqOutput class, with two components:
#' \describe{
#' \item{posterior}{posteior disbribution of all parameters}
#' \item{BIC}{BIC value for this model}
#' }
#' @keywords datasets
"meiotic_trisomy"