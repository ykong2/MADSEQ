#' An S4 class MadSeq object
#' 
#' An MadSeq object returned by the function \code{\link{runMadSeq}}, 
#' the object contains the posterior distribution and deltaBIC value of 
#' a trisomy chromosome 18
#' 
#' @format An MadSeq object
#' @return MadSeq object returned from \code{\link{runMadSeq}} function, 
#' mitotic trisomy has been detected for the chromosome18
#' @docType data
#' @examples 
#' ## to load the data
#' data(aneuploidy_chr18)
#' ## check statistics of the data
#' summary(aneuploidy_chr18)
"aneuploidy_chr18"