## Classes in MADSEQ Package
## 1. MadSeq object which contains result from model
#' The MadSeq class
#'
#' An S4 class contains estimated result returned from \code{\link{runMadSeq}}
#' function
#' @slot posterior A \code{matrix} contains the posterior distribution 
#' from the selected model
#' @slot deltaBIC A \code{numeric vector} contains the deltaBIC value between
#' selected model and other models. The deltaBIC between models indicate the 
#' confidence level that selected model against other models:
#' deltaBIC ~ [0,2]: Not worth more than a bare mention 
#' deltaBIC ~ [2,6]: Positive
#' deltaBIC ~ [6,10]: Strong
#' deltaBIC >10: Very Strong
#' @rdname MadSeq
#' @seealso \code{\link{runMadSeq}}, \code{\link{plotMadSeq}}
#' @section Accessors:
#' \describe{In the code below, \code{x} is a \code{MadSeq} object. \cr\cr
#' \code{posterior(x)}: Get the matrix containing posterior distribution
#' of selected model. \cr\cr
#' \code{deltaBIC(x)}: Get the deltaBIC between selected model and other 
#' models}
#' @section Summary:
#' \describe{In the code below, \code{x} is a \code{MadSeq} object. \cr\cr
#' \code{summary(x)}: summarize the posterior distribution}
#' @section MadSeq Methods:
#' \describe{In the code below, \code{x} is a \code{MadSeq} object. \cr\cr
#' \code{plotMadSeq(x)}: Plot the posterior distribution of all 
#' parameters in selected model. \cr\cr
#' \code{plotFraction(x)}: Plot the estimated distribution of the 
#' fraction of aneuploid sample. \cr\cr
#' \code{plotMixture(x)}: Plot the distribution of AAF estimated from
#' the selected model.}
#' @export
#' @author Yu Kong
MadSeq = setClass( 
    Class="MadSeq",
    representation=representation(
        posterior = "matrix",
        deltaBIC = "numeric"
    ))





