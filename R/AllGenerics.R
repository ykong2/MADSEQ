#### accessor for MadSeq object ####
## 1. posterior
#' @rdname posterior
#' @author Yu Kong
#' @export
setGeneric(
    "posterior", 
    function(object) standardGeneric("posterior"))

## 2. deltaBIC
#' @rdname deltaBIC
#' @export
setGeneric(
    "deltaBIC", 
    function(object) standardGeneric("deltaBIC"))


#### plot functions for MadSeq object ####
## 1. plot all posterior
#' @rdname plotMadSeq
#' @author Yu Kong
#' @export
setGeneric(
    "plotMadSeq", 
    function(object) standardGeneric("plotMadSeq"))

## 2. plot Fraction of aneuploidy cell
#' @rdname plotFraction
#' @author Yu Kong
#' @export
setGeneric(
    "plotFraction", 
    function(object,prob=0.95) standardGeneric("plotFraction"))

## 3. plot Mixtures from posterior distribution
#' @rdname plotMixture
#' @author Yu Kong
#' @export
setGeneric(
    "plotMixture", 
    function(object) standardGeneric("plotMixture"))

