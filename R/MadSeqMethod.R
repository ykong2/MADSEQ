## accessor for MadSeq object
## 1. posterior
#' Accessing posterior distribution of MadSeq object
#' 
#' An S4 method to access the posterior distribution of
#' \code{\link{MadSeq}} object
#' @param object A \code{MadSeq} object returned by \code{\link{runMadSeq}} 
#' function
#' @seealso \code{\link{MadSeq}}, \code{\link{runMadSeq}}
#' @rdname posterior
#' @author Yu Kong
#' @export
setMethod("posterior",
          signature = "MadSeq",
          definition = function(object){
            object@posterior
          })

## 2. deltaBIC
#' Accessing delta BIC of MadSeq object
#' 
#' An S4 method to access the delta BIC values of \code{\link{MadSeq}} object
#' @param object A \code{MadSeq} object returned by \code{\link{runMadSeq}} 
#' function
#' @seealso \code{\link{MadSeq}}, \code{\link{runMadSeq}}
#' @rdname deltaBIC
#' @author Yu Kong
#' @export
setMethod("deltaBIC",
          signature = "MadSeq",
          definition = function(object){
            object@deltaBIC
          })


## summary function for MadSeq object
#' Summarize statistics of the MadSeq object
#' 
#' An S4 method to summarize statistics for \code{\link{MadSeq}} object
#' @param object A \code{MadSeq} object returned by \code{\link{runMadSeq}} 
#' function
#' @author Yu Kong
#' @rdname summary-method
#' @aliases summary
#' @export
setMethod("summary",
          signature = "MadSeq",
          definition = function(object){
            model_selected = names(object@deltaBIC)[1]
            model_selected = substr(model_selected,5,nchar(model_selected))
            cat("model selected:",model_selected,"\n")
            summary(object@posterior)
          })

## show method for MadSeq object
setMethod("show",
          signature = "MadSeq",
          definition = function(object){
            model_selected = names(object@deltaBIC)[1]
            model_selected = substr(model_selected,5,nchar(model_selected))
            cat(paste(class(object),
                      "object with the posterior distribution from",
                      model_selected, "model\n",sep=" "))
            print(head(object@posterior))
            cat("------\n")
            print(object@deltaBIC)
          })