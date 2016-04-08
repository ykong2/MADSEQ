## 1. plot all posterior
#' density plot for posterior distribution of selected model
#' 
#' plot the density plot for each of the parameters in the posterior
#' distribution from selected model
#' @param object A \code{\link{MadSeq}} object returned by 
#' \code{\link{runMadSeq}} function. 
#' @rdname plotMadSeq
#' @return the density plot for parameters in the posterior distribution of 
#' selected model.
#' @examples 
#' ## load the example MadSeq object come with the package
#' data("aneuploidy_chr18")
#' 
#' ## plot the posterior distribution
#' plotMadSeq(aneuploidy_chr18)
#' @seealso \code{\link{runMadSeq}}, \code{\link{plotFraction}}, 
#' \code{\link{plotMixture}}
#' @author Yu Kong
#' @importFrom graphics plot hist par
#' @importFrom stats density
#' @export
setMethod(
    "plotMadSeq",
    signature = "MadSeq",
    definition = function(object){
        data = as.data.frame(object@posterior)
        n_total = dim(data)[2]
        all_name = names(data)
        par(mfrow=c(ceiling(n_total/4),4))
        par(mar=c(2.5,4,3,2))
        for (i in all_name){
            tmp_dat = data[,i]
            if (min(tmp_dat)==max(tmp_dat)){
                hist(tmp_dat, prob=TRUE,
                    main=paste("posterior",i),xlab="")
            }
            else{
                plot(stats::density(tmp_dat), xlab="", lwd=2, col="blue",
                    main=paste("posterior",i))
            }
        }
    })



## 2. plot Fraction of aneuploidy cell
#' histgram for the fraction of aneuploid cells estimated by MadSeq model
#' 
#' histgram of the posterior distribution of the fraction of aneuploid cells 
#' estimated by the selected model.
#' @param object A \code{\link{MadSeq}} object returned by 
#' \code{\link{runMadSeq}} function. 
#' @param prob A \code{numeric} value between 0~1 specify the highest posterior
#' interval (similar to credible interval) for the distribution. Default: 0.95.
#' @return the histgram of posterior distribution of the fraction
#' @examples 
#' ## load the example MadSeq object come with the package
#' data("aneuploidy_chr18")
#' 
#' ## plot estimated fraction of aneuploid cells
#' plotFraction(aneuploidy_chr18)
#' @rdname plotFraction
#' @note If normal model has been selected by \code{\link{runMadSeq}} function,
#' no fraction plot will be produced by this function.
#' @seealso \code{\link{runMadSeq}}, \code{\link{plotMadSeq}}, 
#' \code{\link{plotMixture}}
#' @author Yu Kong
#' @export
#' @importFrom graphics lines hist text par
#' @importFrom grDevices rgb
#' @importFrom stats density
setMethod(
    "plotFraction",
    signature = "MadSeq",
    definition = function(object,prob=0.95){
        if(is.element("f",colnames(object@posterior))){
            par(mfrow=c(1,1))
            f = object@posterior[,"f"]
            f_density = stats::density(f)
            y_max = max(hist(f,breaks=40,plot=FALSE)$density)
            f_mean = signif(mean(f),2)
            f_HDI = signif(HDIofMCMC(f,credMass = prob),2)
            hist(f, prob=TRUE, breaks=50,
                xlab="", ylab="", main="posterior estimate of the fraction",
                col=rgb(0,1,0,0.2), mgp=c(2,1,0), ylim=c(0,1.1*y_max), 
                xlim=c(min(f),max(f)))
            lines(f_density,col="black",lty=2,lwd=3)
            lines(x=c(f_mean,f_mean),y=c(0,1.01*y_max),col="red",lwd=3)
            text((max(f)+min(f))/2,1.03*y_max,
                paste("mean =",as.character(f_mean)),cex=1.1)
            lines(x=c(f_HDI[1],f_HDI[1]),y=c(0,0.25*y_max),col="red",lwd=3)
            lines(x=c(f_HDI[2],f_HDI[2]),y=c(0,0.25*y_max),col="red",lwd=3)
            text((max(f)+min(f))/2, 1.09*y_max,
                paste(as.character(100*prob),
                        "% HPD = [",as.character(f_HDI[1]),
                        ",",as.character(f_HDI[2]),"]",sep=""),
                cex=1)
        }
        else
            message("selected model is normal, f is estimated to be 0%")
    })

## 3. plot Mixtures from posterior distribution
#' density plot for the posterior distribution of alternative allele frequency
#' estimated from the selected model
#' 
#' density plot presents the posterior distribution of alternative allele
#' frequency (AAF) estimated from selected model 
#' @param object A \code{\link{MadSeq}} object returned by 
#' \code{\link{runMadSeq}} function. 
#' @return density plot for the posterior distribution of AAF
#' @examples 
#' ## load the example MadSeq object come with the package
#' data("aneuploidy_chr18")
#' 
#' ## plot the distribution of estimated AAF
#' plotMixture(aneuploidy_chr18)
#' @rdname plotMixture
#' @seealso \code{\link{runMadSeq}}, \code{\link{plotMadSeq}}, 
#' \code{\link{plotFraction}}
#' @author Yu Kong
#' @importFrom graphics curve text par
#' @importFrom stats dbeta
#' @export
setMethod(
    "plotMixture",
    signature = "MadSeq",
    definition = function(object){
        par(mfrow=c(1,1))
        model_selected = names(object@deltaBIC)[1]
        model_selected = substr(model_selected,5,nchar(model_selected))
        data = object@posterior
        data_variables = colnames(data)
        data_median = apply(data,2,median)
        ## set x for futher plot
        x = seq(0,1,0.01)
        ## check if the output is from one mixture model
        if(model_selected == "normal"){
            a = data_median["mu"]*data_median["kappa"]
            b = (1-data_median["mu"])*data_median["kappa"]
            p1 = 0.99
            p2 = 0.01
            y_max = max(dbeta(seq(0,1,0.01),a,b))
            mu = signif(data_median["mu"],2)
            curve(p1*dbeta(x,a,b), col="blue", lwd=2,
                    xlab="alternative allele frequency", ylab="density",
                    main="posterior distribution of mixtures (normal)", 
                    mgp=c(2,1,0), ylim=c(0,1.1*y_max), bty="l",cex.main=0.8)
            curve(p2*dbeta(x,1,1),col="grey4",lwd=2,lty=2,add=TRUE)
            text(x=mu,1.05*y_max,paste("m =",as.character(mu)))
        }
        ## check if the output is from meiotic trisomy mixture model
        else if(model_selected == "meiotic_trisomy"){
            a1 = data_median["mu[1]"]*data_median["kappa"]
            b1 = (1-data_median["mu[1]"])*data_median["kappa"]
            a2 = data_median["mu[2]"]*data_median["kappa"]
            b2 = (1-data_median["mu[2]"])*data_median["kappa"]
            a3 = data_median["mu[3]"]*data_median["kappa"]
            b3 = (1-data_median["mu[3]"])*data_median["kappa"]
            a4 = data_median["mu[4]"]*data_median["kappa"]
            b4 = (1-data_median["mu[4]"])*data_median["kappa"]
            a5 = 1
            b5 = 1
            p1 = data_median["p[1]"]
            p2 = data_median["p[2]"]
            p3 = data_median["p[3]"]
            p4 = data_median["p[4]"]
            p5 = 0.01
            m1 = signif(data_median["mu[1]"],2)
            m2 = signif(data_median["mu[2]"],2)
            m3 = signif(data_median["mu[3]"],2)
            m4 = signif(data_median["mu[4]"],2)
            y_max = max(c(p1*dbeta(seq(0,1,0.01),a1,b1),
                            p2*dbeta(seq(0,1,0.01),a2,b2),
                            p3*dbeta(seq(0,1,0.01),a3,b3),
                            p4*dbeta(seq(0,1,0.01),a4,b4),
                            p5*dbeta(seq(0,1,0.01),a5,b5)))
            curve(p1*dbeta(x,a1,b1), col="blue", lwd=2,
                    xlab="alternative allele frequency", ylab="density",
                    main="posterior distribution of mixtures(meiotic trisomy)",
                    mgp=c(2,1,0), xlim=c(0,1), bty="l", ylim=c(0,1.25*y_max),
                    cex.main=0.8)
            curve(p2*dbeta(x,a2,b2),col="blue",lwd=2,add=TRUE)
            curve(p3*dbeta(x,a3,b3),col="green4",lwd=2,add=TRUE)
            curve(p4*dbeta(x,a4,b4),col="green4",lwd=2,add=TRUE)
            curve(p5*dbeta(x,a5,b5),col="grey4",lwd=2,add=TRUE,lty=2)
            text(m1,1.05*max(p1*dbeta(seq(0,1,0.01),a1,b1)),
                paste("m1 =",as.character(m1)))
            text(m2,1.05*max(p2*dbeta(seq(0,1,0.01),a2,b2)),
                paste("m2 =",as.character(m2)))
            text(m3,1.25*max(p3*dbeta(seq(0,1,0.01),a3,b3)),
                paste("m3 =",as.character(m3)))
            text(m4,1.25*max(p4*dbeta(seq(0,1,0.01),a4,b4)),
                paste("m4 =",as.character(m4)))
        }
        ## check if the output from LOH mixture model
        else if(model_selected == "LOH"){
            a1 = data_median["mu[1]"]*data_median["kappa"]
            b1 = (1-data_median["mu[1]"])*data_median["kappa"]
            a2 = data_median["mu[2]"]*data_median["kappa"]
            b2 = (1-data_median["mu[2]"])*data_median["kappa"]
            a3 = data_median["mu[3]"]*data_median["kappa"]
            b3 = (1-data_median["mu[3]"])*data_median["kappa"]
            a4 = 1
            b4 = 1
            p1 = data_median["p[1]"]
            p2 = data_median["p[2]"]
            p3 = data_median["p[3]"]
            p4 = 0.01
            m1 = signif(data_median["mu[1]"],2)
            m2 = signif(data_median["mu[2]"],2)
            m3 = signif(data_median["mu[3]"],2)
            y_max = max(c(p1*dbeta(seq(0,1,0.01),a1,b1),
                            p2*dbeta(seq(0,1,0.01),a2,b2),
                            p3*dbeta(seq(0,1,0.01),a3,b3),
                            p4*dbeta(seq(0,1,0.01),a4,b4)))
            curve(p1*dbeta(x,a1,b1), col="green4", lwd=2,
                    xlab="alternative allele frequency", ylab="density",
                    main="posterior distribution of mixtures (LOH)",
                    mgp=c(2,1,0), xlim=c(0,1), bty="l", ylim=c(0,1.25*y_max),
                    cex.main=0.8)
            curve(p2*dbeta(x,a2,b2),col="green4",lwd=2,add=TRUE)
            curve(p3*dbeta(x,a3,b3),col="blue",lwd=2,add=TRUE)
            curve(p4*dbeta(x,a4,b4),col="grey4",lwd=2,add=TRUE,lty=2)
            text(m1,1.05*max(p1*dbeta(seq(0,1,0.01),a1,b1)),
                paste("m1 =",as.character(m1)))
            text(m2,1.05*max(p2*dbeta(seq(0,1,0.01),a2,b2)),
                paste("m2 =",as.character(m2)))
            text(m3,1.25*max(p3*dbeta(seq(0,1,0.01),a3,b3)),
                paste("m3 =",as.character(m3)))
        }
        # if not all of the above, then plot two mixtures
        else{
            a1 = data_median["mu[1]"]*data_median["kappa"]
            b1 = (1-data_median["mu[1]"])*data_median["kappa"]
            a2 = data_median["mu[2]"]*data_median["kappa"]
            b2 = (1-data_median["mu[2]"])*data_median["kappa"]
            a3 = 1
            b3 = 1
            p1 = 0.495
            p2 = 0.495
            p3 = 0.01
            m1 = signif(data_median["mu[1]"],2)
            m2 = signif(data_median["mu[2]"],2)
            y_max = max(c(p1*dbeta(seq(0,1,0.01),a1,b1),
                            p2*dbeta(seq(0,1,0.01),a2,b2)))
            curve(p1*dbeta(x,a1,b1), col="blue", lwd=2,
                    xlab="alternative allele frequency", ylab="density",
                    main=paste("posterior distribution of mixtures (",
                            model_selected,")", sep=""),cex.main=0.8,
                    mgp=c(2,1,0), xlim=c(0,1), bty="l", ylim=c(0,1.1*y_max))
            curve(p2*dbeta(x,a2,b2),col="blue",lwd=2,add=TRUE)
            curve(p3*dbeta(x,a3,b3),col="grey4",lwd=2,add=TRUE,lty=2)
            text(m1+0.1,y_max,paste("m1 =",as.character(m1)))
            text(m2-0.1,y_max,paste("m2 =",as.character(m2)))
        }
    })