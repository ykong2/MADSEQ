#####################
## Functions to run models (MCMC) given the data
#  1. function to run the models to fit the data
#' Model to detect and quantify mosaic aneuploidy
#' 
#' Take in the heterozygous sites and coverage information, 
#' use different models (normal, monosomy, mitotic trisomy, meiotic trisomy, 
#' loss of heterozygosity) to fit the data, 
#' and select the model fit the data best according to BIC value and return
#' estimation of the fraction of aneuploid cells.
#' @param hetero A \code{character} specify the location of processed 
#' heterozygous table returned by \code{\link{prepareHetero}} function, or A 
#' \code{GRanges} object returned by \code{\link{prepareHetero}} function
#' @param coverage A \code{character} specify the location of normalized 
#' coverage table returned by \code{\link{normalizeCoverage}} function, or A 
#' \code{GRanges} object from the \code{GRangesList} returned by 
#' \code{\link{normalizeCoverage}} function. Look up your sample by 
#' \code{names(GRangesList)}, and subset your the normalized coverage for your
#' sample by \code{GRangesList["sample_name"]}. For more details, please check
#' the example.
#' @param target_chr A \code{character} specify the chromosome number you want
#' to detect. \bold{Note:} Please check your assembly, use contig name "chr1" 
#' or "1" accordingly.
#' @param adapt A \code{integer} indicate the adaption steps for the MCMC 
#' sampling. Default: 10000
#' @param burnin A \code{integer} indicate burnin steps for the MCMC sampling.
#' Default: 10000. If the posterior distribution is not converged, increasing 
#' burnin steps can be helpful.
#' @param nChain A \code{integer} indicate the number of chains for 
#' the MCMC sampling. Default: 2. \bold{Note:} More than 1 chain is required 
#' if \code{checkConvergence} is set to \code{TRUE}.
#' @param nStep A \code{integer} indicate the number of steps to be recorded 
#' for the MCMC sampling. Default: 10000. Generally, the more steps you record,
#' the more accurate the estimation is. 
#' @param thinSteps A \code{integer} indicate the number of steps to "thin" 
#' (thinSteps=1) means save everystep. Default: 2.
#' @param checkConvergence A \code{Boolean} indicate whether to check the 
#' convergence of independent MCMC chains. If your data is not converged, you
#' may increase adaption step and burnin step. Default: \code{FALSE}
#' @param plot A \code{Boolean}. If \code{TRUE}, the alternative allele 
#' frequency (AAF) for each heterozygous site along the target chromosome will 
#' be plotted.
#' @return An S4 object of class \code{MadSeq} containing the posterior 
#' distribution for the selected model, and deltaBIC between five models.
#' @note 1.If you didn't write normalized coverage into file, please subset the
#' normalized coverage \code{GRanges} object from the \code{GRangesList} 
#' object returned from the \code{\link{normalizeCoverage}} funciton. \cr
#' 2. When specify \code{target_chr}, please make sure it consist with the 
#' contig names in your sequencing data, example: "chr1" and "1". \cr
#' 3. If \code{checkConvergence} set to TRUE, the \code{nChain} has to be >2\cr
#' 4. If it shows that your chains are not converged, 
#' helpful options are increasing the adapt and burnin steps.\cr
#' 5. Because the model is an MCMC sampling process, it can take a very long 
#' time to finish. Running in the background or HPC is recommended.
#' @examples 
#' ## ------------------------------------------------------------------------
#' ## The following example is for the case that normalized coverage and 
#' ## processed heterozygous sites have not been written to files. For more 
#' ## examples, please check the documentation.
#' ## ------------------------------------------------------------------------
#' ##------Prepare Heterozygous Sites
#' ## specify the path to the vcf.gz file for the aneuploidy sample
#' aneuploidy_vcf = system.file("extdata","aneuploidy.vcf.gz",package="MADSEQ")
#' ## specify the path to the bed file containing targeted region
#' target = system.file("extdata","target.bed",package="MADSEQ")
#' ## prepare heterozygous sites
#' aneuploidy_hetero = prepareHetero(aneuploidy_vcf,target, writeToFile=FALSE)
#' 
#' ##------Prepare Normalized Coverage
#' ## specify the path to the bam file
#' aneuploidy_bam = system.file("extdata","aneuploidy.bam",package="MADSEQ")
#' normal_bam = system.file("extdata","normal.bam",package="MADSEQ")
#' 
#' ## prepare coverage data for the samples
#' aneuploidy_cov_gc = prepareCoverageGC(target,aneuploidy_bam,"hg19")
#' normal_cov_gc = prepareCoverageGC(target,normal_bam,"hg19")
#' 
#' ## normalize the coverage
#' normed = normalizeCoverage(aneuploidy_cov_gc,
#'                            control=normal_cov_gc,writeToFile=FALSE)
#'                            
#' ##------subset normalized coverage GRanges object
#' aneuploidy_normed_cov = normed[["aneuploidy_cov_gc"]]
#' ## check chromosome18
#' ## (to speed up the example, we only run one chain and less steps here, 
#' ##  but default settings are recommended in real case)
#' aneuploidy_chr18 = runMadSeq(aneuploidy_hetero, aneuploidy_normed_cov,
#'                              target_chr="chr18", adapt=100, burnin=200,
#'                              nChain =1, nStep = 1000, thinSteps=1)
#' @references Martyn Plummer (2016). rjags: Bayesian Graphical Models using 
#' MCMC. R package version 4-6. \cr
#' \url{https://CRAN.R-project.org/package=rjags}
#' @import rjags
#' @import coda
#' @importFrom VGAM dbetabinom.ab
#' @import GenomicRanges
#' @importFrom stats dnbinom
#' @importFrom utils read.table
#' @importFrom graphics plot abline par mtext
#' @export
#' @seealso \code{\link{MadSeq}}, \code{\link{plotMadSeq}},
#'  \code{\link{plotFraction}}, \code{\link{plotMixture}}
#' @author Yu Kong
runMadSeq = function(
    hetero,
    coverage, 
    target_chr,
    adapt=10000, 
    burnin=10000, 
    nChain=2, 
    nStep=10000, 
    thinSteps=2,
    checkConvergence=FALSE, 
    plot=TRUE){
    
    ## prepare data for MCMC
    ## check if input is the path to file or GRanges object
    if(class(hetero)=="GRanges"){
        sample = as.character(substitute(hetero))
        target_AAF=as.data.frame(hetero[seqnames(hetero)==target_chr])
    }
    else if(class(hetero)=="character"){
        AAF = read.table(hetero,header=TRUE,sep="\t")
        target_AAF = AAF[AAF$seqnames==target_chr,]
    }
    else stop("the file type of 'hetero' is wrong, 
                please provide the path to filtered heterozygous sites 
                or the GRanges object containing filtered heterozygous sites!")
    if(class(coverage)=="GRanges"){
        target_cov=as.data.frame(coverage[seqnames(coverage)==target_chr])
        data_coverage = target_cov$normed_depth
        data_width=sum(width(IRanges(start=target_cov$start,end=target_cov$end)))
        control_coverage = mean(target_cov$ref_depth,na.rm=TRUE)
        auto_cov = coverage[seqnames(coverage)!="chrX"&
                                seqnames(coverage)!="chrY"&
                                seqnames(coverage)!="X"&
                                seqnames(coverage)!="Y"]
        auto_cov_mean = mean(mcols(auto_cov)$ref_depth,na.rm=TRUE)
    }
    else if(class(coverage)=="character"){
        cov = read.table(coverage,header=TRUE,sep="\t")
        target_cov = cov[cov$seqnames==target_chr,]
        data_coverage = target_cov$normed_depth
        data_width=sum(width(IRanges(start=target_cov$start,end=target_cov$end)))
        control_coverage = mean(target_cov$ref_depth,na.rm=TRUE)
        sample = strsplit(coverage,"_")[[1]][1]
        auto_cov = cov[cov$seqnames!="chrX" & cov$seqnames!="chrY" & 
                        cov$seqnames !="X" & cov$seqnames !="Y",]
        auto_cov_mean = mean(auto_cov$ref_depth,na.rm=TRUE)
    }
    
    ## if target chr is sex chromosomes, check the copy of sex chromosome
    if(target_chr=="chrX"|target_chr=="X"){
        chrX_cov_mean = mean(data_coverage)
        ## if there is lack heterozygous sites on chrX (according to poplation
        ## genetics, heterozygosity ~ 0.001, considering filter and genotyping
        ## capacity, if heterozygosity < 0.0001), 
        ## and the coverage for chrX is much smaller than autosome
        ## only one X chromosome
        if(nrow(target_AAF)<(data_width*0.0001) & 
           (chrX_cov_mean/auto_cov_mean)<0.75){
            message("there is only one X chromosome")
            return()
        } 
    }
    else if (target_chr=="chrY"|target_chr=="Y"){
        chrY_cov_mean = mean(data_coverage)
        chrY_cov_len = length(data_coverage)
        ## if there is no heterozygous sites on chrY
        ## and no reads on chrY or coverage for chrY is extremely low
        ## no Y chromosome
        if((nrow(target_AAF)<50&chrY_cov_len<50)|
           (nrow(target_AAF)<50&chrY_cov_mean<(0.2*auto_cov_mean))){
            message("there is NO Y chromosome")
            return()
        }
        else if ((nrow(target_AAF)<(data_width*0.001)&
                (chrY_cov_mean/auto_cov_mean)<0.75)|
                ((nrow(target_AAF)<10)&
                 (chrY_cov_mean/auto_cov_mean)<0.9)){
            message("there is one Y chromosome")
            return()
        }
    }
    
    ## in some cases the coverage data is much longer than AAF data points
    ## so we will downsize the coverage data to speed up, it can also help 
    ## reduce FDR of low confident detection
    if(length(data_coverage)>2*nrow(target_AAF)){
        data_coverage = sample(data_coverage,2*nrow(target_AAF))
    }
    ## print the number of SNP and coverage information
    message (paste("total number of heterozygous site:",nrow(target_AAF)))
    message (paste("total number of coverage",length(data_coverage)))
    
    ## if plot is requested
    if (plot == TRUE){
        par(mfrow=c(1,1))
        par(mar = c(5,4,4,2))
        plotSites(target_AAF)
        mtext(paste("data_coverage: ", round(mean(data_coverage)),
                    ";control_coverage: ", round(control_coverage),sep=""),
                side=1, line=-1)
    }
    
    ## basic settings for MCMC
    load.module("mix")
    
    ## check arch of the platform
#     if(length(grep("i386",R.Version()$arch))>0){
#         message("please check the JAGS model, 
#                 make sure you are using JAGS-4.2.0-Rtools33.exe")
#         return(NULL)
#     }
    ## run normal model
    normal = runNormal(
        target_AAF,
        data_coverage,
        control_coverage,
        checkConvergence=checkConvergence,
        adapt=adapt,
        burnin=burnin,
        nChain=nChain,
        nStep=nStep,
        thinSteps=thinSteps
    )
    
    ## run monosomy model
    monosomy = runMonosomy(
        target_AAF,
        data_coverage,
        control_coverage,
        checkConvergence=checkConvergence,
        adapt=adapt,
        burnin=burnin,
        nChain=nChain,
        nStep=nStep,
        thinSteps=thinSteps
        )
    
    ## run mitotic trisomy model
    mitotic_trisomy = runMitoticTrisomy(
        target_AAF,
        data_coverage,
        control_coverage,
        checkConvergence=checkConvergence,
        adapt=adapt,
        burnin=burnin,
        nChain=nChain,
        nStep=nStep,
        thinSteps=thinSteps
        )
    
    ## run meiotic trisomy model
    meiotic_trisomy = runMeioticTrisomy(
        target_AAF,
        data_coverage,
        control_coverage,
        checkConvergence=checkConvergence,
        adapt=adapt,
        burnin=burnin,
        nChain=nChain,
        nStep=nStep,
        thinSteps=thinSteps
        )
    
    ## run loss of heterozygosity
    LOH = runLOH(
        target_AAF,
        data_coverage,
        control_coverage,
        checkConvergence=checkConvergence,
        adapt=adapt,burnin=burnin,nChain=nChain,
        nStep=nStep,thinSteps=thinSteps
        )
    
    ## model comparison
    # after running all models to fit the data, then compare models by BIC
    message("models done, comparing models")
    BIC = c(normal[[2]], monosomy[[2]], mitotic_trisomy[[2]], 
            meiotic_trisomy[[2]], LOH[[2]])
    # add 10 penalty other comlicated models to get confident results
    BIC = c(BIC[1],BIC[2:5]+10)
    BIC = sort(BIC,decreasing = FALSE)
    best_model = names(which.min(BIC))
    selected = substr(best_model,5,nchar(best_model))
    cat("Order and delta BIC of the preference of models\n")
    delta_BIC = BIC-min(BIC)
    print(delta_BIC)
    cat(paste("model selected:",selected))
    res = get(selected)
    madseq = MadSeq(posterior=res[[1]],deltaBIC=delta_BIC)
    madseq
}