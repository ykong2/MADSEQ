## function to prepare GC content and average coverage for targeted region
#' get sequencing coverage and GC content for targeted regions
#' 
#' Given a bam file and a bed file containing targeted regions, 
#' return sequencing coverage and GC content for each targeted region
#' @param target_bed A \code{character}, specify the path to the location of
#' bed file containing targeted regions.
#' @param bam character, path to the bam file. 
#' Please make sure that bam file is sorted, and the index bam is present
#' @param genome_assembly A \code{character}, indicating the assembly number 
#' of your genome. Default:"hg19". 
#' To see available genome_assembly, 
#' use \code{\link{available.genomes}} from \link{BSgenome} package
#' @return a GRanges object with at least two mcols: depth and GC, 
#' each range indicating a targeted region
#' @note The bam file should be sorted and indexed.
#' @examples 
#' ## specify the path to the location of bed file
#' target = system.file("extdata","target.bed",package="MADSEQ")
#' 
#' ## specify the path to the bam file
#' aneuploidy_bam = system.file("extdata","aneuploidy.bam",package="MADSEQ")
#' normal_bam = system.file("extdata","normal.bam",package="MADSEQ")
#' 
#' ## prepare coverage data for the samples
#' aneuploidy_cov_gc = prepareCoverageGC(target,aneuploidy_bam,"hg19")
#' normal_cov_gc = prepareCoverageGC(target,normal_bam,"hg19")
#' @seealso \code{\link{normalizeCoverage}}
#' @author Yu Kong
#' @import Rsamtools
#' @import GenomicAlignments
#' @import GenomicRanges
#' @import IRanges
#' @import S4Vectors
#' @import Biostrings
#' @import BSgenome
#' @import BSgenome.Hsapiens.UCSC.hg19
#' @import methods
#' @importFrom utils read.table
#' @export
prepareCoverageGC = function(
    target_bed, 
    bam, 
    genome_assembly="hg19"){
    ## prepare coverage
    target_gr = getCoverage(bam=bam, target_bed=target_bed, 
                            genome_assembly = genome_assembly)
    
    ## calculate GC
    gc = calculateGC(range=target_gr, genome_assembly=genome_assembly)
    
    ## check that length of gc equals to length of target
    if(length(gc) == length(target_gr)){
        mcols(target_gr)$GC = gc
    }
    else stop(paste("with", length(target_gr), "target regions, only", 
                    length(gc), 
                    "GC content calculated.Please check your input.",
                    sep=" "))
    
    ## filter out regions around gaps
    target_gr_deGAP = removeGap(target_gr,genome_assembly) 
    
    ## filter out HLA regions
    target_gr_deHLA = removeHLA(target_gr_deGAP,genome_assembly)
    
    ## filter out other highly polymorphic regions
    target_gr_deAQP = removeAQP(target_gr_deHLA,genome_assembly)
    
    ## filter out regions overlap with repeats
    #target_gr_final = removeRE(target_gr_deAQP,genome_assembly)
    target_gr_final = target_gr_deAQP
    target_gr_final
}

##############################################################################
## function to normalize coverage 
#' correct coverage bias due to GC content
#' 
#' function to normalize coverage by GC content and quantile normalization
#' @param object A \code{GRanges} object returned 
#' from \code{\link{prepareCoverageGC}} function. 
#' @param ... additional \code{GRanges} object to pass. \strong{Note1:} If
#' there is only one Granges object given, then coverage will be corrected by 
#' GC content.
#' If there are more than one GRanges object from multiple samples are given, 
#' the function  will first quantile normalize coverage across samples, 
#' then correct coverage by GC content in each sample. \strong{Note2:} If more
#' than one GRanges object provided, make sure they are different samples 
#' sequenced by the same protocol, which means the targeted region is the same
#' \strong{Note3:} If your input samples contain female and male, we suggest
#' you separate them to get a more accurate normalization.
#' @param control A \code{GRanges} object returned from 
#' \code{\link{prepareCoverageGC}} function. \strong{Default value: NULL}. If
#' you have a control normal sample, then put it here
#' @param writeToFile \code{Boolean} Default: \code{TRUE}. If \code{TRUE}, 
#' normalized coverage table for each sample provided will be written to 
#' \code{destination} specified, the file will be named as 
#' "sample_normed_depth.txt". If set to \code{FALSE}, 
#' a \code{\link{GRangesList}} object will be returned
#' @param destination A \code{character}, specify the path to the location 
#' where the normalized coverage table will be written. Default: \code{NULL},
#' the file will be written to current working directory
#' @param plot \code{Boolean} Default: \code{TRUE}. If \code{TRUE}, the 
#' coverage vs. GC content plot before and after normalization will be plotted
#' And the average coverage for each chromosome before and after normalization
#'  will be plotted
#' @return If \code{writeToFile} is set to TRUE, normalized coverage will be 
#' written to the \code{destination}. Otherwise, a \code{\link{GRangesList}}
#' object containing each of input sample will be returned.
#' @note The normalize function works better when you have multiple samples
#' sequenced using the same protocol, namely have the same targeted regions. 
#' And if you have female sample and male sample, the best way is to normalize
#' them separately.
#' @examples 
#' ##------------------------------------
#' ##if you deal with single sample
#' ##------------------------------------
#' ## 1. prepare coverage and gc
#' ## specify the path to the location of bed file
#' target = system.file("extdata","target.bed",package="MADSEQ")
#' 
#' ## specify the path to the bam file
#' aneuploidy_bam = system.file("extdata","aneuploidy.bam",package="MADSEQ")
#' 
#' ## prepare coverage data for the aneuploidy sample
#' aneuploidy_cov_gc = prepareCoverageGC(target,aneuploidy_bam,"hg19")
#' 
#' ## normalize the coverage
#' ##---- if not write to file ----
#' aneuploidy_norm = normalizeCoverage(aneuploidy_cov_gc,writeToFile=FALSE)
#' ## check the GRangesList and subset your sample
#' aneuploidy_norm
#' names(aneuploidy_norm)
#' aneuploidy_norm["aneuploidy_cov_gc"]
#' 
#' ##---- if write to file ----
#' normalizeCoverage(aneuploidy_cov_gc,writeToFile=TRUE,destination=".")
#' 
#' ##-----------------------------------------------------------
#' ##if you deal with multiple samples without normal control
#' ##-----------------------------------------------------------
#' ## specify the path to the location of bed file
#' target = system.file("extdata","target.bed",package="MADSEQ")
#' 
#' ## specify the path to the bam file
#' aneuploidy_bam = system.file("extdata","aneuploidy.bam",package="MADSEQ")
#' normal_bam = system.file("extdata","normal.bam",package="MADSEQ")
#' 
#' ## prepare coverage data for the samples
#' aneuploidy_cov_gc = prepareCoverageGC(target,aneuploidy_bam,"hg19")
#' normal_cov_gc = prepareCoverageGC(target,normal_bam,"hg19")
#' 
#' ## normalize the coverage
#' normed=normalizeCoverage(aneuploidy_cov_gc,normal_cov_gc,writeToFile=FALSE)
#' names(normed)
#' normed["aneuploidy_cov_gc"]
#' normed["normal_cov_gc"]
#' ## or
#' normalizeCoverage(aneuploidy_cov_gc,normal_cov_gc,
#'                   writeToFile=TRUE,destination=".")
#' 
#' ##-----------------------------------------------------------
#' ##if you deal with multiple samples with a normal control
#' ##-----------------------------------------------------------
#' ## specify the path to the location of bed file
#' target = system.file("extdata","target.bed",package="MADSEQ")
#' 
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
#' ## or
#' normalizeCoverage(aneuploidy_cov_gc,control=normal_cov_gc,
#'                   writeToFile=TRUE,destination=".")
#' 
#' @seealso \code{\link{prepareCoverageGC}}
#' @author Yu Kong
#' @references C. Alkan, J. Kidd, T. Marques-Bonet et al (2009).
#' Personalized copy number and segmental duplication maps using 
#' next-generation sequencing. Nature Genetics, 41(10):1061-7.
#' @importFrom preprocessCore normalize.quantiles
#' @importFrom graphics abline axis curve hist legend lines mtext par plot text
#' @importFrom grDevices colors rgb
#' @importFrom stats loess
#' @importFrom utils write.table
#' @export
normalizeCoverage = function(
    object,...,control=NULL,
    writeToFile=TRUE,
    destination=NULL,
    plot=TRUE){
    if(is.null(control)) {
        cat("no control provided")
        data = GRangesList(object, ...)
        control_name = NULL
        sample_name = as.character(substitute(list(object,...)))[-1]
    }
    else {
        data = GRangesList(control,object,...)
        control_name = as.character(substitute(control))
        cat(paste("control:",control_name))
        sample_name=c(control_name,
                        as.character(substitute(list(object,...)))[-1])
    }
    nSample = length(data)
    names(data) = sample_name
    message(paste("there are",nSample,"samples"))
    
    ## check if all the samples have the same number of targeted region
    if (length(unique(sapply(data,length)))>1){
        cat(elementNROWS(data))
        stop("the number of targeted region is different in your samples, 
            please check your input")
    }
    
    ## check if there is >1 samples
    if (nSample > 1) {
        if(plot == TRUE) par(mfrow=c(2,2))
        ## 1. quantile normalize data if there are more than one samples
        quantiled_data = coverageQuantile(data)
        ## 2. correct GC bias for each sample
        corrected = NULL
        for (i in 1:nSample){
            message(paste("correct GC bias in sample'",sample_name[i],"'...",
                            sep=" "))
            sub_corrected = correctGCBias(quantiled_data[[i]],plot=plot)
            if(is.null(corrected)) corrected = GRangesList(sub_corrected)
            else corrected = c(corrected,GRangesList(sub_corrected))
        }
    }
    else if(nSample==1){
        if(plot == TRUE) par(mfrow=c(1,2))
        ## only correct coverage by GC bias
        message(paste("correct GC bias in sample '",sample_name[1],"' ...",
                        sep=""))
        corrected = GRangesList(correctGCBias(data[[1]],plot=plot))
    }
    names(corrected) = sample_name
    
    ## use normed coverage to generate reference coverage
    after_chr = calculateNormedCoverage(corrected,plot=plot)
    ## check how many sample
    if(nSample>1){
        ## check if normal control is provided
        if(is.null(control_name)){
            ## if there are more than 1 sample and there is no control
            ## take median of normed average coverage for each chromosome as 
            ## reference
            ref_cov = apply(after_chr,2,function(x)median(x,na.rm = TRUE))
        }
        ## if control is provided, 
        ## take average coverage for control as reference
        else ref_cov = after_chr[control_name,]
        res = NULL
        for (i in 1:nSample){
            sub_res = corrected[[i]]
            chr_length = sapply(colnames(after_chr),
                                function(x)length(sub_res[seqnames(sub_res)==x]))
            mcols(sub_res)$ref_depth = rep(ref_cov,chr_length)
            if(is.null(res)) res = GRangesList(sub_res)
            else res = c(res,GRangesList(sub_res))
        }
        names(res) = sample_name
    }
    ## if there is only one sample, take the median coverage for chromosome as
    ## reference
    else if(nSample==1){
        res = corrected[[1]]
        mcols(res)$ref_depth = rep(median(after_chr),length(res))
        res = GRangesList(res)
        names(res) = sample_name
    }
    
    ## if write to file requested, then write the normalized coverage table 
    ## file, otherwise return it as a GRangesList 
    if(writeToFile == TRUE){
        ## check if path to write file is provided,
        ## if not write to current working directory
        if(is.null(destination)) path = "."
        else path = destination
        for (i in 1:length(res)){
            tmp_data = as.data.frame(res[[i]])
            write.table(tmp_data,
                        file=paste(path, "/" , sample_name[i],
                                    "_normed_depth.txt", sep=""),
                        quote=FALSE, row.names=FALSE, col.names=TRUE,sep="\t")
            message(paste("normalized depth for sample ",sample_name[i],
                            " is written to ",path, "/" , sample_name[i],
                            "_normed_depth.txt",sep=""))
        }
    }
    else res
}
    
####################### AAF  ####################
#' prepare heterozygous sites for aneuploidy detection
#' 
#' given the vcf file and bed file containing targeted region, generate 
#' processed heterozygous sites for furthur analysis
#' 
#' @param vcffile A \code{character}, specify the path to the location of the 
#' vcf.gz file of your sample. \bold{Note:} the vcf file need to be compressed
#' by \code{bgzip}. The tool is part of \code{tabix} package, 
#' can be download from \url{http://www.htslib.org/}
#' @param target_bed A \code{character}, specify the path to the location of
#' bed file containing targeted regions.
#' @param genome A \code{character}, specify the assembly of your genome. 
#' Default: hg19. To see available genome assembly, 
#' use \code{\link{available.genomes}} from \link{BSgenome} package
#' @param writeToFile \code{Boolean} Default: \code{TRUE}. If \code{TRUE}, 
#' processed table containing heterozygous sites will be written to 
#' \code{destination} specified, the file will be named as 
#' "sample_filtered_heterozygous.txt". If set to \code{FALSE}, 
#' a \code{\link{GRanges}} object containing processed heterozygous sites
#' will be returned
#' @param destination A \code{character}, specify the path to the location 
#' where the processed heterozygous sites table will be written. 
#' Default: \code{NULL}, the file will be written to current working directory
#' @param plot A \code{Boolean} Default: \code{FALSE}. If \code{TRUE},
#' A plot showing AAF before and after filtering for problematic regions will 
#' be generated
#' @return If \code{writeToFile} is set to TRUE, processed table will be 
#' written to the \code{destination}. Otherwise, a \code{\link{GRanges}}
#' object containing each of input sample will be returned.
#' @note 1. The vcf file you provided need to be compressed by bgzip\cr
#' 2. The vcf file should contain depth and allelic depth for variants in the 
#' FORMAT field
#' @examples 
#' ## specify the path to the vcf.gz file for the aneuploidy sample
#' aneuploidy_vcf=system.file("extdata","aneuploidy.vcf.gz",package="MADSEQ")
#' target = system.file("extdata","target.bed",package="MADSEQ")
#' ##------ if not write to file ------
#' aneuploidy_hetero=prepareHetero(aneuploidy_vcf,target,writeToFile=FALSE)
#' 
#' ##------ if write to file ------
#' prepareHetero(aneuploidy_vcf, target,writeToFile=TRUE, destination=".")
#' @seealso \code{\link{runMadSeq}}
#' @author Yu Kong
#' @import VariantAnnotation
#' @import BSgenome
#' @import GenomicRanges
#' @import IRanges
#' @importFrom vcfR read.vcfR
#' @importFrom GenomeInfoDb keepStandardChromosomes
#' @importFrom graphics points
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom utils write.table
#' @importFrom rtracklayer import
#' @importMethodsFrom  GenomeInfoDb seqlevels<-
#' @export
prepareHetero = function(
    vcffile,
    target_bed,
    genome="hg19",
    writeToFile=TRUE,
    destination=NULL,
    plot=FALSE){
    ## first to see if tabix file exist
    try(indexTabix(vcffile,"vcf"))
    
    ## set the target region
    ## read in target bed table
    target_gr = rtracklayer::import(target_bed)
    target_gr = keepStandardChromosomes(target_gr,pruning.mode = "coarse")
    
    ## get sample name
    sample_name = strsplit(vcffile,"/")[[1]]
    sample_name = sample_name[length(sample_name)]
    
    message("reading vcf file")
    vcf = read.vcfR(vcffile)
    
    # get basic pos
    message("processing vcf file")
    basic_factor = c("CHROM","POS","REF","ALT","QUAL","FILTER")
    basic = vcf@fix
    basic = basic[,colnames(basic)%in%basic_factor]
    
    # get GT and AD
    genotype = vcf@gt
    gt_info = genotype[,2]
    gt = sapply(gt_info,function(x)strsplit(x,split=":",fixed=TRUE)[[1]][1])
    ad = sapply(gt_info,function(x)strsplit(x,split=":",fixed=TRUE)[[1]][2])
    ref_d = as.numeric(sapply(ad,function(x)strsplit(x,split=",",fixed=TRUE)[[1]][1]))
    alt_d = as.numeric(sapply(ad,function(x)strsplit(x,split=",",fixed=TRUE)[[1]][2]))
    
    res = GRanges(seqnames=basic[,"CHROM"],
                  ranges=IRanges(start=as.numeric(basic[,"POS"]),end=as.numeric(basic[,"POS"])))
    mcols = data.frame(basic[,-c(1,2)],GT=gt,Ref_D=ref_d,Alt_D=alt_d,DP=ref_d+alt_d)
    mcols(res) = mcols
    
    ## keep only SNPs
    res = res[mcols(res)$GT=="0/1"|mcols(res)$GT=="1/0"]
    
    res = res[!is.na(mcols(res)$Alt_D)&
                  !is.na(mcols(res)$Ref_D)&
                  !is.na(mcols(res)$REF)&
                  !is.na(mcols(res)$ALT)]
    names(res) = seq(1:length(res))
    
    ## keep only SNPs in target regions
    ov = findOverlaps(res,target_gr)
    res = unique(res[queryHits(ov)])
    mcols(res)$ALT = gsub("\\,<NON_REF>","",mcols(res)$ALT)
    ## filter:
    message("filtering vcf file")
    ## 1. filter out non standard chromosome
    res = keepStandardChromosomes(res,pruning.mode = "coarse")
    
    ## 2.filter low DP/Ref_D/Alt_D by 5% quantile
    DP_cut = quantile(mcols(res)$DP,0.05)
    Ref_D_cut = quantile(mcols(res)$Ref_D,0.05)
    Alt_D_cut = quantile(mcols(res)$Alt_D,0.05)
    res = res[mcols(res)$DP>DP_cut&mcols(res)$Ref_D>Ref_D_cut&mcols(res)$Alt_D>Alt_D_cut]
    
    ## 3. remove SNPs around gap
    res = removeGap(res,genome)
    
    ## 4. filter out HLA region on chr6 due to the polymorphics of HLA, the calling of heterozygous sites are ambiguous most of the time
    res = removeHLA(res,genome)
    
    ## 5. the AQP7 region on chr9 also has ambigous heterozygous calling because of the polymorphics, we will remove them for downstream analysis
    res = removeAQP(res,genome)
    
    ## 6. further treat regions most of the SNP are biased due to sequencing issues
    res = filter_hetero(res,binsize=10,plot=plot)
    
    if(writeToFile == TRUE){
        ## check if path to write file is provided,
        ## if not write to current working directory
        if(is.null(destination)) path = "."
        else path = destination
        data = as.data.frame(res)
        write.table(data,
                    file=paste(path, "/" , sample_name,
                               "_filtered_heterozygous.txt", sep=""),
                    quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
        message(paste("filtered heterozygous sites for sample ",sample_name,
                      " is written to ",path, "/" , sample_name,
                      "_filtered_heterozygous.txt",sep=""))
    }
    else
        res
}
