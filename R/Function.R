#-----------------1. Functions to prepare dataList for the analysis----------------

# 1. function to check if the data is valid for further analysis

# data considered as valid if all the sites:
# 1). come from the same chromosome 2). all of them are biallelic sites (with both ref/alt alleles) 3). it has the four fields: chr, position, Ref_D,Alt_D
#' Check if the data is valid for MAD-SEQ analysis
#' 
#' Take in the data frame and check if it meets all the requirements
#' @param data The data frame to be used for analysis
#' @return Ture or error code indicate why the data fails the requirements
#' @export
validateData = function(data){
  # 1. check four required columns
  # check if there is chr column
  if(is.na(match("chr",colnames(data)))){
    stop("Error: No field named as 'chr', please check your data") 
  }
  # check if there is position column
  if(is.na(match("position",colnames(data)))){
    stop("Error: No field named as 'position', please check your data") 
  }
  # check if there is Ref_D column
  if(is.na(match("Ref_D",colnames(data)))){
    stop("Error: No field named as 'Ref_D', please check your data") 
  }
  # check if there is Alt_D column
  if(is.na(match("Alt_D",colnames(data)))){
    stop("Error: No field named as 'Alt_D', please check your data") 
  }
  
  # 2. check if all the sites come from the same chromosome
  if(length(table(data$chr))>1) {
    stop("there are more than one chromosomes, please check your data")
  }
  
  # 3. check if all the sites are biallelic sites
  for (i in 1:nrow(data)){
    if(data[i,"Alt_D"]/(data[i,"Alt_D"]+data[i,"Ref_D"])==0|data[i,"Alt_D"]/(data[i,"Alt_D"]+data[i,"Ref_D"])==1){
      stop(paste("the site in the",i,"row is not biallelic, please check your data"))
    }
  }  
  return(TRUE)
}

#--------------------2. Functions to run MCMC on the data----------------
#  1. function to run the One Mixture Model
#' Use the one mixture model to fit the data
#' 
#' Take in the heterozygous sites and coverage information, use one mixture model to fit the data
#' @param data dataframe: include at least four fields: chr (chromosome), position (position), Ref_D (read depth for reference allele) and Alt_D (read depth for alternative allele). Check if your input dataframe meets all the requirements by function validateData(data) 
#' @param data_coverage numeric vector: the coverage information for all the sites on this chromosome (sites with all genotypes)
#' @param control_coverage numeric vector: control coverage information: it can be the average coverage of the whole genome. Or the mean of the normalized coverage.
#' @param adapt integer: the number of steps to "tune" the samplers. (default=10000)
#' @param burnin integer: the number of steps to "burn-in" the samplers. (default=10000)
#' @param nChain integer: the number of chains to run. (default=2)
#' @param nStep integer: total number of steps in chains to save (default=20000)
#' @param thinSteps interger: save data every thinSteps (default=2)
#' @param checkConvergence Boolean: check if the posterior estimation is converged by gelman and rubin diagnostic. Only available when you have more than 2 chains. (default=FALSE) If chains are not converged, then you should increase adapt and burnin steps.
#' @param plot Boolean: plot the alternative allele frequencies along the chromosome. (default=TRUE)
#' @return A MadSeqOutput object with two components
#' @return \code{posterior} the posterior output from MCMC 
#' @return \code{BIC} the BIC value for this model
#' @import rjags
#' @import coda
#' @export
#' @seealso \code{\link{runModelMonosomy}} \code{\link{runModelTrisomy}} \code{\link{runModelFour}} \code{\link{validateData}}
runModelOne = function(data,data_coverage,control_coverage,adapt=10000,burnin=10000,nChain=2,nStep=20000,thinSteps=2,checkConvergence=FALSE,plot=TRUE){
  # basic settings
  adaptSteps = adapt              # Number of steps to "tune" the samplers.
  burnInSteps = burnin            # Number of steps to "burn-in" the samplers.
  nChains = nChain               # Number of chains to run.
  numSavedSteps=nStep          # Total number of steps in chains to save.
  thinSteps=thinSteps                  # Number of steps to "thin" (1=keep every step).
  nIter = ceiling( ( numSavedSteps * thinSteps ) / nChains )
  parameters_one = c("kappa","p_cov","r_cov","theta","mu","p1","p2")
  record_one = c("kappa","p_cov","r_cov","mu","p1","p2")
  
  # generate dataList for the MCMC
  dataList = creatDataList_One(data=data,N1=control_coverage,N2=data_coverage,plot=plot)
  
  # run model
  fpath = system.file("rjags","OneMix.txt",package="MADSEQ")
  jagsModel_one = jags.model( fpath, data=dataList ,n.chains=nChains , n.adapt=adaptSteps )
  update( jagsModel_one , n.iter=burnInSteps )
  codaSamples_one = coda.samples( jagsModel_one , variable.names=parameters_one , n.iter=nIter , thin=thinSteps )
  mat = as.matrix(codaSamples_one)
  theta = theta_median(mat)
  mat = mat[,record_one]
  
  if(checkConvergence==TRUE){
    diag = gelman.diag(codaSamples_one,multivariate = FALSE)
    upper = max(diag$psrf[,2],na.rm=TRUE)
    if (upper<=2) print("converged!")
    else if (upper > 2 & upper <= 10) print("not fulled converged")
    else if (upper>10) print("not converge")
  }
  post = apply(mat,2,median)
  BIC = calculateBIC1(post,theta,dataList)
  names(BIC) = "BIC_one"
  res = list(posterior=mat,BIC=BIC)
  res = as.MadSeqOutput(res)
  return(res)
}


# 2. function to run Monosomy Mixture Model
#' Use the two mixture monosomy model to fit the data
#' 
#' Take in the heterozygous sites and coverage information, use two mixture monosomy model to fit the data
#' @param data dataframe: include at least four fields: chr (chromosome), position (position), Ref_D (read depth for reference allele) and Alt_D (read depth for alternative allele). Check if your input dataframe meets all the requirements by function validateData(data) 
#' @param data_coverage numeric vector: the coverage information for all the sites on this chromosome (sites with all genotypes)
#' @param control_coverage numeric vector: control coverage information: it can be the average coverage of the whole genome. Or the mean of the normalized coverage.
#' @param adapt integer: the number of steps to "tune" the samplers. (default=10000)
#' @param burnin integer: the number of steps to "burn-in" the samplers. (default=10000)
#' @param nChain integer: the number of chains to run. (default=2)
#' @param nStep integer: total number of steps in chains to save (default=20000)
#' @param thinSteps interger: save data every thinSteps (default=2)
#' @param checkConvergence Boolean: check if the posterior estimation is converged. Only available when you have more than 2 chains. (default=FALSE)
#' @param plot Boolean: plot the alternative allele frequencies along the chromosome. (default=TRUE)
#' @return A MadSeqOutput object with two components
#' @return \code{posterior} the posterior output from MCMC 
#' @return \code{BIC} the BIC value for this model
#' @import rjags
#' @import coda
#' @export
#' @seealso \code{\link{runModelOne}} \code{\link{runModelTrisomy}} \code{\link{runModelFour}} \code{\link{validateData}}
runModelMonosomy = function(data,data_coverage,control_coverage,adapt=10000,burnin=10000,nChain=2,nStep=20000,thinSteps=2,checkConvergence=FALSE,plot=TRUE){
  # basic settings
  adaptSteps = adapt              # Number of steps to "tune" the samplers.
  burnInSteps = burnin            # Number of steps to "burn-in" the samplers.
  nChains = nChain               # Number of chains to run.
  numSavedSteps=nStep          # Total number of steps in chains to save.
  thinSteps=thinSteps                  # Number of steps to "thin" (1=keep every step).
  nIter = ceiling( ( numSavedSteps * thinSteps ) / nChains )
  parameters_two = c("p_cov","m_cov","r_cov","mu","kappa","d1","d2","f","theta","p1","p2","p3")
  record_two = c("f","kappa","mu[1]","mu[2]","d1","d2","p_cov","r_cov","m_cov","p1","p2","p3")
  
  # generate dataList for the MCMC
  dataList = creatDataList_Two(data=data,N1=control_coverage,N2=data_coverage,plot=plot)
  
  # run model
  fpath = system.file("rjags","TwoMixMonosomy.txt",package="MADSEQ")
  jagsModel3 = jags.model( fpath , data=dataList , n.chains=nChains , n.adapt=adaptSteps )
  update( jagsModel3 , n.iter=burnInSteps )
  codaSamples3 = coda.samples( jagsModel3 , variable.names=parameters_two , n.iter=nIter , thin=thinSteps )
  mat = as.matrix(codaSamples3)
  theta = theta_median(mat)
  mat = mat[,record_two]
  
  if(checkConvergence==TRUE){
    diag = gelman.diag(codaSamples3,multivariate = FALSE)
    upper = max(diag$psrf[,2],na.rm=TRUE)
    if (upper<=2) print("converged!")
    else if (upper > 2 & upper <= 10) print("not fulled converged")
    else if (upper>10) print("not converge")
  }
  
  post = apply(mat,2,median)
  BIC = calculateBIC2(post,theta,dataList)
  names(BIC) = "BIC_monosomy"
  res = list(posterior=mat,BIC=BIC)
  res = as.MadSeqOutput(res)
  return(res)
}


# 3. function to run Trisomy Two Mixture Model
#' Use the two mixture trisomy model to fit the data
#' 
#' Take in the heterozygous sites and coverage information, use two mixture trisomy model to fit the data
#' @param data dataframe: include at least four fields: chr (chromosome), position (position), Ref_D (read depth for reference allele) and Alt_D (read depth for alternative allele). Check if your input dataframe meets all the requirements by function validateData(data) 
#' @param data_coverage numeric vector: the coverage information for all the sites on this chromosome (sites with all genotypes)
#' @param control_coverage numeric vector: control coverage information: it can be the average coverage of the whole genome. Or the mean of the normalized coverage.
#' @param adapt integer: the number of steps to "tune" the samplers. (default=10000)
#' @param burnin integer: the number of steps to "burn-in" the samplers. (default=10000)
#' @param nChain integer: the number of chains to run. (default=2)
#' @param nStep integer: total number of steps in chains to save (default=20000)
#' @param thinSteps interger: save data every thinSteps (default=2)
#' @param checkConvergence Boolean: check if the posterior estimation is converged. Only available when you have more than 2 chains. (default=FALSE)
#' @param plot Boolean: plot the alternative allele frequencies along the chromosome. (default=TRUE)
#' @return A MadSeqOutput object with two components
#' @return \code{posterior} the posterior output from MCMC 
#' @return \code{BIC} the BIC value for this model
#' @import rjags
#' @import coda
#' @export
#' @seealso \code{\link{runModelOne}} \code{\link{runModelMonosomy}} \code{\link{runModelFour}} \code{\link{validateData}}
runModelTrisomy = function(data,data_coverage,control_coverage,adapt=10000,burnin=10000,nChain=2,nStep=20000,thinSteps=2,checkConvergence=FALSE,plot=TRUE){
  # basic settings
  adaptSteps = adapt              # Number of steps to "tune" the samplers.
  burnInSteps = burnin            # Number of steps to "burn-in" the samplers.
  nChains = nChain               # Number of chains to run.
  numSavedSteps=nStep          # Total number of steps in chains to save.
  thinSteps=thinSteps                  # Number of steps to "thin" (1=keep every step).
  nIter = ceiling( ( numSavedSteps * thinSteps ) / nChains )
  parameters_two = c("p_cov","m_cov","r_cov","mu","kappa","d1","d2","f","theta","p1","p2","p3")
  record_two = c("f","kappa","mu[1]","mu[2]","d1","d2","p_cov","r_cov","m_cov","p1","p2","p3")
  
  # generate dataList for the MCMC
  dataList = creatDataList_Two(data=data,N1=control_coverage,N2=data_coverage,plot=plot)
  
  # run model
  fpath = system.file("rjags","TwoMixTrisomy.txt",package="MADSEQ")
  jagsModel2 = jags.model(fpath , data=dataList , n.chains=nChains , n.adapt=adaptSteps )
  update( jagsModel2 , n.iter=burnInSteps )
  codaSamples2 = coda.samples( jagsModel2 , variable.names=parameters_two , n.iter=nIter , thin=thinSteps )
  mat = as.matrix(codaSamples2)
  theta = theta_median(mat)
  mat = mat[,record_two]
  
  if(checkConvergence==TRUE){
    diag = gelman.diag(codaSamples2,multivariate = FALSE)
    upper = max(diag$psrf[,2],na.rm=TRUE)
    if (upper<=2) print("converged!")
    else if (upper > 2 & upper <= 10) print("not fulled converged")
    else if (upper>10) print("not converge")
  }
  
  post = apply(mat,2,median)
  BIC = calculateBIC2(post,theta,dataList)
  names(BIC) = "BIC_trisomy"
  res = list(posterior=mat,BIC=BIC)
  res = as.MadSeqOutput(res)
  return(res)
}


# 4. function to run Four Mixture Model
#' Use the four mixture model to fit the data
#' 
#' Take in the heterozygous sites and coverage information, use four mixture model to fit the data
#' @param data dataframe: include at least four fields: chr (chromosome), position (position), Ref_D (read depth for reference allele) and Alt_D (read depth for alternative allele). Check if your input dataframe meets all the requirements by function validateData(data) 
#' @param data_coverage numeric vector: the coverage information for all the sites on this chromosome (sites with all genotypes)
#' @param control_coverage numeric vector: control coverage information: it can be the average coverage of the whole genome. Or the mean of the normalized coverage.
#' @param adapt integer: the number of steps to "tune" the samplers. (default=10000)
#' @param burnin integer: the number of steps to "burn-in" the samplers. (default=10000)
#' @param nChain integer: the number of chains to run. (default=2)
#' @param nStep integer: total number of steps in chains to save (default=20000)
#' @param thinSteps interger: save data every thinSteps (default=2)
#' @param checkConvergence Boolean: check if the posterior estimation is converged. Only available when you have more than 2 chains. (default=FALSE)
#' @param plot Boolean: plot the alternative allele frequencies along the chromosome. (default=TRUE)
#' @return A MadSeqOutput object with two components
#' @return \code{posterior} the posterior output from MCMC 
#' @return \code{BIC} the BIC value for this model
#' @import rjags
#' @import coda
#' @export
#' @seealso \code{\link{runModelOne}} \code{\link{runModelMonosomy}} \code{\link{runModelTrisomy}} \code{\link{validateData}}
runModelFour = function(data,data_coverage,control_coverage,adapt=10000,burnin=10000,nChain=2,nStep=20000,thinSteps=2,checkConvergence=FALSE,plot=TRUE){
  # basic settings
  adaptSteps = adapt              # Number of steps to "tune" the samplers.
  burnInSteps = burnin            # Number of steps to "burn-in" the samplers.
  nChains = nChain               # Number of chains to run.
  numSavedSteps=nStep          # Total number of steps in chains to save.
  thinSteps=thinSteps                  # Number of steps to "thin" (1=keep every step).
  nIter = ceiling( ( numSavedSteps * thinSteps ) / nChains )
  parameters_four = c("p_cov","m_cov","r_cov","mu","kappa","d1","d2","d3","d4","p1","p2","p3","p4","p5","f","theta")
  record_four = c("f","p1","p2","p3","p4","p5","kappa","mu[1]","mu[2]","mu[3]","mu[4]","d1","d2","d3","d4","p_cov","r_cov","m_cov")
  
  # generate dataList for the MCMC
  dataList = creatDataList_Four(data=data,N1=control_coverage,N2=data_coverage,plot=plot)
  
  # run model
  fpath = system.file("rjags","FourMix.txt",package="MADSEQ")
  jagsModel_four = jags.model( fpath , data=dataList ,n.chains=nChains , n.adapt=adaptSteps ) 
  update(jagsModel_four , n.iter=burnInSteps)
  codaSamples_four = coda.samples( jagsModel_four , variable.names=parameters_four ,  n.iter=nIter , thin=thinSteps )
  mat = as.matrix(codaSamples_four)
  theta = theta_median(mat)
  mat = mat[,record_four]
  
  if(checkConvergence==TRUE){
    diag = gelman.diag(codaSamples_four,multivariate = FALSE)
    upper = max(diag$psrf[,2],na.rm=TRUE)
    if (upper<=2) print("converged!")
    else if (upper > 2 & upper <= 10) print("not fulled converged")
    else if (upper>10) print("not converge")
  }
  
  post = apply(mat,2,median)
  BIC = calculateBIC4(post,theta,dataList)
  names(BIC) = "BIC_four"
  res = list(posterior=mat,BIC=BIC)
  res = as.MadSeqOutput(res)
  return(res)
}


# 5. function to compare models, and select the best model
#' select the best MADSEQ model for input data
#' 
#' compare multiple MadSeq models using BIC. Print out the order of preference
#' @param object1 MadSeqOutput object, which is produced by one of the runModel functions
#' @param object2 MadSeqOutput object
#' @param ... MadSeqOutput object
#' @return a character indicating which is the best model
#' @details the function compare different models using BIC, the model with lowest BIC value is the best model. The deltaBIC between models indicate the evidence against higher BIC
#' deltaBIC ~ [0,2]: Not worth more than a bare mention                                                                                                                                                    
#' deltaBIC ~ [2,6]: Positive
#' deltaBIC ~ [6,10]: Strong
#' deltaBIC >10: Very Strong
#' @export
modelSelection = function(object1,object2,...){
  data = list(object1,object2,...)
  n_items = length(data)
  type = NULL
  BIC = NULL
  for (i in 1:n_items){
    tmp = data[[i]]
    tmp_BIC = tmp$BIC
    tmp_type = substr(names(tmp_BIC),5,nchar(names(tmp_BIC)))
    BIC = c(BIC,tmp_BIC)
    type = c(type,tmp_type)
  }
  names(BIC) = type
  BIC = sort(BIC,decreasing = F)
  cat("Order of the preference of models\n")
  print(BIC)
  best = min(BIC)
  best_name = names(BIC)[match(min(BIC),BIC)]
  names(best) = best_name
  cat("\nBest Model\n")
  print(best)
  return(names(best))
}

# 6. function to return deltaBIC for different models
#' compare different MADSEQ models, and calculate deltaBIC between models.
#' 
#' compare multiple MadSeq models using BIC. Print out the deltaBIC between current model and the best model
#' @param object1 MadSeqOutput object, which is produced by one of the runModel functions
#' @param object2 MadSeqOutput object
#' @param ... MadSeqOutput object
#' @return a named numeric vector indicating the deltaBIC between current model and the best model 
#' @details the function calculated deltaBIC between the best model and given models. The model with lowest BIC value is the best model. The deltaBIC between models indicate the evidence against higher BIC
#' deltaBIC ~ [0,2]: Not worth more than a bare mention                                                                                                                                                    
#' deltaBIC ~ [2,6]: Positive
#' deltaBIC ~ [6,10]: Strong
#' deltaBIC >10: Very Strong
#' @export
modelComparison = function(object1,object2,...){
  data = list(object1,object2,...)
  n_items = length(data)
  type = NULL
  BIC = NULL
  for (i in 1:n_items){
    tmp = data[[i]]
    tmp_BIC = tmp$BIC
    tmp_type = substr(names(tmp_BIC),5,nchar(names(tmp_BIC)))
    BIC = c(BIC,tmp_BIC)
    type = c(type,tmp_type)
  }
  names(BIC) = type
  BIC = sort(BIC,decreasing = F)
  
  deltaBIC = BIC-min(BIC)
  return(deltaBIC)
}

#-------------------------3. Function to process output from Model-------------------
# 1. plot the density plots for the posterior distribution
#' Plot the posterior density estimates from runModel output
#' 
#' Display density plots of the posterior estimates of choosen variables, the density is produced by density function.
#' @param object A MadSeqOutput object, which is produced by one of the runModel functions
#' @param variable A character vector indicates the posterior density of which variable to be plot. (default is plot all the variables)
#' @export
plotPosterior = function(object,variable="all"){
  data = as.data.frame(object$posterior)
  n_total = dim(data)[2]
  all_name = names(data)
  
  # if want to check all the posterior distribution
  if(paste(variable,collapse = "")=="all"){
    par(mfrow=c(3,2))
    par(mar=c(2.5,4,3,2))
    for (i in all_name){
      tmp_dat = data[,i]
      if (min(tmp_dat)==max(tmp_dat)){
        hist(tmp_dat,prob=TRUE,main=paste("posterior density of",i),xlab="")
      }
      else if (i=="p1"|i=="p2"|i=="p3"|i=="p4"|i=="p5"){
        par(lwd=2)
        hist(tmp_dat,prob=TRUE,main=paste("posterior density of",i),xlab="",border="blue")
      }
      else{
        plot(density(tmp_dat),xlab="",lwd=2,col="blue",main=paste("posterior density of",i))
      }
    }
  }
  
  # only want to check some of the posterior distribution
  else{
    n_plot = length(variable)
    if (n_plot>=6) par(mfrow=c(3,2))
    else if (n_plot>1) par(mfrow=c(ceiling(n_plot/2),2))
    else if (n_plot==1) par(mfrow=c(1,1))
    par(mar=c(2.5,4,3,2))
    for (i in variable){
      tmp_dat = data[,i]
      if (min(tmp_dat)==max(tmp_dat)){
        hist(tmp_dat,prob=TRUE,main=paste("posterior density of",i),xlab="")
      }
      else{
        plot(density(tmp_dat),xlab="",lwd=2,col="blue",main=paste("posterior density of",i))
      }
    }
  }
}


# 2. plot the posterior estimated distribution of the fraction, marked with mean and HPD intervals you choose
#' plot the posterior distribution of Fraction estimate
#' 
#' display the histogram of posterior distribution of the fraction of aneuploidy cells, with mean and HPD (highest posterior density) intervals indicated on the plot
#' 
#' @param object A MadSeqOutput object, which is produced by one of the runModel functions
#' @param prob A numeric number in the interval(0,1) giving the probability content of the intervals
#' @export
#' @seealso \code{\link{plotPosterior}} \code{\link{plotMixture}}
plotFraction = function(object,prob=0.95){
  f = object$posterior[,"f"]
  f_density = density(f)
  y_max = max(hist(f,breaks=40,plot=FALSE)$density)
  f_mean = signif(mean(f),2)
  f_HDI = signif(HDIofMCMC(f,credMass = prob),2)
  hist(f,prob=TRUE,breaks=40,xlab="",ylab="",main="posterior estimate of the fraction",col=rgb(0,1,0,0.2),mgp=c(2,1,0),ylim=c(0,1.1*y_max),xlim=c(min(f),max(f)))
  lines(f_density,col="black",lty=2,lwd=3)
  lines(x=c(f_mean,f_mean),y=c(0,1.01*y_max),col="red",lwd=3)
  text((max(f)+min(f))/2,1.03*y_max,paste("mean =",as.character(f_mean)),cex=1.1)
  lines(x=c(f_HDI[1],f_HDI[1]),y=c(0,0.25*y_max),col="red",lwd=3)
  lines(x=c(f_HDI[2],f_HDI[2]),y=c(0,0.25*y_max),col="red",lwd=3)
  text((max(f)+min(f))/2,1.09*y_max,paste(as.character(100*prob),"% HPD = [",as.character(f_HDI[1]),",",as.character(f_HDI[2]),"]",sep=""),cex=1)
}


# 3. plot the estimated mixture posterior distribution 
#' plot the posterior estimates of mixture disbribution
#' 
#' display the posterior mixture disbribution of alternative allele frequencies, indicating the mean of each mixture
#' @param object A MadSeqOutput object, which is produced by one of the runModel functions
#' @param labels logical; if TRUE, will label each mixture and give the mean value on the plot
#' @export
#' @seealso  \code{\link{plotPosterior}} \code{\link{plotFraction}}
plotMixture = function(object,labels=TRUE){
  data = object$posterior
  data_variables = colnames(data)
  data_median = apply(data,2,median)
  
  # check if the output is from one mixture model
  if(is.element("mu",data_variables)){
    a = data_median["mu"]*data_median["kappa"]
    b = (1-data_median["mu"])*data_median["kappa"]
    p1 = data_median["p1"]
    p2 = data_median["p2"]
    y_max = max(dbeta(seq(0,1,0.01),a,b))
    mu = signif(a/(a+b),2)
    curve(p1*dbeta(x,a,b),col="blue",lwd=3,xlab="alternative allele frequency",ylab="density",main="posterior distribution of mixtures",mgp=c(2,1,0),ylim=c(0,1.1*y_max),bty="l")
    curve(p2*dbeta(x,1,1),col="grey4",lwd=2,lty=2,add=TRUE)
    if(labels==TRUE){
      text(x= mu,1.05*y_max,paste("m =",as.character(mu)))
    }
  }
  
  # check if the output is from four mixture model
  else if(is.element("mu[3]",data_variables)){
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
    p1 = data_median["p1"]
    p2 = data_median["p2"]
    p3 = data_median["p3"]
    p4 = data_median["p4"]
    p5 = data_median["p5"]
    
    m1 = signif(a1/(a1+b1),2)
    m2 = signif(a2/(a2+b2),2)
    m3 = signif(a3/(a3+b3),2)
    m4 = signif(a4/(a4+b4),2)
    
    y_max = max(c(p1*dbeta(seq(0,1,0.01),a1,b1),p2*dbeta(seq(0,1,0.01),a2,b2),p3*dbeta(seq(0,1,0.01),a3,b3),p4*dbeta(seq(0,1,0.01),a4,b4),p5*dbeta(seq(0,1,0.01),a5,b5)))
    print(y_max)
    curve(p1*dbeta(x,a1,b1),col="blue",lwd=3,xlab="alternative allele frequency",ylab="density",main="posterior distribution of mixtures",mgp=c(2,1,0),xlim=c(0,1),bty="l",ylim=c(0,1.25*y_max))
    curve(p2*dbeta(x,a2,b2),col="blue",lwd=3,add=TRUE)
    curve(p3*dbeta(x,a3,b3),col="green4",lwd=3,add=TRUE)
    curve(p4*dbeta(x,a4,b4),col="green4",lwd=3,add=TRUE)
    curve(p5*dbeta(x,a5,b5),col="grey4",lwd=2,add=TRUE,lty=2)
    
    if(labels==TRUE){
      text(m1,1.05*max(p1*dbeta(seq(0,1,0.01),a1,b1)),paste("m1 =",as.character(m1)))
      text(m2,1.05*max(p2*dbeta(seq(0,1,0.01),a2,b2)),paste("m2 =",as.character(m2)))
      text(m3,1.25*max(p3*dbeta(seq(0,1,0.01),a3,b3)),paste("m3 =",as.character(m3)))
      text(m4,1.25*max(p4*dbeta(seq(0,1,0.01),a4,b4)),paste("m4 =",as.character(m4)))
    }
  }
  
  # then plot two mixtures
  else{
    a1 = data_median["mu[1]"]*data_median["kappa"]
    b1 = (1-data_median["mu[1]"])*data_median["kappa"]
    a2 = data_median["mu[2]"]*data_median["kappa"]
    b2 = (1-data_median["mu[2]"])*data_median["kappa"]
    a3 = 1
    b3 = 1
    
    p1 = data_median["p1"]
    p2 = data_median["p2"]
    p3 = data_median["p3"]
    
    m1 = signif(a1/(a1+b1),2)
    m2 = signif(a2/(a2+b2),2)
    
    y_max = max(c(p1*dbeta(seq(0,1,0.01),a1,b1),p2*dbeta(seq(0,1,0.01),a2,b2)))
    
    curve(p1*dbeta(x,a1,b1),col="blue",lwd=3,xlab="alternative allele frequency",ylab="density",main="posterior distribution of mixtures",mgp=c(2,1,0),xlim=c(0,1),bty="l",ylim=c(0,1.1*y_max))
    curve(p2*dbeta(x,a2,b2),col="blue",lwd=3,add=TRUE)
    curve(p3*dbeta(x,a3,b3),col="grey4",lwd=2,add=TRUE,lty=2)
    
    if(labels==TRUE){
      text(m1+0.1,y_max,paste("m1 =",as.character(m1)))
      text(m2-0.1,y_max,paste("m2 =",as.character(m2)))
    }
  }
}

#------------------------MADSEQOUTPUT OBJECT--------------------------
#-----------------------Constructor--------------------
# constructor function for the class "MadSeqOutput"
#' Produce Object of MadSeqOutput class, which is the output from runModel functions
#' 
#' The function MadSeqOutput is used to represent the output from runModel functions. A MadSeqOutput object is a list with two components: 
#' 1. posterior: the posterior estimation from the MCMC run of the model. 2. BIC: the Bayesian Information Criteria for current model.
#' 
#' @param x A MadSeqOutput object
#' @return An object of class MadSeqOutput. 
#' @return The \code{MadSeqOutput} class has its own \code{print} method
#' @export
MadSeqOutput = function(x){
  # check if x is a list
  if(is.list(x)==FALSE) stop("X must be a list")
  # check if x consists of two components
  if(length(x)!=2) stop("X must have two components")
  # check the elements of X
  if(!is.matrix(x$posterior)) stop("posterior must contain a matrix")
  if(!is.numeric(x$BIC)) stop("BIC must be numeric")
  class(x) = "MadSeqOutput"
  return(x)
}

#' 
#' @rdname MadSeqOutput
as.MadSeqOutput = function(x){
  x = MadSeqOutput(x)
  return(x)
}

#' @rdname MadSeqOutput
is.MadSeqOutput = function(x){
  if(class(x)=="MadSeqOutput") return(TRUE)
  else return(FALSE)
}

#-------------------------Method-----------------------
# 1. print MadSeqOutput object. Print out the mean and sd of posterior distribution and the BIC value
#' print method for MadSeqOutput object
#'
#' print an object of MadSeqOutput class. It will print out the mean and 95% Highest density inverval of the posterior distribution estimated by the runModels function.
#' @param x An object of MadSeqOutput class
#' @param ... further arguments passed to or from other methods.
#' @seealso \code{\link{runModelOne}} \code{\link{runModelMonosomy}} \code{\link{runModelTrisomy}} \code{\link{runModelFour}} \code{\link{MadSeqOutput}} 
#' @method print MadSeqOutput
#' @export 
print.MadSeqOutput = function(x,...){
  mat = x$posterior
  mat = mat[,-match(c("p_cov","r_cov"),colnames(mat))]
  BIC = x$BIC
  cat("Posterior Distribution for Parameters\n")
  post_mean = signif(apply(mat,2,mean),3)
  post_HDI = signif(apply(mat,2,HDIofMCMC),3)
  HDI = paste("[",signif(post_HDI[1,],3),", ",signif(post_HDI[2,],3),"]",sep="")
  post = data.frame(post_mean,HDI)
  names(post) = c("mean","95% HDI")
  print(post)
  cat("\n")
  cat("BIC\n")
  print(unname(BIC))  
}

# 2. summary. Summarize the posterior distribution as a data.frame
#' summary statistics for MadSeqOutput object
#'
#' summary.MadSeqOutput produce the summary statistics for the posterior distribution produced by runModels function. 
#' @param object An object of MadSeqOutput class
#' @param ... further arguments passed to or from other methods.
#' @seealso \code{\link{runModelOne}} \code{\link{runModelMonosomy}} \code{\link{runModelTrisomy}} \code{\link{runModelFour}} \code{\link{MadSeqOutput}}
#' @method summary MadSeqOutput
#' @export 
summary.MadSeqOutput = function(object,...){
  mat = as.data.frame(object$posterior)
  print(summary(mat))
}


