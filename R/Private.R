# this file contains internal function, which are not intended to be used explicitly by users.

# function to create dataList for analysis from the input data
creatDataList_One = function(data,N1,N2,plot=T){
  name = deparse(substitute(data))
  nSNP = dim(data)[1]
  z = data$Alt_D
  N = data$Alt_D+data$Ref_D
  Nclust = 2
  mixture = rep(NA,nSNP)
  m = mean(z/N)
  N1 = N1[N1!=0]
  N2 = N2[N2!=0]
  nSites = length(N2)
  # N2 is the coverage for current chrom
  N_cov = N2
  # N1 is the coverage for the control
  m0 = median(N1)
  if (plot==T){
    plotSites(data,name)
    mtext(paste("m = ",as.character(round(mean(N_cov),1)),", m0 = ",as.character(round(m0,1)),sep=""),side=1,line=-1)
  }
  dataList = list(nSNP = nSNP,
                  nSites = nSites,
                  N_cov = N_cov,
                  m = m,
                  m0 = m0,
                  z = z,
                  N = N,
                  mixture = mixture,
                  Nclust = Nclust)
  print(c(nSNP,nSites))
  return(dataList)
}

creatDataList_Two = function(data,N1,N2,plot=1){
  name = deparse(substitute(data))
  nSNP = dim(data)[1]
  z = data$Alt_D
  N = data$Alt_D+data$Ref_D
  Nclust = 3
  mixture = rep(NA,nSNP)
  mixture[match(max(z/N),z/N)]=1
  mixture[match(min(z/N),z/N)]=2
  m = mean(z/N)
  N1 = N1[N1!=0]
  N2 = N2[N2!=0]
  nSites = length(N2)
  N_cov = N2
  m0 = median(N1)
  
  if (plot==1){
    plotSites(data,name)
    mtext(paste("m = ",as.character(round(mean(N_cov),1)),", m0 = ",as.character(round(m0,1)),sep=""),side=1,line=-1)
  }
  
  dataList = list(nSNP = nSNP,
                  nSites = nSites,
                  N_cov = N_cov,
                  m0 = m0,
                  z = z,
                  N = N,
                  m = m,
                  mixture = mixture,
                  Nclust = Nclust)
  #print(c(nSNP,nSites))
  return(dataList)
}

creatDataList_Four = function(data,N1,N2,plot=1){
  name = deparse(substitute(data))
  nSNP = dim(data)[1]
  z = data$Alt_D
  N = data$Alt_D+data$Ref_D
  Nclust = 5
  m = mean(z/N)
  BAF = z/N
  mixture = rep(NA,nSNP)
  mixture[match(sort(BAF[BAF<m],decreasing=T)[5],BAF)]=2
  mixture[match(sort(BAF[BAF>m],decreasing=F)[5],BAF)]=1  
  N1 = N1[N1!=0]
  N2 = N2[N2!=0]
  nSites = length(N2)
  N_cov = N2
  m0 = median(N1)
  
  if (plot==1){
    plotSites(data,name)
    mtext(paste("m = ",as.character(round(mean(N_cov),1)),", m0 = ",as.character(round(m0,1)),sep=""),side=1,line=-1)
  }
  dataList = list(nSNP = nSNP,
                  nSites = nSites,
                  N_cov = N_cov,
                  m0 = m0,
                  z = z,
                  N = N,
                  m = m,
                  mixture = mixture,
                  Nclust = Nclust)
  print(c(nSNP,nSites))
  print(c(m0=m0,m=mean(N_cov)))
  return(dataList)
}

# get the posterior median of each theta from the runModels output
theta_median = function(data){
  start=match("theta[1]",colnames(data))
  end = dim(data)[2]
  theta = data[,start:end]
  theta_m = apply(theta,2,median)
  return(theta_m)
}

# plot the AAF along chromosome given the data
plotSites = function(dat,name){
  x = dat$position
  y = dat$Alt_D/(dat$Alt_D+dat$Ref_D)
  m = mean(y)
  plot(x,y,pch=20,col="blue",ylim=c(0,1),xlab=paste(dat[1,1],"position"),ylab="Alternative Allele Frequency",main=name)
  abline(h=m,lwd=3,lty=2)
}

# calculate the BIC value for one mixture model
calculateBIC1 = function(post,theta,datalist){
  N = datalist$N
  z = datalist$z
  mu = post["mu"]
  kappa = post["kappa"]
  a = mu*kappa
  b = (1-mu)*kappa
  p1 = post["p1"]
  p2 = post["p2"]
  
  log.aaf = 0
  for (i in 1:length(theta)){
    p_theta = p1*dbeta(theta[i],a,b)+p2*dbeta(theta[i],1,1)
    log.like.aaf = log(p_theta*dbinom(z[i],N[i],prob=theta[i]))
    log.aaf = log.aaf+log.like.aaf
  }
  # print(unname(log.aaf))
  # the log likelihood for coverage
  N_cov = datalist$N_cov
  r_cov = post["r_cov"]
  p_cov = post["p_cov"]
  log.cov = 0
  for (i in 1:length(N_cov)){
    log.like.cov = log(dnbinom(N_cov[i],size=r_cov,prob=p_cov))
    log.cov = log.cov + log.like.cov
  }
  # print(log.cov)
  log.sum = log.aaf+log.cov
  deviance = -2*log.sum
  
  # calculate BIC
  # datapoint
  datapoint = datalist$nSites + sum(datalist$N)
  BIC = 3*log(datapoint)+deviance
  return(BIC)
}

# calculate the BIC value for two mixture model
calculateBIC2 = function(post,theta,datalist){
  N = datalist$N
  z = datalist$z
  a1 = post["mu[1]"]*post["kappa"]
  b1 = (1-post["mu[1]"])*post["kappa"]
  a2 = post["mu[2]"]*post["kappa"]
  b2 = (1-post["mu[2]"])*post["kappa"]
  p1 = post["p1"]
  p2 = post["p2"]
  p3 = post["p3"]
  
  log.aaf = 0
  for (i in 1:length(theta)){
    p_theta = p1*dbeta(theta[i],a1,b1)+p2*dbeta(theta[i],a2,b2)+p3*dbeta(theta[i],1,1)
    log.like.aaf = log(p_theta*dbinom(z[i],N[i],prob=theta[i]))
    log.aaf = log.aaf+log.like.aaf
  }
  # print(unname(log.aaf))
  # the log likelihood for coverage
  N_cov = datalist$N_cov
  r_cov = post["r_cov"]
  p_cov = post["p_cov"]
  log.cov = 0
  for (i in 1:length(N_cov)){
    log.like.cov = log(dnbinom(N_cov[i],size=r_cov,prob=p_cov))
    log.cov = log.cov + log.like.cov
  }
  # print(log.cov)
  log.sum = log.aaf+log.cov
  deviance = -2*log.sum
  
  # calculate BIC
  # datapoint
  datapoint = datalist$nSites + sum(datalist$N)
  BIC = 5*log(datapoint)+deviance
  return(BIC)
}

# calculate the BIC value for four mixture model
calculateBIC4 = function(post,theta,datalist){
  N = datalist$N
  z = datalist$z
  a1 = post["mu[1]"]*post["kappa"]
  b1 = (1-post["mu[1]"])*post["kappa"]
  a2 = post["mu[2]"]*post["kappa"]
  b2 = (1-post["mu[2]"])*post["kappa"]
  a3 = post["mu[3]"]*post["kappa"]
  b3 = (1-post["mu[3]"])*post["kappa"]
  a4 = post["mu[4]"]*post["kappa"]
  b4 = (1-post["mu[4]"])*post["kappa"]
  p1 = post["p1"]
  p2 = post["p2"]
  p3 = post["p3"]
  p4 = post["p4"]
  p5 = post["p5"]
  
  log.aaf = 0
  for (i in 1:length(theta)){
    p_theta = p1*dbeta(theta[i],a1,b1)+p2*dbeta(theta[i],a2,b2)+p3*dbeta(theta[i],a3,b3)+p4*dbeta(theta[i],a4,b4)+p5*dbeta(theta[i],1,1)
    log.like.aaf = log(p_theta*dbinom(z[i],N[i],prob=theta[i]))
    log.aaf = log.aaf+log.like.aaf
  }
  # print(unname(log.aaf))
  # the log likelihood for coverage
  N_cov = datalist$N_cov
  r_cov = post["r_cov"]
  p_cov = post["p_cov"]
  log.cov = 0
  for (i in 1:length(N_cov)){
    log.like.cov = log(dnbinom(N_cov[i],size=r_cov,prob=p_cov))
    log.cov = log.cov + log.like.cov
  }
  # print(log.cov)
  log.sum = log.aaf+log.cov
  deviance = -2*log.sum
  
  # calculate BIC
  # datapoint
  datapoint = datalist$nSites + sum(datalist$N)
  BIC = 7*log(datapoint)+deviance
  return(BIC)
}

# get the HDI of the posterior distribution from MCMC chain at a given interval
HDIofMCMC = function( sampleVec , credMass=0.95 ) {
  # Computes highest density interval from a sample of representative values,
  #   estimated as shortest credible interval.
  # Arguments:
  #   sampleVec
  #     is a vector of representative values from a probability distribution.
  #   credMass
  #     is a scalar between 0 and 1, indicating the mass within the credible
  #     interval that is to be estimated.
  # Value:
  #   HDIlim is a vector containing the limits of the HDI
  sortedPts = sort( sampleVec )
  ciIdxInc = floor( credMass * length( sortedPts ) )
  nCIs = length( sortedPts ) - ciIdxInc
  ciWidth = rep( 0 , nCIs )
  for ( i in 1:nCIs ) {
    ciWidth[ i ] = sortedPts[ i + ciIdxInc ] - sortedPts[ i ]
  }
  HDImin = sortedPts[ which.min( ciWidth ) ]
  HDImax = sortedPts[ which.min( ciWidth ) + ciIdxInc ]
  HDIlim = c( HDImin , HDImax )
  return( HDIlim )
}
