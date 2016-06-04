## plot the AAF along chromosome given the data
plotSites = function(
    dat){
    x = dat$start
    y = dat$Alt_D/(dat$Alt_D+dat$Ref_D)
    m = mean(y)
    plot(x, y, pch=20, col="blue", ylim=c(0,1), 
        xlab=paste(dat[1,1],"position"), ylab="Alternative Allele Frequency")
    abline(h=m,lwd=3,lty=2)
}

## get the HDI of the posterior distribution from MCMC at a given interval
HDIofMCMC = function(
    sampleVec, 
    credMass=0.95) {
    sortedPts = sort(sampleVec)
    ciIdxInc = floor(credMass*length(sortedPts))
    nCIs = length(sortedPts) - ciIdxInc
    ciWidth = rep(0,nCIs)
    for (i in 1:nCIs ) {
        ciWidth[i] = sortedPts[i+ciIdxInc] - sortedPts[i]
    }
    HDImin = sortedPts[which.min(ciWidth)]
    HDImax = sortedPts[which.min(ciWidth) + ciIdxInc]
    HDIlim = c(HDImin,HDImax)
    return(HDIlim)
}
