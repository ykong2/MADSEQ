######################## create datalist for MCMC ##################
# function to create dataList for analysis from the input data
creatDataList_One = function(
    data,
    control_coverage,
    data_coverage){
    nSNP = dim(data)[1]
    z = data$Alt_D
    N = data$Alt_D+data$Ref_D
    Nclust = 2
    mixture = rep(NA,nSNP)
    m = mean(z/N)
    control_coverage = control_coverage[control_coverage!=0]
    data_coverage = data_coverage[data_coverage!=0]
    nSites = length(data_coverage)
    ## coverage for current chrom
    N_cov = data_coverage
    ## coverage for the control
    m0 = median(control_coverage)
    dataList = list(nSNP = nSNP,
                    nSites = nSites,
                    N_cov = N_cov,
                    m = m,
                    m0 = m0,
                    z = z,
                    N = N,
                    mixture = mixture,
                    Nclust = Nclust)
    return(dataList)
}


creatDataList_Two = function(
    data,
    control_coverage,
    data_coverage){
    nSNP = dim(data)[1]
    z = data$Alt_D
    N = data$Alt_D+data$Ref_D
    Nclust = 3
    mixture = rep(NA,nSNP)
    mixture[match(max(z/N),z/N)]=1
    mixture[match(min(z/N),z/N)]=2
    m = mean(z/N)
    control_coverage = control_coverage[control_coverage!=0]
    data_coverage = data_coverage[data_coverage!=0]
    nSites = length(data_coverage)
    ## coverage for current chrom
    N_cov = data_coverage
    ## coverage for the control
    m0 = median(control_coverage)
    dataList = list(nSNP = nSNP,
                    nSites = nSites,
                    N_cov = N_cov,
                    m0 = m0,
                    z = z,
                    N = N,
                    m = m,
                    mixture = mixture,
                    Nclust = Nclust)
    return(dataList)
}


creatDataList_LOH = function(
    data,
    control_coverage,
    data_coverage){
    nSNP = dim(data)[1]
    z = data$Alt_D
    N = data$Alt_D+data$Ref_D
    mixture = rep(NA,nSNP)
    m = mean(z/N)
    control_coverage = control_coverage[control_coverage!=0]
    data_coverage = data_coverage[data_coverage!=0]
    nSites = length(data_coverage)
    ## coverage for current chrom
    N_cov = data_coverage
    ## coverage for the control
    m0 = median(control_coverage)
    
    dataList = list(nSNP = nSNP,
                    nSites = nSites,
                    N_cov = N_cov,
                    m0 = m0,
                    z = z,
                    N = N,
                    m = m,
                    mixture = mixture)
    return(dataList)
}

creatDataList_Four = function(
    data,
    control_coverage,
    data_coverage){
    nSNP = dim(data)[1]
    z = data$Alt_D
    N = data$Alt_D+data$Ref_D
    Nclust = 5
    m = mean(z/N)
    BAF = z/N
    mixture = rep(NA,nSNP)
    mixture[match(sort(BAF[BAF<m],decreasing=TRUE)[5],BAF)]=2
    mixture[match(sort(BAF[BAF>m],decreasing=FALSE)[5],BAF)]=1  
    control_coverage = control_coverage[control_coverage!=0]
    data_coverage = data_coverage[data_coverage!=0]
    nSites = length(data_coverage)
    ## coverage for current chrom
    N_cov = data_coverage
    ## coverage for the control
    m0 = median(control_coverage)
    dataList = list(nSNP = nSNP,
                    nSites = nSites,
                    N_cov = N_cov,
                    m0 = m0,
                    z = z,
                    N = N,
                    m = m,
                    mixture = mixture,
                    Nclust = Nclust)
    return(dataList)
}
##################### invoke JAGS to run MCMC #################
## running normal model to fit the data
runNormal = function(
    data,
    data_coverage,
    control_coverage,
    checkConvergence,
    adapt,
    burnin,
    nChain,
    nStep,
    thinSteps){
    ## parameters for normal model
    parameters_one = c("kappa","p_cov","r_cov","mu","m_cov")
    ## generate dataList for normal model MCMC
    dataList = creatDataList_One(data,control_coverage,data_coverage)
    
    ## run normal model
    message("1. running normal model")
    fpath = system.file("rjags","OneMix.txt",package="MADSEQ")
    jagsModel_one = jags.model(fpath, data=dataList,
                                n.chains=nChain, n.adapt=adapt)
    update(jagsModel_one, n.iter=burnin)
    codaSamples_one = coda.samples(jagsModel_one, 
                                    variable.names=parameters_one,
                                    n.iter=ceiling((nStep*thinSteps)/nChain), 
                                    thin=thinSteps)
    
    ## check convergence of posterior distribution
    if(checkConvergence==TRUE){
        diag = gelman.diag(codaSamples_one,multivariate = FALSE)
        upper = max(diag$psrf[,2],na.rm=TRUE)
        if (upper<=1.5) message("normal model is converged!")
        else if (upper > 1.5 & upper <= 10){
            message("normal model is not fully converged")
            cat(diag)
        } 
        else if (upper>10){
            message("normal model is not converged")
            cat(diag)
        }
    }
    
    ## convert posterior distribution to a matrix
    one = as.matrix(codaSamples_one)
    ## calculate the BIC value for this model
    BIC = calculate_BIC(one,dataList,type="one")
    names(BIC) = "BIC_normal"
    res = list(one,BIC)
    res
}

## running monosomy model to fit the data
runMonosomy = function(
    data,
    data_coverage,
    control_coverage,
    checkConvergence,
    adapt,
    burnin,
    nChain,
    nStep,
    thinSteps){
    ## parameters for monosomy model
    parameters_two = c("p_cov","m_cov","r_cov","mu","kappa","f")
    # generate dataList for monosomy model MCMC
    dataList = creatDataList_Two(data,control_coverage,data_coverage)
    
    ## run monosomy model
    message("2. running monosomy model")
    fpath = system.file("rjags","TwoMixMonosomy.txt",package="MADSEQ")
    jagsModel_Monosomy = jags.model(fpath, data=dataList,
                                    n.chains=nChain, n.adapt=adapt)
    update(jagsModel_Monosomy, n.iter=burnin)
    codaSamples_Monosomy=coda.samples(jagsModel_Monosomy, 
                                    variable.names=parameters_two, 
                                    n.iter=ceiling((nStep*thinSteps)/nChain),
                                    thin=thinSteps)
    ## check convergence of posterior distribution
    if(checkConvergence == TRUE){
        diag = gelman.diag( codaSamples_Monosomy ,multivariate = FALSE)
        upper = max(diag$psrf[,2],na.rm=TRUE)
        if (upper<=1.5) message("monosomy model is converged!")
        else if (upper > 1.5 & upper <= 10){
            message("monosomy model is not fully converged")
            cat(diag)
        } 
        else if (upper>10){
            message("monosomy model is not converged")
            cat(diag)
        }
    }
    
    ## convert posterior distribution to a matrix
    monosomy = as.matrix(codaSamples_Monosomy)
    ## calculate the BIC value for this model
    BIC = calculate_BIC(monosomy,dataList,type="two")
    names(BIC) = "BIC_monosomy"
    res = list(monosomy,BIC)
    res
}

## running mitotic trisomy model to fit the data
runMitoticTrisomy = function(
    data,
    data_coverage,
    control_coverage,
    checkConvergence,
    adapt,
    burnin,
    nChain,
    nStep,
    thinSteps){
    ## parameters for mitotic trisomy model
    parameters_two = c("p_cov","m_cov","r_cov","mu","kappa","f")
    # generate dataList for monosomy model MCMC
    dataList = creatDataList_Two(data,control_coverage,data_coverage)
    
    # run mitotic trisomy model
    message("3. running mitotic trisomy model")
    fpath = system.file("rjags","TwoMixTrisomy.txt",package="MADSEQ")
    jagsModel_Trisomy = jags.model(fpath, data=dataList,
                                    n.chains=nChain, n.adapt=adapt)
    update(jagsModel_Trisomy, n.iter=burnin )
    codaSamples_Trisomy=coda.samples(jagsModel_Trisomy, 
                                    variable.names=parameters_two, 
                                    n.iter=ceiling((nStep*thinSteps)/nChain),
                                    thin=thinSteps)
    
    # check convergence of posterior distribution
    if(checkConvergence==TRUE){
        diag = gelman.diag( codaSamples_Trisomy ,multivariate = FALSE)
        upper = max(diag$psrf[,2],na.rm=TRUE)
        if (upper<=1.5) message("mitotic trisomy model is converged!")
        else if (upper > 1.5 & upper <= 10){
            message("mitotic trisomy model is not fully converged")
            cat(diag)
        } 
        else if (upper>10){
            message("mitotic trisomy model is not converged")
            cat(diag)
        }
    }
    
    # convert posterior distribution to a matrix
    mitotic_trisomy = as.matrix(codaSamples_Trisomy)
    
    # calculate the BIC value for this model
    BIC = calculate_BIC(mitotic_trisomy,dataList,type="two")
    names(BIC) = "BIC_mitotic_trisomy"
    res = list(mitotic_trisomy,BIC)
    res
}

## running meiotic trisomy model to fit the data
runMeioticTrisomy = function(
    data,
    data_coverage,
    control_coverage,
    checkConvergence,
    adapt,
    burnin,
    nChain,
    nStep,
    thinSteps){
    ## parameters for meiotic trisomy model
    parameters_four = c("p_cov","m_cov","r_cov","mu","kappa","p","f")
    ## generate dataList for the MCMC
    dataList = creatDataList_Four(data,control_coverage,data_coverage)
    ## run meiotic trisomy model
    message("4. running meiotic trisomy model")
    fpath = system.file("rjags","FourMix.txt",package="MADSEQ")
    jagsModel_four = jags.model(fpath, data=dataList,
                                n.chains=nChain, n.adapt=adapt)
    update(jagsModel_four, n.iter=burnin)
    codaSamples_four = coda.samples(jagsModel_four, 
                                    variable.names=parameters_four, 
                                    n.iter=ceiling((nStep*thinSteps)/nChain),
                                    thin=thinSteps)
    
    # check convergence of posterior distribution
    if(checkConvergence==TRUE){
        diag = gelman.diag( codaSamples_four ,multivariate = FALSE)
        upper = max(diag$psrf[,2],na.rm=TRUE)
        if (upper<=1.5) message("meiotic trisomy model is converged!")
        else if (upper > 1.5 & upper <= 10){
            message("meiotic trisomy model is not fully converged")
            cat(diag)
        } 
        else if (upper>10){
            message("meiotic trisomy model is not converged")
            cat(diag)
        }
    }
    
    # convert posterior distribution to a matrix
    meiotic_trisomy = as.matrix(codaSamples_four)
    # calculate the BIC value for this model
    BIC = calculate_BIC(meiotic_trisomy,dataList,type="four")
    names(BIC) = "BIC_meiotic_trisomy"
    res = list(meiotic_trisomy,BIC)
    res
}

## running LOH model to fit the data
runLOH = function(
    data,
    data_coverage,
    control_coverage,
    checkConvergence,
    adapt,
    burnin,
    nChain,
    nStep,
    thinSteps){
    ## parameters for LOH model
    parameters_LOH = c("p_cov","m_cov","r_cov","kappa","f","cgp1","cgp2","d1","d2","m")
    ## generate dataList for the MCMC
    dataList = creatDataList_LOH(data,control_coverage,data_coverage)
    ## run LOH model
    message("5. running loss of heterozygosity model")
    fpath = system.file("rjags","LOH.txt",package="MADSEQ")
    jagsModel_LOH = jags.model(fpath, data=dataList,
                                n.chains=nChain, n.adapt=adapt)
    update(jagsModel_LOH, n.iter=burnin)
    codaSamples_LOH = coda.samples(jagsModel_LOH, 
                                    variable.names=parameters_LOH, 
                                    n.iter=ceiling((nStep*thinSteps)/nChain),
                                    thin=thinSteps)
    
    # check convergence of posterior distribution
    if(checkConvergence==TRUE){
        diag = gelman.diag( codaSamples_LOH ,multivariate = FALSE)
        upper = max(diag$psrf[,2],na.rm=TRUE)
        if (upper<=1.5) message("LOH model is converged!")
        else if (upper > 1.5 & upper <= 10){
            message("LOH model is not fully converged")
            cat(diag)
        } 
        else if (upper>10){
            message("LOH model is not converge")
            cat(diag)
        }
    }
    
    # convert posterior distribution to a matrix
    LOH = as.matrix(codaSamples_LOH)
    
    # calculate the BIC value for this model
    BIC = calculate_BIC(LOH,dataList,type="UPD")
    names(BIC) = "BIC_LOH"
    res = list(LOH,BIC)
    res
}

################# Likelihood and BIC ################
# calculate log maximum likelihood for normal model from posterior distribution
log_likelihood_one = function(
    posterior,
    dataList){
    parameter = apply(posterior,2,median)
    z = dataList$z
    N = dataList$N
    a = parameter["mu"]*parameter["kappa"]
    b = (1-parameter["mu"])*parameter["kappa"]
    p1 = 0.99
    p2 = 0.01
    
    # likelihood for aaf
    loglike_aaf = 0
    for (i in 1:length(z)){
        sub_loglike =
            log(p1*dbetabinom.ab(z[i],N[i],a,b)+
                    p2*dbetabinom.ab(z[i],N[i],1,1))
        loglike_aaf = loglike_aaf + sub_loglike
    }
    # likelihood for coverage
    N_cov = dataList$N_cov
    r_cov = parameter["r_cov"]
    p_cov = parameter["p_cov"]
    log_cov = 0
    for (i in 1:length(N_cov)){
        sub_log_cov = stats::dnbinom(N_cov[i],size=r_cov,prob=p_cov,log=TRUE)
        log_cov = log_cov + sub_log_cov
    }
    #print(log_cov)
    #print(log.cov)
    log_sum = loglike_aaf +log_cov
    return(log_sum)  
}

# calculate log maximum likelihood for twomix model from posterior distribution
log_likelihood_two = function(
    posterior,
    dataList){
    parameter = apply(posterior,2,median)
    z = dataList$z
    N = dataList$N
    a1 = parameter["mu[1]"]*parameter["kappa"]
    b1 = (1-parameter["mu[1]"])*parameter["kappa"]
    a2 = parameter["mu[2]"]*parameter["kappa"]
    b2 = (1-parameter["mu[2]"])*parameter["kappa"]
    p1 = 0.495
    p2 = 0.495
    p3 = 0.01
    # likelihood for aaf
    loglike_aaf = 0
    for (i in 1:length(z)){
        sub_loglike = log(p1*dbetabinom.ab(z[i],N[i],a1,b1)+
                            p2*dbetabinom.ab(z[i],N[i],a2,b2)+
                            p3*dbetabinom.ab(z[i],N[i],1,1))
        loglike_aaf = loglike_aaf + sub_loglike
    }
    # likelihood for coverage
    N_cov = dataList$N_cov
    r_cov = parameter["r_cov"]
    p_cov = parameter["p_cov"]
    log_cov = 0
    for (i in 1:length(N_cov)){
        sub_log_cov = stats::dnbinom(N_cov[i],size=r_cov,prob=p_cov,log=TRUE)
        log_cov = log_cov + sub_log_cov
    }
    log_sum = loglike_aaf +log_cov
    return(log_sum)  
}

# calculate log maximum likelihood for FOUR model from posterior distribution
log_likelihood_four = function(
    posterior,
    dataList){
    parameter = apply(posterior,2,median)
    z = dataList$z
    N = dataList$N
    a1 = parameter["mu[1]"]*parameter["kappa"]
    b1 = (1-parameter["mu[1]"])*parameter["kappa"]
    a2 = parameter["mu[2]"]*parameter["kappa"]
    b2 = (1-parameter["mu[2]"])*parameter["kappa"]
    a3 = parameter["mu[3]"]*parameter["kappa"]
    b3 = (1-parameter["mu[3]"])*parameter["kappa"]
    a4 = parameter["mu[4]"]*parameter["kappa"]
    b4 = (1-parameter["mu[4]"])*parameter["kappa"]
    p1 = parameter["p[1]"]
    p2 = parameter["p[2]"]
    p3 = parameter["p[3]"]
    p4 = parameter["p[4]"]
    p5 = 0.01
    
    # likelihood for aaf
    loglike_aaf = 0
    for (i in 1:length(z)){
        sub_loglike =
            log(p1*dbetabinom.ab(z[i],N[i],a1,b1)+
                    p2*dbetabinom.ab(z[i],N[i],a2,b2)+
                    p3*dbetabinom.ab(z[i],N[i],a3,b3)+
                    p4*dbetabinom.ab(z[i],N[i],a4,b4)+
                    p5*dbetabinom.ab(z[i],N[i],1,1))
        loglike_aaf = loglike_aaf + sub_loglike
    }
    #print(loglike_aaf)
    
    # likelihood for coverage
    N_cov = dataList$N_cov
    r_cov = parameter["r_cov"]
    p_cov = parameter["p_cov"]
    log_cov = 0
    for (i in 1:length(N_cov)){
        sub_log_cov = stats::dnbinom(N_cov[i],size=r_cov,prob=p_cov,log=TRUE)
        log_cov = log_cov + sub_log_cov
    }
    
    log_sum = loglike_aaf +log_cov
    return(log_sum)  
}

# calculate log maximum likelihood for UPD model using posterior distribution
log_likelihood_UPD = function(
    posterior,
    dataList){
    parameter = apply(posterior,2,median)
    z = dataList$z
    N = dataList$N
    a1 = (parameter["m"]+parameter["d1"])*parameter["kappa"]
    b1 = (1-parameter["m"]-parameter["d1"])*parameter["kappa"]
    a2 = (parameter["m"]-parameter["d2"])*parameter["kappa"]
    b2 = (1-parameter["m"]+parameter["d2"])*parameter["kappa"]
    a3 = parameter["m"]*parameter["kappa"]
    b3 = (1-parameter["m"])*parameter["kappa"]
    cgp1 = round(parameter["cgp1"])
    cgp2 = round(parameter["cgp2"])
    
    # likelihood for aaf
    loglike_aaf = 0
    for (i in 1:length(z)){
        if(i<cgp1){
            sub_loglike = log(0.99*dbetabinom.ab(z[i],N[i],a3,b3)+
                                0.01*dbetabinom.ab(z[i],N[i],1,1))
        }
        else if (i>=cgp1 & i<=cgp2){
            sub_loglike =log(0.495*dbetabinom.ab(z[i],N[i],a1,b1)+
                            0.495*dbetabinom.ab(z[i],N[i],a2,b2)+
                            0.01*dbetabinom.ab(z[i],N[i],1,1))
        }
        else if (i >cgp2){
            sub_loglike = log(0.99*dbetabinom.ab(z[i],N[i],a3,b3)+
                                0.01*dbetabinom.ab(z[i],N[i],1,1))
        }
        loglike_aaf = loglike_aaf + sub_loglike
    }
    
    # likelihood for coverage
    N_cov = dataList$N_cov
    r_cov = parameter["r_cov"]
    p_cov = parameter["p_cov"]
    log_cov = 0
    for (i in 1:length(N_cov)){
        sub_log_cov = dnbinom(N_cov[i],size=r_cov,prob=p_cov,log=TRUE)
        log_cov = log_cov + sub_log_cov
    }
    log_sum = loglike_aaf +log_cov
    return(log_sum)
}

# calculate BIC value for different models
calculate_BIC = function(
    posterior,
    dataList,
    type="two"){
    datapoint = dataList$nSites+sum(dataList$N)
    if (type=="one"){
        loglike = log_likelihood_one(posterior,dataList)
        BIC = -2*loglike + 1*log(datapoint)
    }
    else if (type=="two"){
        loglike = log_likelihood_two(posterior,dataList)
        BIC = -2*loglike + 2*log(datapoint)
    }
    else if (type=="four"){
        loglike = log_likelihood_four(posterior,dataList)
        BIC = -2*loglike + 3*log(datapoint)
    }
    else if (type=="UPD"){
        loglike = log_likelihood_UPD(posterior,dataList)
        BIC = -2*loglike + 4*log(datapoint)
    }
    return(BIC)
}