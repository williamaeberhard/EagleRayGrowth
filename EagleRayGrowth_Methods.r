#///////////////////////////////////////////////////////////////////////////////
#### EagleRayGrowth v0.4.3 ####
#///////////////////////////////////////////////////////////////////////////////


#///////////////////////////////////////////////////////////////////////////////
#### simFabensFixedCov: simulate Lrecap according to Fabens (1965) Gaussian ####
#///////////////////////////////////////////////////////////////////////////////

# FixedCov means that Lcap and deltaT are supplied and fixed, only iid Gaussian
# errors (mean 0 and sd supplied by user) are generated and added to a mean
# vector computed from vB equation given Lcap, deltaT and supplied Linf and K
# parameters. Use in LRT_2pop_fa65: Linf,K,sigma are estimates under constrained
# model under H0 from real/suplied data

simFabensFixedCov <- function(Linf,K,sdeps,Lcap,deltaT){
  # use supplied Lcap and deltaT and only sim Lrecap with given Linf,K,sigma
  # following Fabens/vB with additive iid Gaussian meas err
  
  meanLrecap <- Linf-(Linf-Lcap)*exp(-K*deltaT) # from FabensBayesian v0.2.2b
  Lrecap <- rnorm(n=length(Lcap), mean=meanLrecap, sd=sdeps) # add iid meas err
  
  return(list('Lcap'=Lcap,'Lrecap'=Lrecap,'deltaT'=deltaT, # actual data
              'meanLrecap'=meanLrecap # for reference
  ))
}


#///////////////////////////////////////////////////////////////////////////////
#### fa65: Fabens (1965), frequentist (least squares) ####
#///////////////////////////////////////////////////////////////////////////////

# least squares estimation: no assumption on distribution of errors

fa65 <- function(par,L1,L2,deltaT,meth='nlminb',compute.se=T){
  ### objective function and gradient wrt (Linf,K)
  rss <- function(par,L1,L2,deltaT){
    Linf <- par[1]
    K <- par[2]
    dL.fitted <- (Linf-L1)*(1-exp(-K*deltaT))
    dL <- L2-L1
    RSS <- sum((dL-dL.fitted)^2)
    return(RSS) # residual sum of squares
  }
  gr.rss <- function(par,L1,L2,deltaT){
    Linf <- par[1]
    K <- par[2]
    dL.fitted <- (Linf-L1)*(1-exp(-K*deltaT))
    dL <- L2-L1
    resid <- dL-dL.fitted
    expKdt <- exp(-K*deltaT)
    gr <- -2*c(sum(resid*(1-expKdt)),sum(resid*(Linf-L1)*deltaT*expKdt))
    return(gr)
  }
  ### optimization
  if (meth=='nlminb'){
    opt <- nlminb(start=par,objective=rss,gradient=gr.rss,
                  L1=L1,L2=L2,deltaT=deltaT)
    res <- list('par'=opt$par,'value'=opt$obj,'conv'=opt$mess)
  } else {
    opt <- optim(par=par,fn=rss,gr=gr.rss,
                 L1=L1,L2=L2,deltaT=deltaT,method=meth)
    res <- list('par'=opt$par,'value'=opt$val,'conv'=opt$mess)
  }
  ### optional: compute standard errors
  if (compute.se){
    # setup
    n <- length(L1)
    Linf <- opt$par[1] # plug-in estimated values
    K <- opt$par[2] # plug-in estimated values
    # estimate residual variance under Fabens model (implicit Gaussian error)
    y <- L2-L1 # deltaL
    resvec <- y-(Linf-L1)*(1-exp(-K*deltaT)) # residuals
    sigma2 <- sum(resvec^2)/(n-2) # res sum of squares, df=2 for Linf+K
    # compute M=Q
    expKdt <- exp(-K*deltaT)
    M11 <- sum((1-expKdt)^2)
    M12 <- sum(deltaT*expKdt*(Linf-L1)*(1-expKdt))
    M22 <- sum((deltaT*expKdt*(Linf-L1))^2)
    Minv <- cbind(c(M22,-M12),c(-M12,M11))/(M11*M22-M12^2)
    # compute var-cov sandwich matrix (based on expectations)
    varcov <- sigma2*Minv # (2x2) variance-covariance matrix
    res$se <- c(sqrt(diag(varcov)), NA) # 3rd = NA for se.sigma
    # ^ standard errors, don't bother now computing for sigma
    res$par <- c(res$par,sqrt(sigma2)) # append est for sigma
    res$residuals <- resvec
    res$covmat <- varcov # variance-covariance matrix for (Linf,K)
    names(res$par) <- c('Linf','K','sigma')
    names(res$se) <- c('Linf','K','sigma')
  } else {
    res$se <- rep(NA,2) # neglect entry for sigma if compute.se==F
    res$residuals <- rep(NA,n) # neglect if compute.se==F
    res$covmat <- NA           # neglect if compute.se==F
    names(res$par) <- c('Linf','K') # basic output
    names(res$se) <- c('Linf','K') # basic output
  }
  
  ### output
  return(res)
}


#///////////////////////////////////////////////////////////////////////////////
#### Bfa65: Bayesian Fabens (1965) ####
#///////////////////////////////////////////////////////////////////////////////

# Gaussian distribution assumed for error on lengths at recapture

Bfa65 <- function(par,L1,L2,deltaT,
                  priordist.Linf='lognormal',
                  priordist.K='uniform',
                  priordist.sigma='uniform',
                  hyperpar=NULL,
                  meth='nlminb',compute.se=T,
                  onlyTMB=F,output.post.draws=F,output.mcmc.obj=F,
                  mcmc.control=list('nchains'=5,'iter'=20000,'warmup'=10000,
                                    'adapt_delta'=0.8,'max_treedepth'=20)){
  
  ### setup priors
  
  priordist.code <- integer(3) # codes for Linf, K, and sigma
  
  # prior for Linf
  
  if (priordist.Linf=='uniform'){
    priordist.code[1] <- 1L
    if (!is.list(hyperpar)){ # is.null(hyperpar) # NULL is not a list
      stop('hyperpar must be supplied as a list of 3 elements (Linf, K, sigma)',
           ', where for priordist.Linf="uniform" the 1st one is a vector of ',
           'length 2 of lower and upper bounds.')
    }
  } else if (priordist.Linf=='normal' | priordist.Linf=='Gaussian'){
    priordist.code[1] <- 2L
    if (!is.list(hyperpar)){ # is.null(hyperpar) # NULL is not a list
      stop('hyperpar must be supplied as a list of 3 elements (Linf, K, sigma)',
           ', where for priordist.Linf="normal" the 1st one is a vector of ',
           'length 2 of mean and sd.')
    }
  } else if (priordist.Linf=='lognormal'){
    priordist.code[1] <- 3L
    # hyperpar[[1]] <- c(2*log(hyperpar[[1]][1]) +
    #                      - log(hyperpar[[1]][2]^2+hyperpar[[1]][1]^2)/2,
    #                    sqrt(-2*log(hyperpar[[1]][1]) +
    #                           + log(hyperpar[[1]][2]^2+hyperpar[[1]][1]^2)))
    # #  ^ mean and sd on log scale, user supplies mean and sd on exp scale
    if (!is.list(hyperpar)){ # is.null(hyperpar) # NULL is not a list
      stop('hyperpar must be supplied as a list of 3 elements (Linf, K, sigma)',
           ', where for priordist.Linf="lognormal" the 1st one is a vector of ',
           'length 2 of mean and sd on the log scale.')
    }
  } else { # then no prior specified, straight frequentist Fabens, no MCMC
    priordist.code[1] <- 0L
  }
  
  # prior for K
  
  if (priordist.K=='uniform'){
    priordist.code[2] <- 1L
    if (!is.list(hyperpar)){ # is.null(hyperpar) # NULL is not a list
      stop('hyperpar must be supplied as a list of 3 elements (Linf, K, sigma)',
           ', where for priordist.K="uniform" the 2nd one is a vector of ',
           'length 2 of lower and upper bounds.')
    }
  } else if (priordist.K=='normal' | priordist.K=='Gaussian'){
    priordist.code[2] <- 2L
    if (!is.list(hyperpar)){ # is.null(hyperpar) # NULL is not a list
      stop('hyperpar must be supplied as a list of 3 elements (Linf, K, sigma)',
           ', where for priordist.K="normal" the 2nd one is a vector of ',
           'length 2 of mean and sd.')
    }
  } else if (priordist.K=='lognormal'){
    priordist.code[2] <- 3L
    # hyperpar[[2]] <- c(2*log(hyperpar[[2]][1]) +
    #                      - log(hyperpar[[2]][2]^2+hyperpar[[2]][1]^2)/2,
    #                    sqrt(-2*log(hyperpar[[2]][1]) +
    #                           + log(hyperpar[[2]][2]^2+hyperpar[[2]][1]^2)))
    # # ^ mean and sd on log scale, user supplies mean and sd on exp scale
    if (!is.list(hyperpar)){ # is.null(hyperpar) # NULL is not a list
      stop('hyperpar must be supplied as a list of 3 elements (Linf, K, sigma)',
           ', where for priordist.K="lognormal" the 2nd one is a vector of ',
           'length 2 of mean and sd on the log scale.')
    }
  } else { # then no prior specified, straight frequentist Fabens, no MCMC
    priordist.code[2] <- 0L
  }
  
  # prior for sigma
  
  if (priordist.sigma=='uniform'){
    priordist.code[3] <- 1L
    if (!is.list(hyperpar)){ # is.null(hyperpar) # NULL is not a list
      stop('hyperpar must be supplied as a list of 3 elements (Linf, K, sigma)',
           ', where for priordist.sigma="uniform" the 3rd one is a vector of ',
           'length 2 of lower and upper bounds.')
    }
  } else if (priordist.sigma=='normal' | priordist.sigma=='Gaussian'){
    priordist.code[3] <- 2L
    if (!is.list(hyperpar)){ # is.null(hyperpar) # NULL is not a list
      stop('hyperpar must be supplied as a list of 3 elements (Linf, K, sigma)',
           ', where for priordist.sigma="normal" the 3rd one is a vector of ',
           'length 2 of mean and sd.')
    }
  } else if (priordist.sigma=='lognormal'){
    priordist.code[3] <- 3L
    # hyperpar[[3]] <- c(2*log(hyperpar[[3]][1]) +
    #                      - log(hyperpar[[3]][2]^2+hyperpar[[3]][1]^2)/2,
    #                    sqrt(-2*log(hyperpar[[3]][1]) +
    #                           + log(hyperpar[[3]][2]^2+hyperpar[[3]][1]^2)))
    # # ^ mean and sd on log scale, user supplies mean and sd on exp scale
    if (!is.list(hyperpar)){ # is.null(hyperpar) # NULL is not a list
      stop('hyperpar must be supplied as a list of 3 elements (Linf, K, sigma)',
           ', where for priordist.sigma="lognormal" the 3rd one is a vector of ',
           'length 2 of mean and sd on the log scale.')
    }
  } else { # then no prior specified, straight frequentist Fabens, no MCMC
    priordist.code[3] <- 0L
  }
  
  if (any(priordist.code==0L) & !all(priordist.code==0L)){
    warning('At least one prior distribution specified, but not for all three ',
            'Linf, K, and sigma, so reverting to TMB estimation only and no MCMC.')
  }
  
  ### setup
  if (meth!='nlminb'){warning('Only meth="nlminb" supported for now.')}
  # if (!compute.se){warning('se will be computed anyway.')}
  
  # n <- length(L1)
  
  datalist <- list('L1'=L1,'L2'=L2,
                   'deltaT'=deltaT,
                   'hp_Linf'=hyperpar[[1]], # user-supplied
                   'hp_K'=hyperpar[[2]], # user-supplied
                   'hp_sigma'=hyperpar[[3]], # user-supplied
                   'priordist'=priordist.code)
  parlist <- list('logLinf'=log(par[1]),
                  'logK'=log(par[2]),
                  'logsigma'=0)
  
  ### TMB estimation (not necessary but good benchmark)
  obj <- MakeADFun(data=datalist,parameters=parlist,
                   DLL="FabensBayesian",silent=T)
  opt <- nlminb(start=obj$par,obj=obj$fn,gr=obj$gr,
                control=list(eval.max=5000,iter.max=5000))
  theta.tmb <- exp(opt$par[1:3])
  # ^ naive estimates without sdreport(), incl sigma
  
  ### MCMC based on tmbstan::tmbstan, by default uses NUTS
  if (!onlyTMB & all(priordist.code!=0)){
    mcmc.obj <- tmbstan(obj=obj,lower=rep(-Inf,3),upper=rep(Inf,3),
                        silent=T,laplace=F,
                        chains=mcmc.control$nchains,
                        warmup=mcmc.control$warmup,
                        iter=mcmc.control$iter,
                        control=list(adapt_delta=mcmc.control$adapt_delta,
                                     # ^ def=0.8, larger = safer but slower
                                     max_treedepth=mcmc.control$max_treedepth
                                     # ^ def=10, larger helps transitions a lot
                        ),
                        # init='random'
                        init='last.par.best' # start from TMB's MLE above
    )
    
    # traceplot(mcmc.obj, pars=names(obj$par), inc_warmup=TRUE) # check conv
    # pairs(mcmc.obj, pars=names(obj$par)) # post dist and scatterplots
    # # ^ Linf and K typically correlate a lot negatively (though not linearly)
    
    # extract MCMC post draws for derived quantities specified in obj's REPORT
    mcmc.post <- as.matrix(mcmc.obj)
    mcmc.est <- matrix(NA_real_,nrow=nrow(mcmc.post),ncol=3) # Linf, K, sigma
    for (i in 1:nrow(mcmc.post)){
      mcmc.est[i,] <- unlist(obj$report(mcmc.post[i,-ncol(mcmc.post)])[
        c('Linf','K','sigma')]) # need all param for DIC
    }
    # colMeans(mcmc.est[,1:2]) # post means for Linf and K
    
    res <- list('par'=c(median(mcmc.est[,1]),  # post median for Linf 
                        median(mcmc.est[,2]),  # post median foir K
                        median(mcmc.est[,3])), # post median for sigma
                # ^ post median better, dist can be very skewed with small n
                'par.TMB'=theta.tmb
    )
    # ^ MCMC point estimates are posterior median as of v0.3.1
    names(res$par) <- c('Linf','K','sigma')
    names(res$par.TMB) <- c('Linf','K','sigma')
  } else {
    res <- list('par'=c(NA,NA,NA),'par.TMB'=theta.tmb)
    names(res$par.TMB) <- c('Linf','K','sigma')
    # ^ incl sigma by def
  }
  # res$sigma <- median(mcmc.est[,3])
  
  ### optional: compute standard errors
  if (compute.se){
    if (!onlyTMB & all(priordist.code!=0)){
      # res$se <- c(sqrt(var(mcmc.est[,1])),sqrt(var(mcmc.est[,2])))
      # ^ posterior naive se numerically unstable if low n
      res$se <- c(median(abs(mcmc.est[,1]-res$par[1])),
                  median(abs(mcmc.est[,2]-res$par[2])),
                  median(abs(mcmc.est[,3]-res$par[3]))) # MADAM
      names(res$se) <- c('Linf','K','sigma')
      # res$se.sigma <- median(abs(mcmc.est[,3]-res$sigma)) # MADAM
    } else {
      res$se <- c(NA,NA,NA) # incl sigma by def
      # res$se.sigma <- NA
    }
    rep <- sdreport(obj)
    res$se.TMB <- summary(rep)[c('Linf','K','sigma'),2]
    # ^ delta method std errors
    names(res$se.TMB) <- c('Linf','K','sigma')
  } else {
    res$se.TMB <- c(NA,NA,NA)
    if (!onlyTMB){
      res$se <- c(NA,NA,NA)
    }
  }
  
  ### output
  if (output.post.draws){
    if (onlyTMB | any(priordist.code==0)){
      warning('No posterior draws if onlyTMB=TRUE or not all priordist set.')
    } else {
      res$post.draws <- list('Linf'=mcmc.est[,1],
                             'K'=mcmc.est[,2],
                             'sigma'=mcmc.est[,3])
    }
  }
  
  if (!onlyTMB & all(priordist.code!=0)){
    probvec <- c(0.025,0.975) # equal tails 95%
    res$cred.int <- list('Linf'=quantile(mcmc.est[,1],probs=probvec),
                         'K'=quantile(mcmc.est[,2],probs=probvec),
                         'sigma'=quantile(mcmc.est[,3],probs=probvec)
    )
    # ^ equal-tailed 95% credible intervals based on MCMC draws
    
    # DIC, both p_D and p_V versions
    loglikvec <- function(par,L1,L2,deltaT,log=T){
      Linf <- par[1]
      K <- par[2]
      sigma <- par[3]
      meanL2 <- Linf-(Linf-L1)*exp(-K*deltaT)
      return(dnorm(x=L2, mean=meanL2, sd=sigma, log=log))
    }
    loglik <- function(par,L1,L2,deltaT){
      return(sum(loglikvec(par,L1,L2,deltaT)))
    } # common wrapper for sum of log-likelihoods over obs
    
    postdev <- -2*apply(X=mcmc.est,MARGIN=1,FUN=loglik,
                        L1=L1,L2=L2,deltaT=deltaT)
    # ^ post draws of deviance = -2*loglik
    ll.mean <- loglik(par=colMeans(mcmc.est),
                      L1=L1,L2=L2,deltaT=deltaT)
    
    dic1 <- 2*(mean(postdev)+ll.mean)
    # ^ original DIC def with p_D, Spiegelhalter et al. (2002, JRSSB)
    dic2 <- mean((postdev-mean(postdev))^2)-2*ll.mean
    # ^ alt DIC def with p_V = 0.5*post var of deviance, Gelman et al. (2014)
    res$DIC <- c(dic1,dic2)
    names(res$DIC) <- c('pD','pV')
    
    # WAIC with variance-based "bias correction" term
    llmat <- apply(X=mcmc.est, MARGIN=1, FUN=loglikvec,
                   L1=L1, L2=L2, deltaT=deltaT)
    # ^ n rows, iter cols
    pWAIC2 <- mean(apply(X=llmat,MARGIN=1,FUN=var))
    # ^ sum over obs of emp var over MCMC draws, eq. (7.12) Gelman et al. (2014)
    
    likmat <- apply(X=mcmc.est, MARGIN=1, FUN=loglikvec,
                    L1=L1, L2=L2, deltaT=deltaT, log=F) # post density
    lppd <- sum(log(rowMeans(likmat)))
    # ^ log pointwise pred density, eq. (7.5) Gelman et al. (2014)
    
    res$WAIC <- -2*(lppd-pWAIC2)
    
    if (output.mcmc.obj){
      res$mcmc.obj <- mcmc.obj # to run traceplot() and pairs()
    }
  } else {
    res$cred.int <- NA
    res$DIC <- c(NA,NA)
    res$WAIC <- NA
  }
  
  res$AIC <- 2*opt$obj+2*length(opt$par) # output AIC anyway
  res$value.TMB <- opt$obj
  res$conv.TMB <- opt$mess
  
  return(res)
}


#///////////////////////////////////////////////////////////////////////////////
#### LRT_2pop_fa65: LRT comparison two populations, Fabens (1965) ####
#///////////////////////////////////////////////////////////////////////////////

# LRT = likelihood ratio test

# Gaussian distribution assumed for error on lengths at recapture

# require(TMB)
# # compile("FabensTwoPop.cpp")
# dyn.load(dynlib("FabensTwoPop"))

LRT_2pop_fa65 <- function(par=NULL,alpha=0.05,
                          pvalue='both',nboot=1e4,
                          L1.pop1,L2.pop1,deltaT.pop1,
                          L1.pop2,L2.pop2,deltaT.pop2){
  # pvalue arg can be 'chi2', 'boot', or 'both'
  
  ### setup
  if (is.null(par)){
    par.ini <- c(1,0.5) # bad starting values
    warning('No par supplied, using default ini (Linf,K) = c(1,0.5).')
  } else {
    par.ini <- c(par[1],par[2]) # keep only first two for Linf and K in pop1
  }
  
  if (!pvalue%in%c('chi2','boot','both')){
    stop('pvalue must be "chi2", "boot", or "both"')
  }
  
  datalist <- list(
    'L1_1'=L1.pop1,
    'L2_1'=L2.pop1,
    'deltaT_1'=deltaT.pop1,
    'L1_2'=L1.pop2,
    'L2_2'=L2.pop2,
    'deltaT_2'=deltaT.pop2
  )
  parlist <- list(
    'logLinf'=log(par.ini[1]),
    'logK'=log(par.ini[2]),
    'logsigma1'=0, # pop1
    'deltaLinf'=0,
    'deltaK'=0,
    'logsigma2'=0 # pop2
  )
  
  
  ### fit1: unconstrained, full model where all param estimated
  obj1 <- MakeADFun(data=datalist,
                    parameters=parlist,
                    DLL="FabensTwoPop",silent=T)
  
  opt1 <- nlminb(start=obj1$par,obj=obj1$fn,gr=obj1$gr,
                 control=list(eval.max=2000,iter.max=2000))
  # opt1$mess
  
  # summ.rep1 <- summary(sdreport(obj1)) # as of v0.4.3: don't need se
  summ.rep1 <- obj1$rep()
  
  
  ### fit0: constrained under H0: (deltaLinf,deltaK) = (0,0)
  obj0 <- MakeADFun(data=datalist,
                    parameters=parlist,
                    map=list('deltaLinf'=factor(NA), # H0: fixed at 0
                             'deltaK'=factor(NA)),   # H0: fixed at 0
                    DLL="FabensTwoPop",silent=T)
  
  opt0 <- nlminb(start=obj0$par,obj=obj0$fn,gr=obj0$gr,
                 control=list(eval.max=2000,iter.max=2000))
  # opt0$mess
  
  # summ.rep0 <- summary(sdreport(obj0)) # as of v0.4.3: don't need se
  summ.rep0 <- obj0$rep()
  
  par.M0 <- unlist(summ.rep0[c('Linf','K','sigma1','sigma2')])
  
  
  ### compute LRT test stat, chi2 p-value, param boot p-value
  lrtval <- 2*(opt0$obj-opt1$obj) # LRT stat value, compare to chi^2 with df=2
  
  if (pvalue%in%c('chi2','both')){
    chi2quant <- qchisq(p=alpha,df=2,lower.tail=F)
    chi2pval <- pchisq(q=lrtval,df=2,lower.tail=F)
  } else {
    chi2quant <- NA
    chi2pval <- NA
  } # simpler to report same output list regardless of pvalue option
  
  if (pvalue%in%c('boot','both')){
    ### LRT parametric bootstrap, sim obs from Fabens model under H0
    boot.lrtval <- double(nboot)
    
    for (b in 1:nboot){
      dat1boot <- simFabensFixedCov(Linf=par.M0[1],
                                    K=par.M0[2],
                                    sdeps=par.M0[3],
                                    Lcap=L1.pop1,
                                    deltaT=deltaT.pop1
      )
      dat2boot <- simFabensFixedCov(Linf=par.M0[1],
                                    K=par.M0[2],
                                    sdeps=par.M0[4],
                                    Lcap=L1.pop2,
                                    deltaT=deltaT.pop2
      )
      # ^ param boot: sim with est on ori data but constr under H0: Linf and K
      #   common to both pop but sigma allowed different
      
      datalist.boot <- list(
        'L1_1'=dat1boot$Lcap,
        'L2_1'=dat1boot$Lrecap,
        'deltaT_1'=dat1boot$deltaT,
        'L1_2'=dat2boot$Lcap,
        'L2_2'=dat2boot$Lrecap,
        'deltaT_2'=dat2boot$deltaT
      )
      # ^ boot sample 

      # fit1: unconstrained
      obj1.boot <- MakeADFun(data=datalist.boot,
                        parameters=parlist, # same for all boot samples
                        DLL="FabensTwoPop",silent=T)
      opt1.boot <- nlminb(start=obj1.boot$par,obj=obj1.boot$fn,gr=obj1.boot$gr,
                          control=list(eval.max=2000,iter.max=2000))
      # opt1.boot$mess
      
      # fit0: constrained under H0: (deltaLinf,deltaK) = (0,0)
      obj0.boot <- MakeADFun(data=datalist.boot,
                        parameters=parlist,
                        map=list('deltaLinf'=factor(NA), # H0: fixed at 0
                                 'deltaK'=factor(NA)),   # H0: fixed at 0
                        DLL="FabensTwoPop",silent=T)
      
      opt0.boot <- nlminb(start=obj0.boot$par,obj=obj0.boot$fn,gr=obj0.boot$gr,
                          control=list(eval.max=2000,iter.max=2000))
      # opt0.boot$mess
      
      boot.lrtval[b] <- 2*(opt0.boot$obj-opt1.boot$obj) # LRT stat value
    }
    # summary(boot.lrtval) # some neg values can happen if small sample size
    
    boot.lrtval.ok <- boot.lrtval[boot.lrtval>=0] # exclude neg values
    nboot.ok <- sum(boot.lrtval>=0) # <=nboot
    # summary(boot.lrtval.ok)
    # ^ rough clean-up
    
    bootpval <- sum(boot.lrtval.ok > lrtval)/nboot.ok
    # ^ parametric bootstrap p-value
  } else {
    bootpval <- NA
  }
  
  
  ### output
  par.pop1 <- unlist(summ.rep1[c('Linf','K','sigma1')]) # estimates only
  par.pop2 <- c(unlist(summ.rep1[c('Linf','K')]) +
                  + unlist(summ.rep1[c('deltaLinf','deltaK')]),
                unlist(summ.rep1['sigma2'])) # estimates only
  names(par.pop1) <- c('Linf','K','sigma')
  names(par.pop2) <- c('Linf','K','sigma')
  
  res <- list('lrt.teststat'=lrtval,
              'loglik.full'=-opt1$obj,
              'loglik.constr'=-opt0$obj,
              'chi2.critval'=chi2quant,
              'chi2.pvalue'=chi2pval,
              'boot.pvalue'=bootpval,
              'est.M0'=par.M0,
              'est.M1'=c(par.pop1,par.pop2)
  )
  # ^ simple output for now, no se
  
  return(res)
}


#///////////////////////////////////////////////////////////////////////////////
#### BF_2pop_Bfa65: BF comparison two populations, Bayesian Fabens (1965) ####
#///////////////////////////////////////////////////////////////////////////////

# BF = Bayes factor
# M0 = constrained model, M1 = full model

# Gaussian distribution assumed for error on lengths at recapture

# require(TMB)
# compile("growthestimation/FabensTwoPopBayesian_M0.cpp")
# compile("growthestimation/FabensTwoPopBayesian_M1.cpp")
# dyn.load(dynlib("growthestimation/FabensTwoPopBayesian_M0")) 
# dyn.load(dynlib("growthestimation/FabensTwoPopBayesian_M1"))

BF_2pop_Bfa65 <- function(par=NULL,
                          L1.pop1,L2.pop1,deltaT.pop1,
                          L1.pop2,L2.pop2,deltaT.pop2,
                          priordist.M0,priordist.M1.pop1,priordist.M1.pop2,
                          hyperpar.M0,hyperpar.M1.pop1,hyperpar.M1.pop2){
  
  # priordist: 0=no prior, 1=unif, 2=Gaussian, 3=lognormal
  # hyperpar: list of 3 elements (Linf,K,sigma), each a vector of length 2
  
  ### setup
  if (is.null(par)){
    par.ini <- c(1,0.5) # bad starting values
    warning('No par supplied, using default ini (Linf,K) = c(1,0.5).')
  } else {
    par.ini <- c(par[1],par[2]) # keep only first two for Linf and K
  }
  
  ### M0: constrained model, same (Linf,K) but sigma1 != sigma2
  datalist.M0 <- list(
    'L1_1'=L1.pop1,
    'L2_1'=L2.pop1,
    'deltaT_1'=deltaT.pop1,
    'L1_2'=L1.pop2,
    'L2_2'=L2.pop2,
    'deltaT_2'=deltaT.pop2,
    'hp_Linf'=hyperpar.M0[[1]],
    'hp_K'=hyperpar.M0[[2]],
    'hp_sigma1'=hyperpar.M0[[3]], # same prior for both sigma, for simplicity
    'hp_sigma2'=hyperpar.M0[[3]], # same prior for both sigma, for simplicity
    'priordist'=priordist.M0
    # ^ 0=no prior, 1=unif, 2=Gaussian, 3=lognormal
  )
  parlist.M0 <- list(
    'logLinf'=log(par.ini[1]),
    'logK'=log(par.ini[2]),
    'logsigma1'=0, # default value
    'logsigma2'=0  # default value
  )
  
  obj.M0 <- MakeADFun(data=datalist.M0, parameters=parlist.M0,
                      DLL="FabensTwoPopBayesian_M0",silent=T) # no random effect
  opt.M0 <- nlminb(start=obj.M0$par, obj=obj.M0$fn, gr=obj.M0$gr,
                   control=list(eval.max=5000,iter.max=5000))
  # unlist(obj.M0$rep()[c('Linf','K','sigma1','sigma2')])
  # ^ TMB est = posterior modes, used as ini for Laplace approx of int
  
  
  ### M1: full model, different (Linf,K,sigma) between two pop
  datalist.M1 <- list(
    'L1_1'=L1.pop1,
    'L2_1'=L2.pop1,
    'deltaT_1'=deltaT.pop1,
    'L1_2'=L1.pop2,
    'L2_2'=L2.pop2,
    'deltaT_2'=deltaT.pop2,
    'hp_Linf1'=hyperpar.M1.pop1[[1]],
    'hp_K1'=hyperpar.M1.pop1[[2]],
    'hp_sigma1'=hyperpar.M1.pop1[[3]], # allowed to be different here
    'hp_Linf2'=hyperpar.M1.pop2[[1]],
    'hp_K2'=hyperpar.M1.pop2[[2]],
    'hp_sigma2'=hyperpar.M1.pop2[[3]], # allowed to be different here
    'priordist1'=priordist.M1.pop1,
    'priordist2'=priordist.M1.pop2
    # ^ 0=no prior, 1=unif, 2=Gaussian, 3=lognormal
  )
  parlist.M1 <- list(
    'logLinf1'=log(par.ini[1]),
    'logK1'=log(par.ini[2]),
    'logsigma1'=0, # default value
    'logLinf2'=log(par.ini[1]),
    'logK2'=log(par.ini[2]),
    'logsigma2'=0  # default value
  ) # same par.ini for both pop
  
  obj.M1 <- MakeADFun(data=datalist.M1, parameters=parlist.M1,
                      DLL="FabensTwoPopBayesian_M1",silent=T) # no random effect
  opt.M1 <- nlminb(start=obj.M1$par, obj=obj.M1$fn, gr=obj.M1$gr,
                   control=list(eval.max=5000,iter.max=5000))
  # unlist(obj.M1$rep()[c('Linf1','K1','sigma1','Linf2','K2','sigma2')])
  # ^ TMB est = posterior modes, used as ini for Laplace approx of int
  
  
  ### compute BF from Laplace approx of marginal pdf
  obj.M0.rand <- MakeADFun(data=datalist.M0,
                           parameters=as.list(opt.M0$par),
                           random=names(parlist.M0), # all integrated out
                           DLL="FabensTwoPopBayesian_M0",silent=T)
  # exp(-obj.M0.rand$fn(obj.M0.rand$par)) # marginal pdf of data|M0
  
  obj.M1.rand <- MakeADFun(data=datalist.M1,
                           parameters=as.list(opt.M1$par),
                           random=names(parlist.M1), # all integrated out
                           DLL="FabensTwoPopBayesian_M1",silent=T)
  # exp(-obj.M1.rand$fn(obj.M1.rand$par)) # marginal pdf of data|M1
  
  BF <- as.numeric(exp(-obj.M1.rand$fn(obj.M1.rand$par) +
                         + obj.M0.rand$fn(obj.M0.rand$par)))
  # ^ Laplace-approximated BF = P(data|M1)/P(data|M0), where both numerator and
  #   denominator integrals are approximated by Laplace's method (TMB's obj$fn) 
  
  # BF interpretation table from Kass and Raftery (1995, p. 777, JASA), adapted
  # from Jeffreys (1961) Theory of Probability, 3rd Ed., Appendix B:
  # 1 - 3.2  | Not worth more than a bare mention
  # 3.2 - 10 | Substantial
  # 10 - 100 | Strong
  # >100     | Decisive
  
  
  ### output
  par.M0 <- unlist(obj.M0$rep()[c('Linf','K','sigma1','sigma2')]) # est only
  par.M1 <- unlist(obj.M1$rep()[c('Linf1','K1','sigma1',
                                  'Linf2','K2','sigma2')]) # est only
  
  res <- list('BF'=BF,
              'marg.logpdf.M0'=as.numeric(-obj.M0.rand$fn(obj.M0.rand$par)),
              'marg.logpdf.M1'=as.numeric(-obj.M1.rand$fn(obj.M1.rand$par)),
              'estTMB.M0'=par.M0,
              'estTMB.M1'=par.M1)
  # ^ simple output for now, no se
  
  return(res)
}





#///////////////////////////////////////////////////////////////////////////////
#### FitFabensRandef: Fabens (1965), ML with random intercept on Linf ####
#///////////////////////////////////////////////////////////////////////////////

# ML with iid Gaussian errors on Lrecap and Gaussian dist for randef

FitFabensRandef <- function(par.ini,Lmat,deltaTmat,
                            meanLinf,sdLinf,
                            K,sdeps,mindeltaT){
  ### corresponds to FabensRandef.cpp | v0.4.3
  # nrow(Lmat) = nrow(deltaTmat) = n = nb indep indiv
  # ncol(Lmat) = max(n_i) = max nb of meas incl cap and 1st recap (NAs for indiv
  #                         with fewer recap)
  # ncol(deltaTmat) = max(n_i)-1, all deltaT relative to cap
  
  # require(TMB)
  # system.time(compile("FabensRandef.cpp")) # run only once
  # dyn.load(dynlib("FabensRandef")) # run for every new session
  
  ### setup
  nbmeas <- rowSums(!is.na(Lmat))
  # ^ nb of meas per indiv, incl cap and 1st recap (so all >=2)
  
  
  ### fit
  parlist <- list('logmuinf'=log(par.ini[1]),
                  'logsigmainf'=0,
                  'logK'=log(par.ini[2]),
                  'logsigma'=0,
                  'logLinf'=rep(0,dim(Lmat)[1]) # randeff
  )
  datalist <- list('Lmat'=Lmat,
                   'nbmeas'=nbmeas,
                   'deltaTmat'=deltaTmat
  )
  
  obj <- MakeADFun(data=datalist,parameters=parlist,
                   random=c('logLinf'),
                   DLL="FabensRandef",silent=T)
  
  system.time(opt <- nlminb(start=obj$par,obj=obj$fn,gr=obj$gr,
                            control=list(eval.max=5000,iter.max=5000))) # < 1 s
  # opt$mess # ok
  
  fit.rep <- obj$rep()
  
  
  ### output
  res <- list('par'=unlist(fit.rep[c('muinf','K','sigmainf','sigma')]),
              'randef.Linf'=fit.rep$Linf) # predicted random intercepts
  names(res$par) <- c('mean.Linf','K','sd.Linf','sigma')
  
  return(res)
}






# END SERgrowth_Methods
