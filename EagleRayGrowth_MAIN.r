#////////////////////////////////////////////////////////////////////
#### EagleRayGrowth: illustrate methods on simulated data v0.4.3 ####
#////////////////////////////////////////////////////////////////////

# rm(list=ls())

### // Load libraries and functions ----
require(TMB)
require(tmbstan)

compile("FabensBayesian.cpp") # only need to run once
dyn.load(dynlib("FabensBayesian")) # to run for every new R session

source('EagleRayGrowth_Methods.r') # create all necessary functions


### // Simulate data: single population ----
# For the sake of illustration, we simulate data here according to the
# Fabens (1965) model. The point is that the user should then be able to format
# their data accordingly to use the methods implemented here.

n <- 100 # sample size, nb of capture-recapture pairs
trueLinf <- 167 # loosely based on aquarium female Aetobatus narinari
trueK <- 0.36   # loosely based on aquarium female Aetobatus narinari

sigma <- 5 # standard deviation (sd) of Gaussian error, for data simulation

set.seed(1234) # for replicability
Lcap <- rnorm(n=n,mean=100,sd=5)
# ^ lengths at capture, mean length of 100 cm (somewhat arbitrary)
deltaT <- rgamma(n=n,shape=4,scale=0.25)
# ^ times at liberty in years, mean time = 1 year

dat <- simFabensFixedCov(Linf=trueLinf,K=trueK,sdeps=sigma,
                         Lcap=Lcap,deltaT=deltaT)
str(dat)
# ^ simulate lengths at recapture ($Lrecap) from supplied Lcap and deltaT
#   according to Fabens (1965) with iid Gaussian errors with mean 0 and sd
#   sigma. meanLrecap is the mean of Lrecap without the Gaussian errors (exact
#   von Bertalanffy growth equation).


### // setup for estimation ----
Lmax <- 186.7
# ^ max length in historical data (disc width of "Nerina" at Discovery Cove)

hp.Linf <- c(5.23955339, 0.07836825) # aquarium female Aetobatus narinari
# ^ we specify a lognormal prior for Linf so this vector of hyperparameters (hp)
#   contains the lognormal mean and sd on the log scale (as in dlnrom). These hp
#   were derived from Lmax/0.99 set as the prior median.
hp.K <- c(1e-2, 1)
hp.sigma <- c(1e-5, 50)
# ^ for K and sigma we set a uniform prior, so these two vectors of
#   hyperparameters contain the lower and upper bounds defining the uniform
#   support

hplist <- list(hp.Linf, hp.K, hp.sigma)

par.ini <- c(Lmax/0.99, 0.5)
# ^ initial values for Linf and K: for Linf it is based on historical Lmax,
#   while for K it is simply in the right ballpark

mcmc.control <- list('nchains'=5,
                     'iter'=20000,'warmup'=10000,
                     'adapt_delta'=0.8,
                     'max_treedepth'=20
)
# ^ passed on to Stan's NUTS


### // Frequentist estimation ----
est1 <- fa65(par=par.ini,
             L1=dat$Lcap,
             L2=dat$Lrecap,
             deltaT=dat$deltaT
)
# ^ fa65 is a least squares estimation of the Fabens (1965) model
str(est1)
# ^ $par are the estimates, $se are the standard errors, $covmat is the 2x2
#   covariance matrix subset corresponding to (Linf,K), and the $residuals are
#   important for model validation.
# ^ $se=NA for sigma is deliberate, we never computed it (not a parameter of
#   interest)



### // Bayesian estimation ----
est2 <- Bfa65(par=par.ini,
              L1=dat$Lcap,
              L2=dat$Lrecap,
              deltaT=dat$deltaT,
              priordist.Linf='lognormal',
              priordist.K='uniform',
              priordist.sigma='uniform',
              hyperpar=hplist,
              mcmc.control=mcmc.control
)
# ^ Bfa65 is the Bayesian version of fa65
# ^ warning about "Re-cycling inits to match number of chains" can be safely
#   ignored
str(est2)
# ^ $par are the main point estimates (posterior medians), $se are the posterior
#   standard deviations, $par.TMB and $se.TMB are from the TMB MAP estimation
#   used as initial values for parallel chains in NUTS.


### // Compare estimates ----
cbind(c(trueLinf,trueK), est1$par[1:2], est2$par[1:2])
# ^ Bayesian estimate is much closer to the true values used for simulation
#   here, the prior definitely helped

cbind(est1$par[1:2]-qnorm(0.975)*est1$se[1:2],
      est1$par[1:2]+qnorm(0.975)*est1$se[1:2])
# ^ 95% confidence interval (based on asymptotic normality, not guaranteed to
#   cover well in this setting)

rbind(est2$cred.int$Linf, est2$cred.int$K)
# ^ 95% credible interval (equal tails for simplicity here)


# TODO:
# * simulate two populations that are different in their growth parameters,
#   apply LRT and BF
# * secondary analysis of Fabens model with random intercept on Linf, simulate
#   >=1 recaptures and apply FitFabensRandef
# * DW-weight relation analysis, including bootstrap for colored envelope



# END EagleRayGrowth_Main
