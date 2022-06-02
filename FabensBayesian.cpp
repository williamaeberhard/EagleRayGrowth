// Bayesian formulation of Fabens (1965) | v0.4.3
// Gaussian likelihood assumed (original Fabens only based on moments)
// priors: choice for Linf, K, and sigma (Gaussian, uniform, lognormal)
// model length at recap L2, given observed T1, T2 and L1
#include <TMB.hpp>

template <class Type>
Type neglogdunif(Type x, Type a, Type b, int n){
	// neg log density of unif(a,b), a<b
	// n = sample size, so that largecst dominates other negloglik contrib
	Type largecst = Type(100.0)*n+Type(100.0);
	Type halfdensity = Type(0.5)*log(b-a);
	Type res1 = CppAD::CondExpGt(x, a, halfdensity, largecst);
	Type res2 = CppAD::CondExpLt(x, b, halfdensity, largecst);
	return res1+res2; // neg log density if a<x<b
}

template<class Type>
Type objective_function<Type>::operator() () {

	//--------------------------------------------------------------------------
	// Inputs
	//--------------------------------------------------------------------------

	// Data
	DATA_VECTOR(L1); // lengths at cap, dim n
	DATA_VECTOR(L2); // lengths at recap, dim n
	// DATA_VECTOR(T1); // time at cap, dim n
	// DATA_VECTOR(T2); // time at recap, dim n
	DATA_VECTOR(deltaT); // time diff betwee cap and recap, dim n
	// ^ as of v0.3: directly supply deltaT

	// Parameters
	PARAMETER(logLinf); // vB Linf
	PARAMETER(logK); // vB growth rate
	PARAMETER(logsigma); // additive Gaussian error sd

	// Random effects: none

	// Hyperparameters
	DATA_VECTOR(hp_Linf); // hyperparam for prior on Linf, dim 2
	DATA_VECTOR(hp_K); // hyperparam for prior on K, dim 2
	DATA_VECTOR(hp_sigma); // hyperparam for prior on sigma, dim 2
	// ^ hp have different meaning depending on priordist:
	//    - priordist=0: unused
	//    - priordist=1: lower and upper bound of uniform density
	//    - priordist=2: mean and sd of Gaussian density
	//    - priordist=3: mean and sd of Gaussian density on log scale (lognormal)

	// Misc
	// DATA_INTEGER(priordist);
	DATA_IVECTOR(priordist); // code for (Linf,K,sigma), dim 3
	// ^ 0 = no prior, 1 = uniform, 2 = Gaussian, 3 = lognormal


	//--------------------------------------------------------------------------
	// Setup, procedures and init
	//--------------------------------------------------------------------------

	int n = L1.size(); // i = 1, ..., n

	Type Linf = exp(logLinf);
	Type K = exp(logK);
	Type sigma = exp(logsigma);
	
	Type nll = 0.0; // init neg loglik


	//--------------------------------------------------------------------------
	// Priors
	//--------------------------------------------------------------------------

	if (priordist(0)==1){ // uniform prior on Linf
		nll += neglogdunif(Linf, hp_Linf(0), hp_Linf(1), n);
	} else if (priordist(0)==2){ // Gaussian prior on Linf
		nll -= dnorm(Linf, hp_Linf(0), hp_Linf(1), true);
	} else if (priordist(0)==3){ // lognormal prior on Linf
		nll -= dnorm(logLinf, hp_Linf(0), hp_Linf(1), true) - logLinf;
		// ^ lognormal log-pdf evaluated at param on exp scale
	} // else no prior on Linf

	if (priordist(1)==1){ // uniform prior on K
		nll += neglogdunif(K, hp_K(0), hp_K(1), n);
	} else if (priordist(1)==2){ // Gaussian prior on K
		nll -= dnorm(K, hp_K(0), hp_K(1), true);
	} else if (priordist(1)==3){ // lognormal prior on K
		nll -= dnorm(logK, hp_K(0), hp_K(1), true) - logK;
		// ^ lognormal log-pdf evaluated at param on exp scale
	} // else no prior on K

	if (priordist(2)==1){ // uniform prior on sigma
		nll += neglogdunif(sigma, hp_sigma(0), hp_sigma(1), n);
	} else if (priordist(2)==2){ // Gaussian prior on sigma
		nll -= dnorm(sigma, hp_sigma(0), hp_sigma(1), true);
	} else if (priordist(2)==3){ // lognormal prior on sigma
		nll -= dnorm(logsigma, hp_sigma(0), hp_sigma(1), true) - logsigma;
		// ^ lognormal log-pdf evaluated at param on exp scale
	} // else no prior on sigma


	//--------------------------------------------------------------------------
	// Random effects
	//--------------------------------------------------------------------------

	// none
	

	//--------------------------------------------------------------------------
	// Observation equations
	//--------------------------------------------------------------------------

	for (int i = 0; i < n; i++){
		// Type deltaT = T2(i)-T1(i); // time increment between cap and recap
		// Type meanL2 = Linf-(Linf-L1(i))*exp(-K*deltaT); // vB at recap
		Type meanL2 = Linf-(Linf-L1(i))*exp(-K*deltaT(i)); // vB at recap
		nll -= dnorm(L2(i), meanL2, sigma, true); // Gaussian lkhd
		// ^ vB reference point = cap
	}

	//--------------------------------------------------------------------------
	// Outputs
	//--------------------------------------------------------------------------

	REPORT(Linf);
	REPORT(K);
	REPORT(sigma);

	ADREPORT(Linf);
	ADREPORT(K);
	ADREPORT(sigma);

	return nll;
}
