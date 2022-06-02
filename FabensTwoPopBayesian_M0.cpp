// Fabens (1965) estimation for two populations, Bayesian | v0.4.3
// Gaussian likelihood assumed (original Fabens only based on moments)
// * model length at recap L2, given observed T1, T2 and L1
// * priors specified by priordist1 and priordist2 arguments for (Linf,K,sigma)
//   for both pop 1, respectively.
// * M0 model: common (Linf,K) for both pop but distinct sigma, 4 param overall

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
	DATA_VECTOR(L1_1); // lengths at cap, pop 1, dim n1
	DATA_VECTOR(L2_1); // lengths at recap, pop 1, dim n1
	DATA_VECTOR(deltaT_1); // time at liberty, pop 1, dim n1

	DATA_VECTOR(L1_2); // lengths at cap, pop 2, dim n2
	DATA_VECTOR(L2_2); // lengths at recap, pop 2, dim n2
	DATA_VECTOR(deltaT_2); // time at liberty, pop 2, dim n2

	// Parameters
	PARAMETER(logLinf); // vB Linf both pop
	PARAMETER(logK); // vB growth rate both pop
	PARAMETER(logsigma1); // additive Gaussian error sd pop 1
	PARAMETER(logsigma2); // additive Gaussian error sd pop2

	// Random effects: none

	// Hyperparameters
	DATA_VECTOR(hp_Linf); // hyperparam for prior on Linf both pop, dim 2
	DATA_VECTOR(hp_K); // hyperparam for prior on K both pop, dim 2
	DATA_VECTOR(hp_sigma1); // hyperparam for prior on sigma pop 1, dim 2
	DATA_VECTOR(hp_sigma2); // hyperparam for prior on sigma pop 2, dim 2
	// ^ hp have different meaning depending on priordist:
	//   - priordist=0: unused
	//   - priordist=1: lower and upper bound of uniform density
	//   - priordist=2: mean and sd of Gaussian density
	//   - priordist=3: mean and sd of Gaussian density on log scale (lognormal)

	// Misc
	DATA_IVECTOR(priordist); // code for (Linf,K,sigma), dim 3
	// ^ 0 = no prior, 1 = uniform, 2 = Gaussian, 3 = lognormal
	// ^ for simplicity, same prior family enforced for sigma1 and sigma2, with
	//   distinct hp though


	//--------------------------------------------------------------------------
	// Setup, procedures and init
	//--------------------------------------------------------------------------

	int n1 = L1_1.size(); // i = 1, ..., n1
	int n2 = L1_2.size(); // i = 1, ..., n2
	int n = n1+n2;

	Type Linf = exp(logLinf);   // both pop
	Type K = exp(logK);         // both pop
	Type sigma1 = exp(logsigma1); // pop 1
	Type sigma2 = exp(logsigma2); // pop 2
	
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
		// ^ lognormal log-pdf evaluated at param on log scale
	} // else no prior on Linf

	if (priordist(1)==1){ // uniform prior on K
		nll += neglogdunif(K, hp_K(0), hp_K(1), n);
	} else if (priordist(1)==2){ // Gaussian prior on K
		nll -= dnorm(K, hp_K(0), hp_K(1), true);
	} else if (priordist(1)==3){ // lognormal prior on K
		nll -= dnorm(logK, hp_K(0), hp_K(1), true) - logK;
		// ^ lognormal log-pdf evaluated at param on log scale
	} // else no prior on K

	if (priordist(2)==1){ // uniform prior on sigma
		nll += neglogdunif(sigma1, hp_sigma1(0), hp_sigma1(1), n1);
		nll += neglogdunif(sigma2, hp_sigma2(0), hp_sigma2(1), n2);
	} else if (priordist(2)==2){ // Gaussian prior on sigma
		nll -= dnorm(sigma1, hp_sigma1(0), hp_sigma1(1), true);
		nll -= dnorm(sigma2, hp_sigma2(0), hp_sigma2(1), true);
	} else if (priordist(2)==3){ // lognormal prior on sigma
		nll -= dnorm(logsigma1, hp_sigma1(0), hp_sigma1(1), true) - logsigma1;
		nll -= dnorm(logsigma2, hp_sigma2(0), hp_sigma2(1), true) - logsigma2;
		// ^ lognormal log-pdf evaluated at param on log scale
	} // else no prior on sigma


	//--------------------------------------------------------------------------
	// Random effects
	//--------------------------------------------------------------------------

	// none
	

	//--------------------------------------------------------------------------
	// Observation equations
	//--------------------------------------------------------------------------

	for (int i = 0; i < n1; i++){ // pop 1
		Type meanL2 = Linf - (Linf-L1_1(i))*exp(-K*deltaT_1(i)); // vB
		nll -= dnorm(L2_1(i), meanL2, sigma1, true); // Gaussian lkhd, pop 2
	}

	for (int i = 0; i < n2; i++){ // pop 2
		Type meanL2 = Linf - (Linf-L1_2(i))*exp(-K*deltaT_2(i)); // vB
		nll -= dnorm(L2_2(i), meanL2, sigma2, true); // Gaussian lkhd, pop 1
	}


	//--------------------------------------------------------------------------
	// Outputs
	//--------------------------------------------------------------------------

	REPORT(Linf);
	REPORT(K);
	REPORT(sigma1);
	REPORT(sigma2);

	return nll;
}
