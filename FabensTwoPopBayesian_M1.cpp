// Fabens (1965) estimation for two populations, Bayesian | v0.4.3
// Gaussian likelihood assumed (original Fabens only based on moments)
// * model length at recap L2, given observed T1, T2 and L1
// * priors specified by priordist1 and priordist2 arguments for (Linf,K,sigma)
//   for both pop 1 and 2, respectively.
// * M1 model: (Linf,K,sigma) for both pop freely estimated, 6 param overall

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
	PARAMETER(logLinf1); // vB Linf pop 1
	PARAMETER(logK1); // vB growth rate pop 1
	PARAMETER(logsigma1); // additive Gaussian error sd pop 1
	PARAMETER(logLinf2); // vB Linf pop 2
	PARAMETER(logK2); // vB growth rate pop 2
	PARAMETER(logsigma2); // additive Gaussian error sd pop2
	// PARAMETER(deltaLinf); // diff vB Linf, pop 2 - pop 1
	// PARAMETER(deltaK); // diff vB growth rate, pop 2 - pop 1

	// Random effects: none

	// Hyperparameters
	DATA_VECTOR(hp_Linf1); // hyperparam for prior on Linf pop 1, dim 2
	DATA_VECTOR(hp_K1); // hyperparam for prior on K pop 1, dim 2
	DATA_VECTOR(hp_sigma1); // hyperparam for prior on sigma pop 1, dim 2
	DATA_VECTOR(hp_Linf2); // hyperparam for prior on Linf pop 2, dim 2
	DATA_VECTOR(hp_K2); // hyperparam for prior on K pop 2, dim 2
	DATA_VECTOR(hp_sigma2); // hyperparam for prior on sigma pop 2, dim 2
	// ^ hp have different meaning depending on priordist:
	//   - priordist=0: unused
	//   - priordist=1: lower and upper bound of uniform density
	//   - priordist=2: mean and sd of Gaussian density
	//   - priordist=3: mean and sd of Gaussian density on log scale (lognormal)
	// DATA_VECTOR(hp_deltaLinf); // uniform prior lbub on deltaLinf, dim 2
	// DATA_VECTOR(hp_deltaK); // uniform prior lbub on deltaK, dim 2

	// Misc
	DATA_IVECTOR(priordist1); // code for (Linf,K,sigma) pop 1, dim 3
	DATA_IVECTOR(priordist2); // code for (Linf,K,sigma) pop 1, dim 3
	// ^ 0 = no prior, 1 = uniform, 2 = Gaussian, 3 = lognormal


	//--------------------------------------------------------------------------
	// Setup, procedures and init
	//--------------------------------------------------------------------------

	int n1 = L1_1.size(); // i = 1, ..., n1
	int n2 = L1_2.size(); // i = 1, ..., n2

	Type Linf1 = exp(logLinf1);   // pop 1
	Type K1 = exp(logK1);         // pop 1
	Type sigma1 = exp(logsigma1); // pop 1

	Type Linf2 = exp(logLinf2);   // pop 2
	Type K2 = exp(logK2);         // pop 2
	Type sigma2 = exp(logsigma2); // pop 2
	
	Type nll = 0.0; // init neg loglik


	//--------------------------------------------------------------------------
	// Priors
	//--------------------------------------------------------------------------

	// pop 1
	if (priordist1(0)==1){ // uniform prior on Linf1
		nll += neglogdunif(Linf1, hp_Linf1(0), hp_Linf1(1), n1);
	} else if (priordist1(0)==2){ // Gaussian prior on Linf1
		nll -= dnorm(Linf1, hp_Linf1(0), hp_Linf1(1), true);
	} else if (priordist1(0)==3){ // lognormal prior on Linf1
		nll -= dnorm(logLinf1, hp_Linf1(0), hp_Linf1(1), true) - logLinf1;
		// ^ lognormal log-pdf evaluated at param on log scale
	} // else no prior on Linf1

	if (priordist1(1)==1){ // uniform prior on K1
		nll += neglogdunif(K1, hp_K1(0), hp_K1(1), n1);
	} else if (priordist1(1)==2){ // Gaussian prior on K1
		nll -= dnorm(K1, hp_K1(0), hp_K1(1), true);
	} else if (priordist1(1)==3){ // lognormal prior on K1
		nll -= dnorm(logK1, hp_K1(0), hp_K1(1), true) - logK1;
		// ^ lognormal log-pdf evaluated at param on log scale
	} // else no prior on K1

	if (priordist1(2)==1){ // uniform prior on sigma1
		nll += neglogdunif(sigma1, hp_sigma1(0), hp_sigma1(1), n1);
	} else if (priordist1(2)==2){ // Gaussian prior on sigma1
		nll -= dnorm(sigma1, hp_sigma1(0), hp_sigma1(1), true);
	} else if (priordist1(2)==3){ // lognormal prior on sigma1
		nll -= dnorm(logsigma1, hp_sigma1(0), hp_sigma1(1), true) - logsigma1;
		// ^ lognormal log-pdf evaluated at param on log scale
	} // else no prior on sigma1


	// pop 2
	if (priordist2(0)==1){ // uniform prior on Linf2
		nll += neglogdunif(Linf2, hp_Linf2(0), hp_Linf2(1), n2);
	} else if (priordist2(0)==2){ // Gaussian prior on Linf2
		nll -= dnorm(Linf2, hp_Linf2(0), hp_Linf2(1), true);
	} else if (priordist2(0)==3){ // lognormal prior on Linf2
		nll -= dnorm(logLinf2, hp_Linf2(0), hp_Linf2(1), true) - logLinf2;
		// ^ lognormal log-pdf evaluated at param on log scale
	} // else no prior on Linf2

	if (priordist2(1)==1){ // uniform prior on K2
		nll += neglogdunif(K2, hp_K2(0), hp_K2(1), n2);
	} else if (priordist2(1)==2){ // Gaussian prior on K2
		nll -= dnorm(K2, hp_K2(0), hp_K2(1), true);
	} else if (priordist2(1)==3){ // lognormal prior on K2
		nll -= dnorm(logK2, hp_K2(0), hp_K2(1), true) - logK2;
		// ^ lognormal log-pdf evaluated at param on log scale
	} // else no prior on K2

	if (priordist2(2)==1){ // uniform prior on sigma2
		nll += neglogdunif(sigma2, hp_sigma2(0), hp_sigma2(1), n2);
	} else if (priordist2(2)==2){ // Gaussian prior on sigma2
		nll -= dnorm(sigma2, hp_sigma2(0), hp_sigma2(1), true);
	} else if (priordist2(2)==3){ // lognormal prior on sigma2
		nll -= dnorm(logsigma2, hp_sigma2(0), hp_sigma2(1), true) - logsigma2;
		// ^ lognormal log-pdf evaluated at param on log scale
	} // else no prior on sigma2


	//--------------------------------------------------------------------------
	// Random effects
	//--------------------------------------------------------------------------

	// none
	

	//--------------------------------------------------------------------------
	// Observation equations
	//--------------------------------------------------------------------------

	for (int i = 0; i < n1; i++){ // pop 1
		Type meanL2 = Linf1 - (Linf1-L1_1(i))*exp(-K1*deltaT_1(i)); // vB pop 1
		nll -= dnorm(L2_1(i), meanL2, sigma1, true); // Gaussian lkhd
	}

	for (int i = 0; i < n2; i++){ // pop 2
		Type meanL2 = Linf2 - (Linf2-L1_2(i))*exp(-K2*deltaT_2(i)); // vB pop 2
		nll -= dnorm(L2_2(i), meanL2, sigma2, true); // Gaussian lkhd
	}


	//--------------------------------------------------------------------------
	// Outputs
	//--------------------------------------------------------------------------

	REPORT(Linf1);
	REPORT(K1);
	REPORT(sigma1);
	REPORT(Linf2);
	REPORT(K2);
	REPORT(sigma2);

	// ADREPORT(Linf);
	// ADREPORT(K);
	// ADREPORT(sigma1);
	// ADREPORT(deltaLinf);
	// ADREPORT(deltaK);
	// ADREPORT(sigma2);

	return nll;
}
