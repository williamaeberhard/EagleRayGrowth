// Fabens (1965) estimation with random effect on Linf | v0.4.3
// Gaussian likelihood assumed (original Fabens only based on moments)
// frequentist only, no priors
// * model length at (multiple) recap given observed length at cap and
//   corresponding deltaT's
// * allow for arbitrary nb of individuals (n) and arbitrary nb of measurements
//   (cap and multiple recap) for each individual, thus Lmat has max(n_i) cols
//   with NAs after nbmeas(i) cols for row/indiv i
// * reference remains capture so we don't need to involve age (relative to
//   birth or T0) as yet another random effect
#include <TMB.hpp>

// template <class Type> 
// Type sqr(Type x){
// 	return x*x;
// }

// template <class Type>
// Type neglogdunif(Type x, Type a, Type b, int n){
// 	// neg log density of unif(a,b), a<b
// 	// n = sample size, so that largecst dominates other negloglik contrib
// 	Type largecst = Type(100.0)*n+Type(100.0);
// 	Type halfdensity = Type(0.5)*log(b-a);
// 	Type res1 = CppAD::CondExpGt(x, a, halfdensity, largecst);
// 	Type res2 = CppAD::CondExpLt(x, b, halfdensity, largecst);
// 	return res1+res2; // neg log density if a<x<b
// }

template<class Type>
Type objective_function<Type>::operator() () {

	//--------------------------------------------------------------------------
	// Inputs
	//--------------------------------------------------------------------------

	// Data
	DATA_MATRIX(Lmat); // lengths at cap and recap, dim (n x max(n_i))
	DATA_IVECTOR(nbmeas); // nb meas per indiv incl cap (so all >=2), dim n
	DATA_MATRIX(deltaTmat); // times at liberty, dim (n x (max(n_i)-1))

	// Parameters
	PARAMETER(logmuinf); // mean of Gaussian randef in vB Linf
	PARAMETER(logsigmainf); // sd of Gaussian randef in vB Linf
	PARAMETER(logK); // vB growth rate
	PARAMETER(logsigma); // sd of Gaussian meas err, iid for all recap lengths

	// Random effects
	PARAMETER_VECTOR(logLinf); // iid indiv-specific vB Linf, dim n

	// Hyperparameters
	// DATA_VECTOR(hp_muinf); // hyperparam for prior on muinf, dim 2
	// DATA_VECTOR(hp_sigmainf); // hyperparam for prior on sigmainf, dim 2
	// DATA_VECTOR(hp_K); // hyperparam for prior on K, dim 2
	// DATA_VECTOR(hp_sigmaeps); // hyperparam for prior on meas err sigma, dim 2
	// // ^ hp have different meaning depending on priordist:
	// //    - priordist=0: unused
	// //    - priordist=1: lower and upper bound of uniform density
	// //    - priordist=2: mean and sd of Gaussian density
	// //    - priordist=3: mean and sd (log scale) of lognormal density

	// Misc
	// DATA_IVECTOR(priordist); // code for (Linf,K,sigma), dim 6
	// // ^ 0 = no prior, 1 = uniform, 2 = Gaussian, 3 = lognormal


	//--------------------------------------------------------------------------
	// Setup, procedures and init
	//--------------------------------------------------------------------------

	int n = Lmat.rows(); // nb indiv (indep fish)

	Type muinf = exp(logmuinf);
	Type sigmainf = exp(logsigmainf);
	Type K = exp(logK);
	Type sigma = exp(logsigma);

	vector<Type> Linf = exp(logLinf);

	Type nll = 0.0; // init neg loglik


	//--------------------------------------------------------------------------
	// Priors
	//--------------------------------------------------------------------------

	// if (priordist(0)==1){ // uniform prior on muinf
	// 	nll += neglogdunif(muinf, hp_muinf(0), hp_muinf(1), n);
	// } else if (priordist(0)==2){ // Gaussian prior on muinf
	// 	nll -= dnorm(muinf, hp_muinf(0), hp_muinf(1), true);
	// } else if (priordist(0)==3){ // lognormal prior on muinf
	// 	nll -= dnorm(logmuinf, hp_muinf(0), hp_muinf(1), true) - logmuinf;
	// 	// ^ lognormal dist on exp/ori scale but hp on log scale
	// } // else no prior on muinf

	// if (priordist(1)==1){ // uniform prior on sigmainf
	// 	nll += neglogdunif(sigmainf, hp_sigmainf(0), hp_sigmainf(1), n);
	// } else if (priordist(1)==2){ // Gaussian prior on sigmainf
	// 	nll -= dnorm(sigmainf, hp_sigmainf(0), hp_sigmainf(1), true);
	// } else if (priordist(1)==3){ // lognormal prior on sigmainf
	// 	nll -= dnorm(logsigmainf,hp_sigmainf(0),hp_sigmainf(1),true)-logsigmainf;
	// 	// ^ lognormal dist on exp/ori scale but hp on log scale
	// } // else no prior on sigmainf

	// if (priordist(2)==1){ // uniform prior on K
	// 	nll += neglogdunif(K, hp_K(0), hp_K(1), n);
	// } else if (priordist(2)==2){ // Gaussian prior on K
	// 	nll -= dnorm(K, hp_K(0), hp_K(1), true);
	// } else if (priordist(2)==3){ // lognormal prior on K
	// 	nll -= dnorm(logK, hp_K(0), hp_K(1), true) - logK;
	// 	// ^ lognormal dist on exp/ori scale but hp on log scale
	// } // else no prior on K

	// if (priordist(5)==1){ // uniform prior on sigmaeps
	// 	nll += neglogdunif(sigmaeps, hp_sigmaeps(0), hp_sigmaeps(1), n);
	// } else if (priordist(5)==2){ // Gaussian prior on sigmaeps
	// 	nll -= dnorm(sigmaeps, hp_sigmaeps(0), hp_sigmaeps(1), true);
	// } else if (priordist(5)==3){ // lognormal prior on sigmaeps
	// 	nll -= dnorm(logsigmaeps,hp_sigmaeps(0),hp_sigmaeps(1),true)-logsigmaeps;
	// 	// ^ lognormal log-pdf evaluated at param on exp scale
	// } // else no prior on sigmaeps


	//--------------------------------------------------------------------------
	// Random effects
	//--------------------------------------------------------------------------

	for (int i = 0; i < n; i++){
		nll -= dnorm(Linf(i), muinf, sigmainf, true);
		// ^ iid Gaussian randef on scale of Linf
	}
	

	//--------------------------------------------------------------------------
	// Observation equations
	//--------------------------------------------------------------------------

	// cond on Linf, L1 and L2 are indep and ref is T0 with length=0
	for (int i = 0; i < n; i++){
		Type Lcap = Lmat(i,0); // row 1 of Lmat has cap upon which we condition
		for (int j = 1; j < nbmeas(i); j++){ // nbmeas(i)>=2 for all i
			Type meanL = Linf(i)-(Linf(i)-Lcap)*exp(-K*deltaTmat(i,j-1));
			// ^ vB at recap_j, cap is ref
			nll -= dnorm(Lmat(i,j), meanL, sigma, true);
			// ^ iid Gaussian meas err
		}
	}


	//--------------------------------------------------------------------------
	// Outputs
	//--------------------------------------------------------------------------

	REPORT(muinf); // mean Gaussian randef on Linf, iid
	REPORT(sigmainf); // sd Gaussian randef on Linf, iid
	REPORT(K); // vB growth param, constant across individuals
	REPORT(sigma); // sd additive iid Gaussian meas err, both cap and recap

	REPORT(Linf); // randeff

	// only REPORT for now, ADREPORT for later when needed

	return nll;
}
