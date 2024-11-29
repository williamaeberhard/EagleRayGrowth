EagleRayGrowth: Growth Estimation for Eagle Rays
------------------------------------------------

These R scripts provide functions for estimating the von Bertalanffy (vB) growth parameters (Linf and K) from mark-recapture data, as well as comparing two populations (e.g. wild and aquarium) in terms of their growth parameters, and fitting a version of the Fabens (1965) model with a random intercept on the Linf parameter. These scripts were developped for the analyses in Boggio-Pasqua et al. (2024), where animal length is represented by disc width (DW).

### Contents

Files contained in this repository:
* EagleRayGrowth_MAIN.r: the main R script to follow, loading all required libraries, compiling the Template Model Builder (TMB) C++ scripts, simulating data for the sake of illustration, and applying all the methods discussed in Boggio-Pasqua et al. (2024);
* EagleRayGrowth_Methods.r: an R script to source to create all the necesary functions required by the estimation methods and tests;
* FabensBayesian.cpp: a C++ script to compile with the R package TMB for the Bayesian fitting of the Fabens (1965) model;
* FabensRandef.cpp: a C++ script to compile for the Fabens (1965) model with random intercept on the Linf parameter when more than one recaptures are available;
* FabensTwoPop.cpp: a C++ script to compile for the likelihood ratio test (LRT) comparing two populations in terms of the vB growth parameters Linf and K;
* FabensTwoPopBayesian_M0.cpp and FabensTwoPopBayesian_M1.cpp: C++ scripts to compile for the Bayes factor comparing two populations in terms of vB growth parameters;
* this README file.


### Version History

This is EagleRayGrowth version 0.4.3. This is the initial release.


### References

Boggio-Pasqua, A., Bassos-Hull, K., Aeberhard, W. H., Hoopes, L. A., Swider, D. A., Wilkinson, K. A., and Dureuil, M. (2024) Whitespotted eagle ray (Aetobatus narinari) age and growth in wild (in-situ) versus aquarium-housed (ex-situ) individuals: Implications for conservation and management. Frontiers in Marine Science 9, 1–18. DOI: 10.3389/fmars.2022.960822

Fabens, A. J. (1965). Properties and Fitting of the von Bertalanffy Growth Curve. Growth, 29, 265-289


