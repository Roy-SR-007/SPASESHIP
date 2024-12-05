# SPASESHIP

## SPArSe Estimation using SHrInkage Priors (SPASESHIP)

* Course Instructor: Dr. Anirban Bhattacharya (TAMU Statistics)

This repository contains the R codes used in STAT 633 Project Report on: ``Efficient Bayesian Computation using Variational Inference for a Shrinkage Prior``. The details are as follows:

* `Armagan_GDP_Gibbs.R`: performs the data-augmented Gibbs sampler as outlined in (Armagan et al. (2013))[https://www3.stat.sinica.edu.tw/statistica/j23n1/J23N16/J23N16.html], with the GDP prior endowed upon the model regression coefficients.
* `Horseshoe.R`: considers using the `monomvn` package in `R` to compute the posterior mean estimates, by putting a Horseshoe prior over the model regression coefficients, where the prior hierarchy follows from Section 1.1 of [https://academic.oup.com/biomet/article-abstract/97/2/465/219397?redirectedFrom=fulltext](Carvalho et al. (2010)).
* `LAPLACE_BLASSO.R`: once again uses the `monomvn` package in `R` to consider implementing the Bayesian LASSO framework, with the LASSO parameter having a Gamma prior with shape and rate parameters as: 2 and 0.1 respectively.
* `NORMAL_BRIDGE.R`: performs the Bayesian ridge regression, using the `monomvn` package in `R`. The ridge parameter is endowed upon with an inverse-Gamma prior having scale and shape parameters as: 1/1000.
* `SCAD.R`: uses the SCAD penalty from [https://www.jstor.org/stable/3085904](Fan and Li (2001)), which is implemented using the `SIS` package in `R`.
* `TAVIE_GDP.R`: implements the VI (TAVIE) algorithm for the GDP prior.
