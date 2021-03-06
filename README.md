# RoBoT
This R code is the accompanying software for the paper "RoBoT: A Robust Bayesian Hypothesis Testing Method for Basket Trials" by Tianjian Zhou and Yuan Ji, published in *Biostatistics*. [[Link](https://doi.org/10.1093/biostatistics/kxaa005)]

## Details
- The file `fn_RoBoT_MCMC.so` contains the compiled C code for running the Markov chain Monte Carlo sampling for the RoBoT model. 
- The file `fn_RoBoT.R` contains the R functions for running RoBoT. 
- The file `main_RoBoT.R` is the main R file to run RoBoT on a dataset.

To replicate the results reported in the manuscript, just download all these three files. Open `main_RoBoT.R` and execute the code line by line.

The real data for the two case studies are also included in the `main_RoBoT.R` file.

The code is supposed to be run on Mac OS and was tested on Mac OS High Sierra version 10.13.6.

Should you have any questions, please contact Tianjian Zhou at tjzhou@uchicago.edu.
