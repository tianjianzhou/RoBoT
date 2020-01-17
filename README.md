# RoBoT
This R code is the accompanying software for the paper "RoBoT: A Robust Bayesian Hypothesis Testing Method for Basket Trials" by Tianjian Zhou and Yuan Ji

- The file `fn_RoBoT_MCMC.so` contains the compiled C code for running the Markov chain Monte Carlo sampling for the RoBoT model. 
- The file `fn_RoBoT.R` contains the R functions for running RoBoT. 
- The file `main_RoBoT.R` is the main R file to run RoBoT on a dataset.

To replicate the results reported in the manuscript, just download all these three files. Open `main_RoBoT.R` and execute the code line by line.

The real data for the two case studies are also included in the `main_RoBoT.R` file.

Should you have any questions, please contact Tianjian Zhou at tjzhou@uchicago.edu.
