# RoBoT
This R code is the accompanying software for the paper "RoBoT: A Robust Bayesian Hypothesis Testing Method for Basket Trials" by Tianjian Zhou and Yuan Ji, published in *Biostatistics*. [[Link](https://doi.org/10.1093/biostatistics/kxaa005)]

Because of popular demand, the source C code for MCMC sampling has been made publicly available.

## Details
- The file `fn_RoBoT_MCMC.cpp` contains the source C code (to be called by R) for MCMC sampling for the RoBoT model. 
- The file `fn_RoBoT.R` contains the R functions for running RoBoT. 
- The file `main_RoBoT.R` is the main R file to run RoBoT on a dataset.

To replicate the results reported in the manuscript, 
- Download all three files.
- Open Terminal, run `R CMD SHLIB fn_RoBoT_MCMC.cpp` to generate `fn_RoBoT_MCMC.so` (Mac) or `fn_RoBoT_MCMC.dll` (Windows), the compiled C code for MCMC sampling for the RoBoT model.
- Open `main_RoBoT.R` and execute the code line by line. For Windows users, you need to replace the line `dyn.load("fn_RoBoT_MCMC.so")` in `fn_RoBoT.R` by `dyn.load("fn_RoBoT_MCMC.dll")` for the code to operate. Everything else can be kept as is.
- The real data for the two case studies are also included in the `main_RoBoT.R` file.

## Optional
- The file `fn_RoBoT_MCMC.so` contains the compiled C code (for Mac) for MCMC sampling for the RoBoT model.
- The file `fn_RoBoT_MCMC.dll` contains the compiled C code (for Windows) for MCMC sampling for the RoBoT model.
- Should you have trouble compiling the C code yourself, you may download one of the compiled C codes and run `main_RoBoT.R` directly. Again, for Windows users, you need to replace the line `dyn.load("fn_RoBoT_MCMC.so")` in `fn_RoBoT.R` by `dyn.load("fn_RoBoT_MCMC.dll")` for the code to operate.

Should you have any questions, please contact Tianjian Zhou at tianjian.zhou@colostate.edu.

