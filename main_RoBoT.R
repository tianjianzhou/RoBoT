
source("fn_RoBoT.R")

########################################################
## Case Study 1: Imatinib Trial
########################################################

## Data
# n: length J vector, J is the number of baskets, and
#    n[j] is the total number of patients in basket j
# r: length J vector, r[j] is the # of responses in basket j
# pi0: length J vector, pi0[j] is the control response rate for basket j

J = 10
n = c(15, 13, 12 ,28, 29, 29, 26, 5, 2, 20)
r = c(2, 0, 1, 6, 7, 3, 5, 1, 0, 3)
pi0 = rep(0.3, J)

## MCMC sampling for the RoBoT model

MCMC_spls = RoBoT_MCMC(n, r, pi0, 
                       niter = 10000, burnin = 50000, thin = 5, 
                       return_MCMC_spls = TRUE, alpha = 2)

# pi_spls: J * niter matrix, posterior samples of pi
pi_spls = MCMC_spls$pi_spls


# PP: posterior probability of the J alternatives
PP = apply(pi_spls > pi0, 1, sum) / dim(pi_spls)[2]
# PM: posterior mean of the response rates in the J baskets
PM = apply(pi_spls, 1, mean)
# PM_lower and PM_upper: posterior CI bounds
PM_lower = apply(pi_spls, 1, function(x){quantile(x, probs = 0.025)})
PM_upper = apply(pi_spls, 1, function(x){quantile(x, probs = 0.975)})

## NOTE: results may be slightly different from what were reported in
##       the paper due to randomness of Monte Carlo sampling

########################################################
## Case Study 2: Vemurafenib Trial
########################################################

## Data
J = 6
n = c(7, 14, 8 ,26, 10, 19)
r = c(2, 6, 1, 1, 0, 8)
pi0 = rep(0.15, J)

## Do the same thing as Case Study 1
