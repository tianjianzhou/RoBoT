
library(truncdist)
library(gmp)

dyn.load("fn_RoBoT_MCMC.so")


## Some Setup

expit = function(x){
  return( 1 - 1 / (1 + exp(x)) )
}


logit = function(x){
  return( log(x / (1-x)) )
}


get_theta0 = function(pi0, J){
  
  if(length(pi0) != 1 & length(pi0) != J){
    stop("Length of pi0 must either be 1 or match with the length of n and r.")
  }
  
  theta0 = logit(pi0)
  theta0[theta0 < -10] = -10
  theta0[theta0 > 10] = 10

  if(length(pi0) == 1){
    theta0 = rep(theta0, J)
  }

  return(theta0)

}


## Main MCMC function, which calls the compiled C++ file "fn_RoBoT_MCMC.so"

RoBoT_MCMC = function(n, r, pi0, niter, burnin = 5000, thin = 5, 
  return_MCMC_spls = FALSE, alpha = 2){
  
  #################################################################
  # INPUT: 
  # n: J length vector, total number of patients in each group
  # r: J length vector, number of responses in each group
  # niter: number of posterior samples we want 
  #################################################################

  J = length(n)
  
  if(length(pi0) != 1 & length(pi0) != J){
    stop("Length of pi0 must either be 1 or match with the length of n and r.")
  }

  if(length(pi0) == 1){
    pi0 = rep(pi0, J)
  }

  
  #################################################################
  # Set hyperparameters
  #################################################################
  
  # alpha = 2
  pi_prior_shape1 = 0.01
  pi_prior_shape2 = 0.01
  mu_prior_mean = 0
  mu_prior_sd = 2
  tau_prior_loc = 0
  tau_prior_scale = 1


  #################################################################
  # Initialize Markov chains
  #################################################################
  
  pi_spls = matrix(0, J, niter)
  theta_spls = matrix(0, J, niter)
  mu_spls = matrix(0, J, niter)
  tau_spls = matrix(0, J, niter)
  
  
  pi_init = r/n
  pi_init[pi_init < 0.01] = 0.01
  pi_init[pi_init > 0.99] = 0.99
  pi_spls[ , 1] = pi_init
  theta_spls[ , 1] = rep(0, J)
  mu_spls[ , 1] = rep(0, J)
  tau_spls[ , 1] = rep(1, J)

  #################################################################
  # MCMC burnin
  #################################################################
  
  output_C = .C("RoBoT_MCMC", pi = as.double(pi_spls[ , 1]),
                theta = as.double(theta_spls[ , 1]),
                mu = as.double(mu_spls[ , 1]),
                tau = as.double(tau_spls[ , 1]),
                alpha = as.double(alpha),
                J = as.integer(J), 
                n = as.double(n), 
                r = as.double(r), 
                pi0 = as.double(pi0),
                pi_prior_shape1 = as.double(pi_prior_shape1), 
                pi_prior_shape2 = as.double(pi_prior_shape2),
                mu_prior_mean = as.double(mu_prior_mean),
                mu_prior_sd = as.double(mu_prior_sd),
                tau_prior_loc = as.double(tau_prior_loc),
                tau_prior_scale = as.double(tau_prior_scale),
                m = as.integer(5),
                niter = as.integer(burnin))
  
  pi_spls[ , 1] = output_C$pi
  theta_spls[ , 1] = output_C$theta
  mu_spls[ , 1] = output_C$mu
  tau_spls[ , 1] = output_C$tau
  
  #################################################################
  # MCMC iterations
  #################################################################

  for(i in 2:niter){
    
    pi_spls[ , i] = pi_spls[ , i-1]
    theta_spls[ , i] = theta_spls[ , i-1]
    mu_spls[ , i] = mu_spls[ , i-1]
    tau_spls[ , i] = tau_spls[ , i-1]

    output_C = .C("RoBoT_MCMC", pi = as.double(pi_spls[ , i]),
                  theta = as.double(theta_spls[ , i]),
                  mu = as.double(mu_spls[ , i]),
                  tau = as.double(tau_spls[ , i]),
                  alpha = as.double(alpha),
                  J = as.integer(J), 
                  n = as.double(n), 
                  r = as.double(r), 
                  pi0 = as.double(pi0),
                  pi_prior_shape1 = as.double(pi_prior_shape1), 
                  pi_prior_shape2 = as.double(pi_prior_shape2),
                  mu_prior_mean = as.double(mu_prior_mean),
                  mu_prior_sd = as.double(mu_prior_sd),
                  tau_prior_loc = as.double(tau_prior_loc),
                  tau_prior_scale = as.double(tau_prior_scale),
                  m = as.integer(5),
                  niter = as.integer(thin))

    pi_spls[ , i] = output_C$pi
    theta_spls[ , i] = output_C$theta
    mu_spls[ , i] = output_C$mu
    tau_spls[ , i] = output_C$tau
  
  } 
  
  if(return_MCMC_spls){
    
    MCMC_spls = list()
    MCMC_spls$pi_spls = pi_spls
    MCMC_spls$theta_spls = theta_spls
    MCMC_spls$mu_spls = mu_spls
    MCMC_spls$tau_spls = tau_spls
    return(MCMC_spls)

  } else {

    return(pi_spls)
    
  }

}




