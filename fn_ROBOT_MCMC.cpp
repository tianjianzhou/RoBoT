

#include "stdio.h"
#include "R.h"
#include "Rmath.h"
#include "assert.h"
#include "math.h"
#include "stdlib.h"
#include "string.h"


extern "C"{
    
/////////////////////////////////////////////////////////////////
// Functions for MCMC updates
/////////////////////////////////////////////////////////////////



void update_pi(double *pi, double *theta, int J,
  double pi_prior_shape1, double pi_prior_shape2,
  double *n, double *r, double *pi0){

  // Update pi element-wise
  
  int j;
  double lambda_j, pi_j, gamma_j;
  double f1, f2, prob1, prob2, prob;
   
  for(j = 0; j < J; j++){
    
    gamma_j = 1.0 - (1.0 / (1.0 + exp(theta[j])));
    
    f1 = pbeta(pi0[j], r[j] + pi_prior_shape1, n[j] - r[j] + pi_prior_shape2, 1, 0);
    f2 = pbeta(pi0[j], pi_prior_shape1, pi_prior_shape2, 1, 0);
    prob1 = (1 - f1) / (1 - f2);
    prob2 = f1 / f2;
    prob = (prob1 * gamma_j) / (prob1 * gamma_j + prob2 * (1 - gamma_j));

    GetRNGstate();
    lambda_j = rbinom(1, prob);
    PutRNGstate();
    
    // change this to inverse-CDF method??
    // using rejection sampling
    GetRNGstate();
    pi_j = rbeta(r[j] + pi_prior_shape1, n[j] - r[j] + pi_prior_shape2);
    PutRNGstate();

    if(lambda_j > 0.5) {
      
      while(pi_j <= pi0[j]) {
        GetRNGstate();
        pi_j = rbeta(r[j] + pi_prior_shape1, n[j] - r[j] + pi_prior_shape2);
        PutRNGstate();
      }

    }
    else {
      
      while(pi_j > pi0[j]) {
        GetRNGstate();
        pi_j = rbeta(r[j] + pi_prior_shape1, n[j] - r[j] + pi_prior_shape2);
        PutRNGstate();
      }

    }

    pi[j] = pi_j;
    
  }

  return;
}







void update_theta(double *pi, double *theta, double *mu, double *tau,
  int J, double *pi0) {
  
  int j;
  
  double theta_j_cur, theta_j_pro, loglik_cur, loglik_pro, u;
  
  for(j = 0; j < J; j++) {
    
    theta_j_cur = theta[j];
    
    GetRNGstate();
    theta_j_pro = rnorm(mu[j], tau[j]);
    PutRNGstate();
    
    if(pi[j] > pi0[j]) {
      // log1pexp(x) = log(1+exp(x))
      loglik_cur = theta_j_cur - log1pexp(theta_j_cur);
      loglik_pro = theta_j_pro - log1pexp(theta_j_pro);
    }
    else {
      loglik_cur = - log1pexp(theta_j_cur);
      loglik_pro = - log1pexp(theta_j_pro);
    }
    
    GetRNGstate();
    u = runif(0.0, 1.0);
    PutRNGstate();
    
    if(log(u) < loglik_pro - loglik_cur){
      theta[j] = theta_j_pro;
    }
    
  }
  
  return;
 
 
}








void update_mu_tau(double *theta, double *mu, double *tau,
  double *alpha, int J,
  double mu_prior_mean, double mu_prior_sd,
  double tau_prior_loc, double tau_prior_scale,
  int m){
  
  
  int is_unique_j;
  int j, j1, m1;
  
  double *mu_pro, *tau_pro; 
  mu_pro = (double *)calloc(m, sizeof(double));
  tau_pro = (double *)calloc(m, sizeof(double));

  double *prob_c;
  prob_c = (double *)calloc(J + m, sizeof(double));

  double mu_j, tau_j, max_logprob_c, sum_prob_c, tau_pro_temp;
  double u;

  /////////////////////////////////////////////////////////////////////
  // Step 1: Update class label of object j
  
  for(j = 0; j < J; j++){
    
    is_unique_j = 1;
    max_logprob_c = -INFINITY;
    sum_prob_c = 0;
    
    mu_j = mu[j];
    tau_j = tau[j];
    
    for(m1 = 0; m1 < m; m1++) {
      mu_pro[m1] = 0;
      tau_pro[m1] = 0;
    }
    
    for(j1 = 0; j1 < J + m; j1++) {
      prob_c[j1] = 0;
    }
    
    // Start sampling
    // first check if j is a singleton
    for(j1 = 0; j1 < J; j1++){
      
      if(j1 != j && mu[j1] == mu_j && tau[j1] == tau_j) {
        is_unique_j = 0;
        break;
      }

    }
    
    
    for(j1 = 0; j1 < J; j1++) {
        
      if(j1 != j){
        
        prob_c[j1] = dnorm(theta[j], mu[j], tau[j], 1);
        
        if(max_logprob_c < prob_c[j1]) max_logprob_c = prob_c[j1];
        
      }
    
    }
    
    // If j is not a singleton, draw m values from prior
    // If j is a singleton (is_unique_j = 1), draw m - 1 values from prior 
    //    and keep the current value
    if(is_unique_j == 1) {
        
      mu_pro[0] = mu[j];
      tau_pro[0] = tau[j];
      prob_c[J] = dnorm(theta[j], mu_pro[0], tau_pro[0], 1);
      prob_c[J] = prob_c[J] + log((*alpha)) - log((double) m);
      
      if(max_logprob_c < prob_c[J]) max_logprob_c = prob_c[J];
    }
      
    for(m1 = is_unique_j; m1 < m; m1++) {
      
      // draw m proposals from prior
      
      // draw mu from normal
      GetRNGstate();
      mu_pro[m1] = rnorm(mu_prior_mean, mu_prior_sd);
      PutRNGstate();
      
      // draw tau from half-Cauchy
      GetRNGstate();
      tau_pro_temp = rcauchy(tau_prior_loc, tau_prior_scale);
      PutRNGstate();
      
      while(tau_pro_temp <= 0) {
        GetRNGstate();
        tau_pro_temp = rcauchy(tau_prior_loc, tau_prior_scale);
        PutRNGstate();
      }

      tau_pro[m1] = tau_pro_temp;
      
      
      prob_c[J + m1] = dnorm(theta[j], mu_pro[m1], tau_pro[m1], 1);
      prob_c[J + m1] = prob_c[J + m1] + log((*alpha)) - log((double) m);
      
      if(max_logprob_c < prob_c[J + m1]) max_logprob_c = prob_c[J + m1];
    
    }
    
    
    for(j1 = 0; j1 < J + m; j1++) {

      if(j1 != j) {
      
        prob_c[j1] = prob_c[j1] - max_logprob_c;
        prob_c[j1] = exp(prob_c[j1]);
        sum_prob_c = sum_prob_c + prob_c[j1];

      }

    }
    
    prob_c[j] = 0;
    
    
    // sample class label
    GetRNGstate();
    u = runif(0.0, 1.0);
    PutRNGstate(); 
      
    j1 = 0;
    prob_c[j1] = prob_c[j1] / sum_prob_c;
    while(u > prob_c[j1]){
  	  j1++;
  	  prob_c[j1] = prob_c[j1] / sum_prob_c + prob_c[j1 - 1];
    }
    
    // if from an existing cluster
    if(j1 < J) {
      mu[j] = mu[j1];
      tau[j] = tau[j1];
    }
    // if from a new cluster
    else {
      mu[j] = mu_pro[j1 - J];
      tau[j] = tau_pro[j1 - J];
    }
  
  
  }
  
  
  
  
  
  /////////////////////////////////////////////////////////////////////
  // Step 2: Update distinct mu & tau's

  int *mu_tau_updated;
  mu_tau_updated = (int *)calloc(J, sizeof(int));
  
  int n_c;
  double theta_c_sum, norm_prop, logpost_cur, logpost_pro;
  double sd_prop = 0.5;

  double mu_post_mean, mu_post_var;
  double mu_j_new, tau_j_pro, tau_j_new;
  
  for(j = 0; j < J; j++){
    
    // if mu[j] and tau[j] has not been updated yet
    if(mu_tau_updated[j] == 0){
      
      mu_j = mu[j];
      tau_j = tau[j];
      
      GetRNGstate();
      norm_prop = rnorm(0.0, sd_prop);
      PutRNGstate();
      
      tau_j_pro = tau_j * exp(norm_prop);

      n_c = 0;
      theta_c_sum = 0;    // sum theta[j]

      for(j1 = 0; j1 < J; j1++){

        if(mu[j1] == mu_j && tau[j1] == tau_j){
          
          n_c++;
          theta_c_sum = theta_c_sum + theta[j1];

        }

      }
      
      // Sample mu: conjugate Normal
      mu_post_var = 1.0 / ( (1.0 / pow(mu_prior_sd, 2)) + ((double) n_c / pow(tau_j, 2)) );
      mu_post_mean = mu_post_var * ( (mu_prior_mean / pow(mu_prior_sd, 2)) + (theta_c_sum / pow(tau_j, 2)) );
          
      GetRNGstate();
      mu_j_new = rnorm(mu_post_mean, sqrt(mu_post_var));
      PutRNGstate(); 
      
      ///// Update tau, Check this!!!!
      logpost_cur = 0;
      logpost_pro = 0;

      for(j1 = 0; j1 < J; j1++){

        if(mu[j1] == mu_j && tau[j1] == tau_j){
          
          logpost_cur = logpost_cur + dnorm(theta[j1], mu_j_new, tau_j, 1);
          logpost_pro = logpost_pro + dnorm(theta[j1], mu_j_new, tau_j_pro, 1);

        }

      }
      
      // constant for hC cancels out??

      logpost_cur = logpost_cur + dcauchy(tau_j, tau_prior_loc, tau_prior_scale, 1) + log(tau_j);
      logpost_pro = logpost_pro + dcauchy(tau_j_pro, tau_prior_loc, tau_prior_scale, 1) + log(tau_j_pro);
      
      GetRNGstate();
      u = runif(0.0, 1.0);
      PutRNGstate();
      
      if(log(u) < logpost_pro - logpost_cur) {
        tau_j_new = tau_j_pro;
      }
      else {
        tau_j_new = tau_j;
      }
      //// Finish tau
      
      for(j1 = 0; j1 < J; j1++) {
        
        if(mu[j1] == mu_j && tau[j1] == tau_j) {
          
          mu[j1] = mu_j_new;
          tau[j1] = tau_j_new;
          mu_tau_updated[j1] = 1;

        }

      }
      
    }
    

  }

  free(mu_pro);
  free(tau_pro);
  free(prob_c);
  free(mu_tau_updated);

  return;

}








void ROBOT_MCMC(double *pi, double *theta, double *mu, double *tau,
  double *alpha, int *J,
  double *n, double *r, double *pi0,
  double *pi_prior_shape1, double *pi_prior_shape2,
  double *mu_prior_mean, double *mu_prior_sd,
  double *tau_prior_loc, double *tau_prior_scale,
  int *m, int *niter){
  
  // the m here is the D in the supplement of the paper

  int i;

  for(i = 0; i < (*niter); i++){
  
    update_pi(pi, theta, *J, *pi_prior_shape1, *pi_prior_shape2, n, r, pi0);
    
    update_theta(pi, theta, mu, tau, *J, pi0);
    
    update_mu_tau(theta, mu, tau, alpha, *J, *mu_prior_mean, *mu_prior_sd, *tau_prior_loc, *tau_prior_scale, *m);

  }

  return;

}



}



