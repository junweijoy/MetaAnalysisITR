
<<<<<<< HEAD:Data_Sparsity.R
##### R CODE FOR Sparse two-stage Bayesian meta-analysis for individualized treatments
=======
##### R CODE FOR

# Sparse two-stage Bayesian meta-analysis for individualized treatments
<<<<<<< HEAD:Data_Sparsity.R
>>>>>>> bdb6bd94e4ad2b7d24ffc2721f163bc2c69c0676:sparse example.R
=======
>>>>>>> bdb6bd94e4ad2b7d24ffc2721f163bc2c69c0676:sparse example.R
#



rm( list = ls() )

library(tidyverse)
library(dplyr)
library(MASS)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)



#' Generate individual-level data
#'
#' @param k the number of sites; k=10 in simulations
#' @param sample_size_site a vector of length k including sample sizes in each site
#' @param heter_level heterogeneity level; It could be 0, 0.1, 0.2, and 0.3 in simulations
#' @param psi_df a list of blip function parameters
#' @param coeff_outcome_df a list of treatment-free function parameters
#' @param sigma_w within-site variance
#' @export
#'



Individual_Data_sim_f <- function( k = 10, sample_size_site, heter_level,
                                   psi_df, coeff_outcome_df, sigma_w )
{
  # Sites and patients
  
  Tab1 <- data.frame( site = rep( seq(1:k) , sample_size_site ) )
  
  Tab1$patient          <- rep( NA , length(Tab1$site) )
  
  Tab1$sample_size_site <- rep( NA , length(Tab1$site) )
  
  
  
  Tab1$heter_level     <- heter_level
  
  
  for( i in 1:k )
  {
    Tab1[Tab1$site == i,]$patient <- seq( 1 , sample_size_site[i] )
    Tab1[Tab1$site == i,]$sample_size_site <- rep( sample_size_site[i] , sample_size_site[i] )
  }
  
  # Covariates x1, x2
  # x1: binary
  # x2: categorical
  
  
  Tab1$x1 <- rep( NA , length(Tab1$site) )
  Tab1$x2 <- rep( NA , length(Tab1$site) )
  Tab1$x3 <- rep( NA , length(Tab1$site) )
  
  
  for ( i in 1:k ) 
  {
    
    if( unique( Tab1[Tab1$site == i,]$site ) %% 3 == 0)
    {
  
      Tab1[Tab1$site == i,]$x1 <- rep(1, sample_size_site[i] )
      s <- rmultinom(sample_size_site[i], 1, c(0, 0.5, 0.5))
      Tab1[Tab1$site ==i,]$x2 <- s[2,]
      Tab1[Tab1$site ==i,]$x3 <- s[3,]
    }
    else if( unique( Tab1[Tab1$site == i,]$site ) %% 3 == 1 )
    {
      Tab1[Tab1$site == i,]$x1 <- rep(0, sample_size_site[i] )
      s <- rmultinom(sample_size_site[i],1, c(0.5, 0, 0.5))
      Tab1[Tab1$site == i,]$x2 <- s[2,]
      Tab1[Tab1$site == i,]$x3 <- s[3,]
      
    }
    else if( unique( Tab1[Tab1$site == i,]$site ) %% 3 == 2)
    {
      Tab1[Tab1$site == i,]$x1 <- rbinom(sample_size_site[i], 1, 0.5)
      s <- rmultinom(sample_size_site[i], 1, c(1/3, 1/3, 1/3))
      Tab1[Tab1$site == i,]$x2 <- s[2,]
      Tab1[Tab1$site == i,]$x3 <- s[3,]
    }
  }
  
  
  
  
  # Treatment 
  
  Tab1$A <- rep( NA , dim(Tab1)[1] )
  
  
  
  Tab1$A      <- rbinom( sum(sample_size_site) , 1, 0.5 )
  
  # Outcome model
  
  
  sigma_b <- heter_level * sigma_w / ( 1 - heter_level )
  
  
  
  var_matrix <-  diag( sigma_b, length( psi_df ) + length( coeff_outcome_df ) )
  
  coeff_outcome_df_common <- coeff_outcome_df
  
  psi_df_common <- psi_df
  
  coeff <- mvrnorm(k, mu = unlist(c( coeff_outcome_df_common, psi_df_common )),
                   Sigma = var_matrix )
  
  coeff_outcome_df <- data.frame(
    beta_0 = rep( coeff[,1], sample_size_site ) ,
    beta_1 = rep( coeff[,2], sample_size_site ),
    beta_2 = rep( coeff[,3], sample_size_site ),
    beta_3 = rep( coeff[,4], sample_size_site )
  )
  
  
  psi_df <- data.frame( psi_0 = rep( coeff[,5], sample_size_site ),
                        psi_1 = rep( coeff[,6], sample_size_site ),
                        psi_2 = rep( coeff[,7], sample_size_site ),
                        psi_3 = rep( coeff[,8], sample_size_site )
  )

  
  Tab1$beta_0 <- coeff_outcome_df$beta_0 
  
  Tab1$beta_1 <- coeff_outcome_df$beta_1
  
  Tab1$beta_2 <- coeff_outcome_df$beta_2
  
  Tab1$beta_3 <- coeff_outcome_df$beta_3
  
  Tab1$psi_0 <- psi_df$psi_0
  
  Tab1$psi_1 <- psi_df$psi_1
  
  Tab1$psi_2 <- psi_df$psi_2
  
  Tab1$psi_3 <- psi_df$psi_3
  
  
  epsilon <- rnorm( sum(sample_size_site), mean = 0, sd = sqrt(sigma_w) )
  
  tf      <- coeff_outcome_df$beta_0 + coeff_outcome_df$beta_1 * Tab1$x1 + coeff_outcome_df$beta_2 * Tab1$x2+
    coeff_outcome_df$beta_3 * Tab1$x3
  
  blip    <- Tab1$A * ( psi_df$psi_0 + psi_df$psi_1 * Tab1$x1 + 
                          psi_df$psi_2 * Tab1$x2 + psi_df$psi_3 * Tab1$x3 )
  
  Y       <- tf + blip + epsilon 
  
  Tab1$epsilon  <- epsilon
  
  Tab1$Y  <- Y
  
  return(Tab1)
  
}


#' Generate site-specific estimates
#'
#' @param Tab1 a data frame of individual-level data
#' @export
#'

Individual_lm_Data_sim_f <- function(Tab1)
{
  
  Tab2  <-  list()
  
  site  <-  unique( Tab1$site )
  
  for (i in 1:length(site) ) 
  {
    
    sub_Tab1  <-  Tab1[Tab1$site == i, ]
    
    Y         <-  sub_Tab1$Y
    
    x1        <-  sub_Tab1$x1
    x2        <-  sub_Tab1$x2
    x3        <-  sub_Tab1$x3
    A         <-  sub_Tab1$A
    
    mod       <-  lm(Y~ (x1 + x2 + x3) * A)
    
    Tab2[[i]] <- list( estimate     = coefficients(mod),
                       covariance   = vcov(mod)
    )
  }
  
  return(Tab2)
}

# -----------------------------------------------------------------
#  stan code/model for two-stage approach
# -----------------------------------------------------------------


#'
#' @param k the number of sites; k=10 in simulations in a sparse data setting
#' @param ppsi the number of blip function parameters; ppsi = 4 in simulations in a sparse data setting
#' @param pbeta the number of treatment-free function parameters; pbeta = 4 in simulations in a sparse data setting
#' @param index1,index2,index3 site ID in three different sparsity scenarios
#' @param estimate matrix of site-specific blip function parameter estimates
#' @param varest matrix of variances associated with site-specific blip function parameter estimates
#' @param priorsdscale scale parameter for half-cauchy prior
#' @param para common blip function parameters
#' @param sdpara between-site standard deviation
#' @export
#'


stancode_twostage<-"
  data{
  
  int<lower=1> k; 
  
  int<lower=1> ppsi;
  
  int<lower=1> index1[4];
  int<lower=1> index2[3];
  int<lower=1> index3[3];
  
  real estimate[k,ppsi]; 
  
  real<lower=0> varest[k, ppsi]; 
  
  real<lower=0> priorsdscale; 
  }
  
  
  parameters{
  real para[ppsi];  
  real<lower=0> sdpara[ppsi];  
  }
  
  
  model{
  
  
  for(i in index1)
{
estimate[i,1] ~ normal( para[1], sqrt(varest[i,1] + sdpara[1]^2 ) );
estimate[i,4] ~ normal( para[4], sqrt(varest[i,4] + sdpara[4]^2 ) );
}

  for(i in index2)
{
estimate[i,1] ~ normal( para[1], sqrt(varest[i,1] + sdpara[1]^2 ) );
estimate[i,2] ~ normal( para[2], sqrt(varest[i,2] + sdpara[2]^2 ) );
estimate[i,3] ~ normal( para[3], sqrt(varest[i,3] + sdpara[3]^2 ) );
estimate[i,4] ~ normal( para[4], sqrt(varest[i,4] + sdpara[4]^2 ) );
}

  for(i in index3)
{
estimate[i,1] ~ normal( para[1] + para[2] + para[4], 
                 sqrt(varest[i,1] + sdpara[1]^2 + sdpara[2]^2 + sdpara[4]^2) );
estimate[i,3] ~ normal( para[3] - para[4], sqrt(varest[i,3] + sdpara[3]^2 + sdpara[4]^2));
}
      
   
   para ~ normal( 0, 100 );
   sdpara ~ cauchy( 0, priorsdscale );
  }
"

stan_twostage<-stan_model(model_code = stancode_twostage)




twostage_fun <- function(Tab2,k=10, ppsi=4,pbeta=4, priorsdscale)
{
  
  estimate <- matrix(0, nrow=k, ncol = ppsi )
  
  
  for(i in 1:k)
    estimate[i,] <- Tab2[[i]]$estimate[(pbeta+1):(pbeta+ppsi)]
  
  
  estimate[is.na(estimate)] <- 0
  
  
  varest <- matrix(0, nrow = k, ncol = ppsi )
  
  for(i in 1:k)
    varest[i,] <- diag(Tab2[[i]]$covariance)[(pbeta+1):(pbeta+ppsi)]
  
  varest[is.na(varest)] <- 0
  
  data_twostage <- list(
    
    ppsi = ppsi,
    
    k = k,
    
    index1 = c(1,4,7,10),
    index2 = c(2,5,8),
    index3 = c(3,6,9),
    
    estimate = estimate,
    varest = varest,
    
    priorsdscale = priorsdscale
  )
  
  twostage <- sampling(stan_twostage,data = data_twostage, chains = 2, iter =2000,
                         control = list(adapt_delta = 0.95,
                                        max_treedepth = 20 ))
  
  return(twostage)
}


# -----------------------------------------------------------------
#  stan code/model for one-stage approach
# -----------------------------------------------------------------

#'
#' @param k the number of sites; k=10 in simulations in a sparse data setting
#' @param N total sample size
#' @param ppsi the number of blip function parameters; ppsi = 4 in simulations in a sparse data setting
#' @param pbeta the number of treatment-free function parameters; pbeta = 4 in simulations in a sparse data setting
#' @param Xmatrix design matrix
#' @param reward patient outcomes
#' @param site vector of site ID for each patient
#' @param priorsdscale scale parameter for half-cauchy prior
#' @param common common treatment-free and blip function parameters
#' @param betweensd between-site standard deviation
#' @param withinsd within-site standard deviation
#' @export
#'



stancode_onestage<-"
  data{
  
  int<lower=1> k; 
  int<lower=1> N; 
  
  int<lower=1> pbeta; 
  int<lower=1> ppsi; 
  
  matrix[N,pbeta+ppsi] Xmatrix;
  
  vector[N] reward; 
  
  int<lower=1> site[N]; 
  
  real<lower=0> priorsdscale; 
  }
  
  
  parameters{
  matrix[k, pbeta+ppsi] eta;
  
  vector[pbeta+ppsi] common; 
  real<lower=0> betweensd[pbeta+ppsi];  
  real<lower=0> withinsd[k]; 
  }
  
  transformed parameters{
  matrix[k, pbeta+ppsi] sitepar;
  
   for(p in 1:k)
  { 
  for(q in 1:(pbeta+ppsi))
   {
   sitepar[p,q] = common[q] + betweensd[q] * eta[p,q];
   }
  }
  }
  
  
  model{
  
  for(i in 1:N)
     reward[i] ~ normal(Xmatrix[i,] * (sitepar[site[i],])', 
     withinsd[site[i]]);
      
  for(p in 1:k)
    for(q in 1:(pbeta + ppsi))
    {
      eta[p,q] ~ normal(0,1);
    }
        
  for(q in 1:(pbeta + ppsi ) )      
   {
     common[q] ~ normal( 0, 100 );
     betweensd[q] ~ cauchy( 0, priorsdscale );
   }
   
   for(p in 1:k)
   withinsd[p] ~ cauchy( 0, priorsdscale );
  }
"

stan_onestage<-stan_model(model_code = stancode_onestage)





onestage_fun <- function(Tab1,k=10,N,pbeta=4, ppsi = 4,priorsdscale)
{
  
  data_onestage <- list(
    k = k,
    N = N,
    
    pbeta = pbeta,
    ppsi = ppsi,
    
    Xmatrix = model.matrix(~(x1+x2+x3) * A, data = Tab1),
    
    reward = Tab1$Y,
    
    site = Tab1$site,
    
    priorsdscale = priorsdscale
    
  )
  
  onestage <- sampling(stan_onestage, data = data_onestage, chains=2, iter = 2000,
                  control = list(adapt_delta = 0.95,
                                 max_treedepth = 20))
  return(onestage)
}



###########################################################
#           Example for one simulated dataset
###########################################################

k <- 10
N <- 2000

sample_size_site <- c( 200*0.06 , 200*0.07 , 200*0.08 , 200*0.09 , 200*0.10 , 200*0.10 , 200*0.11 , 200*0.12 , 200*0.13 , 200*0.14 )*10 # if large sample

psi_df <- data.frame(psi_0 = 1, 
                     psi_1 = 1, 
                     psi_2 = -2.5, 
                     psi_3 = 2)
coeff_outcome_df <- data.frame(   beta_0 = 4,
                                  beta_1 = 1,
                                  beta_2 = 1,
                                  beta_3 = -1 )

heter_level <- 0.3

sigma_w <- 0.25

Tab1 <-Individual_Data_sim_f( k = k, 
                              sample_size_site = sample_size_site ,
                              heter_level = heter_level,
                              psi_df = psi_df,
                              coeff_outcome_df = coeff_outcome_df,
                              sigma_w = sigma_w )


mod_onestage <- onestage_fun(Tab1 = Tab1,
                             k = 10,
                             N = 2000,
                             pbeta = 4,
                             ppsi = 4,
                             priorsdscale = 1)

Tab2 <- Individual_lm_Data_sim_f( Tab1 = Tab1 )

mod_twostage <- twostage_fun(Tab2 = Tab2,
                             k = 10,
                             ppsi = 4,
                             pbeta = 4,
                             priorsdscale = 1)





