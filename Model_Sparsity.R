
##### R CODE FOR Sparse two-stage Bayesian meta-analysis for individualized treatments
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
#' @param num_covariate the number of covariates; 10 or 20 in simulations
#' @param sample_size_site a vector of length k including sample sizes in each site
#' @param heter_level heterogeneity level; It could be 0, 0.1, 0.2, and 0.3 in simulations
#' @param psi_df a list of blip function parameters
#' @param coeff_outcome_df a list of treatment-free function parameters
#' @param sigma_w within-site variance
#' @export
#'

Individual_Data_sim_f <- function( k , num_covariate, sample_size_site ,  heter_level, 
                                   psi_df , coeff_outcome_df, sigma_w )
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
  
  # Covariates x1, x2, and x3 have tailoring effects treatment assignment, i.e., non-zero treatment-covariate interactions. 
  # x1, x3: continuous
  # x2: binary
  # all remaining covariates: normally distributed, and no tailoring effects
  
  
  
  Tab1$x1 <- rep( NA , length(Tab1$site) )
  Tab1$x2 <- rep( NA , length(Tab1$site) )
  Tab1$x3 <- rep( NA , length(Tab1$site) )
  other_covariate <- matrix(rnorm(length(Tab1$site) * (num_covariate - 3), 0, 1),
                            ncol = num_covariate - 3)
  colnames(other_covariate) <- paste0("x", 4:num_covariate)
  Tab1 <- cbind(Tab1,other_covariate)
  
  for ( i in 1:k ) 
  {
    
    if( unique( Tab1[Tab1$site == i,]$site ) %% 3 == 0)
    {
      Tab1[Tab1$site == i,]$x1 <- rnorm(sample_size_site[i],5,1)
      Tab1[Tab1$site == i,]$x2 <- ifelse( runif( sample_size_site[i] ) < 0.50, 1, 0 )
      Tab1[Tab1$site == i,]$x3 <- rexp(sample_size_site[i], 1)
    }
    else if( unique( Tab1[Tab1$site == i,]$site ) %% 3 == 1 )
    {
      Tab1[Tab1$site == i,]$x1 <- 6 * rbeta(sample_size_site[i], 4, 4)+2
      Tab1[Tab1$site == i,]$x2 <- ifelse( runif( sample_size_site[i] ) < 0.30, 1, 0 )
      Tab1[Tab1$site == i,]$x3 <- rexp(sample_size_site[i],1/0.6)
    }
    else if( unique( Tab1[Tab1$site == i,]$site ) %% 3 == 2)
    {
      Tab1[Tab1$site == i,]$x1 <- runif(sample_size_site[i], 2, 8)
      Tab1[Tab1$site == i,]$x2 <- ifelse( runif( sample_size_site[i] ) < 0.70, 1, 0 )
      Tab1[Tab1$site == i,]$x3 <- rexp(sample_size_site[i],1/1.4)
    }
  }
  
  
  # Treatment 
  
  Tab1$A <- rep( NA , dim(Tab1)[1] )

  
  Tab1$A      <- rbinom( sum(sample_size_site) , 1, 0.5 )

  
  # Generate site-specific parameters
  sigma_b <- heter_level * sigma_w / ( 1 - heter_level )
  
  var_matrix <-  diag( sigma_b, length( psi_df ) + length( coeff_outcome_df ) )
  
  coeff_outcome_df_common <- coeff_outcome_df
  
  psi_df_common <- psi_df
  
  
  coeff <- mvrnorm(k, mu = unlist(c( coeff_outcome_df_common, psi_df_common )),
                   Sigma = var_matrix )
  
  coeff_outcome_df <- as.data.frame(
    sapply(1:(num_covariate+1), function(index)
    {
      rep( coeff[,index], sample_size_site)
    }
    )
  )
  colnames(coeff_outcome_df) <- paste0("beta_", 0:num_covariate)
  rep( coeff[,11], sample_size_site )
  
  
  psi_df <- as.data.frame(
    sapply((num_covariate+2):(2*num_covariate+2), function(index)
    {
      rep( coeff[,index], sample_size_site)
    }
    )
  )
  
  colnames(psi_df) <- paste0("psi_", 0:num_covariate)
  
  
  Tab1 <- cbind(Tab1, coeff_outcome_df, psi_df)
  
  
  
  # Generate outcomes
  epsilon <- rnorm( sum(sample_size_site), mean = 0, sd = sqrt(sigma_w) )
  
  
  tf <- apply(coeff_outcome_df * cbind(1, Tab1[,6:(5+num_covariate)]), 1, sum)
  
  
  blip <- Tab1$A * apply(
    psi_df * cbind(1, Tab1[,6:(5+num_covariate)]), 1, sum
  )
  
  Y       <- tf + blip + epsilon 
  
  Tab1$epsilon  <- epsilon
  
  Tab1$Y  <- Y
  
  
  return(Tab1)
  
}


#' Generate site-specific estimates
#'
#' @param Tab1 a data frame of all individual-level data
#' @param num_covariate number of available covariates
#' @export
#'

Individual_lm_Data_sim_f <- function(Tab1,num_covariate)
{
  
  Tab2  <-  list()
  
  site  <-  unique( Tab1$site )
  
  for (i in 1:length(site) ) 
  {
    
    sub_Tab1  <-  Tab1 %>%
      dplyr::filter(site == i) %>%
      dplyr::select(paste0("x",1:num_covariate),A,Y)
    
    formula <- paste0( "Y ~ (", paste0(colnames(sub_Tab1)[1:num_covariate], collapse = " + "),") * A" )
    
    
    
    mod       <-  lm(formula, data=sub_Tab1 )
    
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
#' @param ppsi the number of blip function parameters; ppsi = 11 or 21 in simulations in many covariates setting
#' @param pbeta the number of treatment-free function parameters; pbeta = 11 or 21 in simulations in many covariates setting
#' @param estimate matrix of site-specific blip function parameter estimates
#' @param varest matrix of variances associated with site-specific blip function parameter estimates
#' @param main main treatment effect that is not shrunk
#' @param para common blip function parameters, i.e., treatment-covariate interactions
#' @param sdpara between-site standard deviation
#' @param tau global shrinkage parameter
#' @param lambda local shrinkage parameters
#' @param aux1_global,aux2_global,aux1_local,aux2_local,aux1_sdpara,aux2_sdpara parameters in half-t prior
#' @param nu_local,nu_global,nu_sdpara parameters in half-t prior; set as 1, so it becomes half-cauchy prior
#' @param scale_global,scale_sd scale parameters for global shrinkage parameters and between-site standard deviation
#' @param priormean mean of normal prior
#' @param priorsd standard deviation of normal prior
#' @export
#'


stancode_twostage<-"
  data{
  
  int<lower=1> k; 
  
  int<lower=1> ppsi;

  
  vector[ppsi] estimate[k]; 
  
  vector<lower=0>[ppsi] varest[k]; 
  
  real priormean; 
  real<lower=0> priorsd; 
  
  real<lower=0> scale_global; // scale for tau
  
  real<lower=0> scale_sd;
  
  
  real<lower=1> nu_local;
  real<lower=1> nu_global;
  real<lower=1> nu_sdpara;
  }
  
  
  parameters{
 real main;
  
 real<lower=0> aux1_global; 
 real<lower=0> aux2_global; 
 
 
 vector<lower=0>[ppsi-1] aux1_local; 
 vector<lower=0>[ppsi-1] aux2_local;
 
 vector<lower=0>[ppsi] aux1_sdpara;
 vector<lower=0>[ppsi] aux2_sdpara;
 
  
  vector[ppsi - 1] z;
  }
  
  transformed parameters{
  vector[ppsi-1] para;
  
  real<lower=0> tau; 
  vector<lower=0>[ppsi-1] lambda;  
  
  vector<lower=0>[ppsi] sdpara;
 
  
  lambda = aux1_local .* sqrt(aux2_local); 
  tau = aux1_global * sqrt(aux2_global) * scale_global; 
  sdpara = aux1_sdpara .* sqrt(aux2_sdpara) * scale_sd;
  

  para = z .* lambda * tau;
  
  }
  
  
  model{
  
  for(i in 1:k)
{
estimate[i,1] ~ normal( main, sqrt(varest[i,1] + sdpara[1]^2 ) );
for(j in 2:ppsi)
{
estimate[i,j] ~ normal( para[j-1], sqrt(varest[i,j] + sdpara[j]^2 ) );
}
}
      
   
   main ~ normal( priormean, priorsd );
   z ~ normal(0, 1);
   
   aux1_local ~ normal(0, 1);
   aux2_local ~ inv_gamma(0.5*nu_local, 0.5*nu_local);
   
   aux1_global ~ normal(0, 1);
   aux2_global ~ inv_gamma(0.5*nu_global, 0.5*nu_global);
  
   aux1_sdpara ~ normal(0,1);
   aux2_sdpara ~ inv_gamma(0.5*nu_sdpara, 0.5*nu_sdpara);

  }
"

stan_twostage<-stan_model(model_code = stancode_twostage)



# Fit the two-stage approach
twostage_fun <- function(Tab2,k=10, ppsi=21,pbeta=21)
{
  
  estimate <- matrix(0, nrow=k, ncol = ppsi )
  
  
  for(i in 1:k)
    estimate[i,] <- Tab2[[i]]$estimate[(pbeta+1):(pbeta+ppsi)]
  
  
  varest <- matrix(0, nrow = k, ncol = ppsi )
  
  for(i in 1:k)
    varest[i,] <- diag(Tab2[[i]]$covariance)[(pbeta+1):(pbeta+ppsi)]
  
  
  data_twostage <- list(
    
    ppsi = ppsi,
    
    k = 10,
    
    estimate = estimate,
    varest = varest,
    
    priormean = 0,
    priorsd = 100,
    
    scale_global = 1,
    
    nu_local = 1,
    nu_global = 1,
    nu_sdpara=1,
    
    scale_sd = 1
  )
  
  
  twostage <- sampling(stan_twostage,data = data_twostage, chains = 2, iter =2000,
                         control = list(adapt_delta = 0.99,
                                        max_treedepth = 20 ))
  
  return(twostage)
}


# -----------------------------------------------------------------
#  stan code/model for one-stage approach
# -----------------------------------------------------------------

#'
#' @param k the number of sites; k=10 in simulations in a sparse data setting
#' @param N total sample size
#' @param ppsi the number of blip function parameters; ppsi = 11 or 21 in simulations in many covariates setting
#' @param pbeta the number of treatment-free function parameters; pbeta = 11 or 21 in simulations in many covariates setting
#' @param Xmatrix design matrix
#' @param reward patient outcomes
#' @param site vector of site ID for each patient
#' @param nu_local,nu_global,nu_betweensd,nu_withinsd parameters in half-cauchy prior
#' @param commonbeta common treatment-free parameters
#' @param main_psi main treatment effect
#' @param common common blip parameters (not including main treatment effect)
#' @param aux1_global,aux2_global,aux1_local,aux2_local,aux1_betweensd,aux2_betweensd,aux1_withinsd,aux2_withinsd parameters in half-t prior
#' @param nu_local,nu_global,nu_betweensd,nu_withinsd parameters in half-t prior; set as 1, so it becomes half-cauchy prior
#' @param scale_global,scale_sd,scale_intercept scale parameters for global shrinkage parameters, between-site standard deviation, and intercept
#' @param tau global shrinkage parameter
#' @param lambda local shrinkage parameters
#' @param betweensd between-site standard deviation
#' @param withinsd within-site standard deviation
#' @param sitepar site-specific parameters
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
  
  real<lower=0> scale_intercept; 
  real<lower=0> scale_global; 
  
  real<lower=0> scale_sd;
  
  
  real<lower=1> nu_local;
  real<lower=1> nu_global;
  real<lower=1> nu_betweensd;
  real<lower=1> nu_withinsd;
  }
  
  
  parameters{
 
  vector[pbeta] commonbeta;
  real main_psi;
  
  matrix[k,pbeta+ppsi] eta;
  
 real<lower=0> aux1_global; 
 real<lower=0> aux2_global; 
 
 
 vector<lower=0>[ppsi-1] aux1_local; 
 vector<lower=0>[ppsi-1] aux2_local;
 
 vector<lower=0>[pbeta+ppsi] aux1_betweensd;
 vector<lower=0>[pbeta+ppsi] aux2_betweensd;
 
 vector<lower=0>[k] aux1_withinsd;
 vector<lower=0>[k] aux2_withinsd;
  
  vector[ppsi - 1] z;
  }
  
  transformed parameters {

  matrix[k,pbeta+ppsi] sitepar;
  vector[ppsi -1] common;
  
  real<lower=0> tau; 
  vector<lower=0>[ppsi-1] lambda;  
  
  vector<lower=0>[pbeta+ppsi] betweensd;
  vector<lower=0>[k] withinsd;
  
  lambda = aux1_local .* sqrt(aux2_local); 
  tau = aux1_global * sqrt(aux2_global) * scale_global; 
  withinsd = aux1_withinsd .* sqrt(aux2_withinsd) * scale_sd;
  betweensd = aux1_betweensd .* sqrt(aux2_betweensd) * scale_sd;
  

  common = z .* lambda * tau;
  
  for(s in 1:k)
  { 
  for(p in 1:pbeta)
   {
   sitepar[s,p] = commonbeta[p] + betweensd[p] * eta[s,p];
   }
  }
  
  for(s in 1:k)
{
  sitepar[s, pbeta+1] = main_psi + betweensd[pbeta+1] * eta[s, pbeta+1];
}
  
  for(s in 1:k)
  for(p in (pbeta+2):(pbeta+ppsi))
  sitepar[s,p] = common[p-1-pbeta] + betweensd[p] * eta[s,p];
  
  }
  
  
  model{
  
  for(i in 1:N)
     reward[i] ~ normal(Xmatrix[i,] * (sitepar[site[i],])', 
     withinsd[site[i]]);
  
   for(s in 1:k)
   for(p in 1:(pbeta+ppsi))
   eta[s,p] ~ normal(0,1);
        
   commonbeta ~ normal( 0, scale_intercept );
   
   z ~ normal(0, 1);
   
   aux1_local ~ normal(0, 1);
   aux2_local ~ inv_gamma(0.5*nu_local, 0.5*nu_local);
   
   aux1_global ~ normal(0, 1);
   aux2_global ~ inv_gamma(0.5*nu_global, 0.5*nu_global);
  
   aux1_betweensd ~ normal(0,1);
   aux2_betweensd ~ inv_gamma(0.5*nu_betweensd, 0.5*nu_betweensd);
   
   aux1_withinsd ~ normal(0,1);
   aux2_withinsd ~ inv_gamma(0.5*nu_withinsd,0.5*nu_withinsd);
   
  }
"



stan_onestage<-stan_model(model_code = stancode_onestage)





onestage_fun <- function(Tab1,k=10,N,pbeta=21, ppsi = 21)
{
  
  formula_char <- paste0( "~(", paste0( paste0("x",1:(ppsi-1)), collapse = " + "),") * A" )
  
  data_onestage <- list(
    k = k,
    N = N,
    
    pbeta = pbeta,
    ppsi = ppsi,
    
    Xmatrix = model.matrix(as.formula(formula_char), data = Tab1),
    
    reward = Tab1$Y,
    
    site = Tab1$site,
    
    scale_intercept = 10,
    scale_global = 1,
    
    nu_local = 1,
    nu_global = 1,
    nu_withinsd=1,
    nu_betweensd=1,
    
    
    scale_sd = 1
    
  )
  
  
  
  onestage <- sampling(stan_onestage, data = data_onestage, chains=2, iter = 20,
                  control = list(adapt_delta = 0.99,
                                 max_treedepth = 2000))
  return(onestage)
}



###########################################################
#           Example for one simulated dataset
###########################################################

k <- 10
N <- 2000

sample_size_site <- c( 200*0.06 , 200*0.07 , 200*0.08 , 200*0.09 , 200*0.10 , 200*0.10 , 200*0.11 , 200*0.12 , 200*0.13 , 200*0.14 )*10 # if large sample

psi_df <- data.frame(psi_0 = 2.5, psi_1 = -0.5 , psi_2 = 2, psi_3 = -1, psi_4 = 0,
                     psi_5 = 0, psi_6 = 0, psi_7 = 0, psi_8 = 0,
                     psi_9 = 0, psi_10 = 0, 
                     psi_11 = 0 , psi_12 = 0, psi_13 = 0, psi_14 = 0,
                     psi_15 = 0, psi_16 = 0, psi_17 = 0, psi_18 = 0,
                     psi_19 = 0, psi_20 = 0
)

coeff_outcome_df <- data.frame(   beta_0 = 4      ,
                                  beta_1 = 1    ,
                                  beta_2 = 1,
                                  beta_3 = 1,
                                  beta_4 = 1,
                                  beta_5 = 1,
                                  beta_6 = 1,
                                  beta_7 = 1,
                                  beta_8 = 1,
                                  beta_9 = 1,
                                  beta_10 = 1,
                                  beta_11 = 1    ,
                                  beta_12 = 1,
                                  beta_13 = 1,
                                  beta_14 = 1,
                                  beta_15 = 1,
                                  beta_16 = 1,
                                  beta_17 = 1,
                                  beta_18 = 1,
                                  beta_19 = 1,
                                  beta_20 = 1
)




heter_level <- 0.1

sigma_w <- 0.25

Tab1 <-Individual_Data_sim_f( k = k, 
                              num_covariate = 20,
                              sample_size_site = sample_size_site ,
                              heter_level = heter_level,
                              psi_df = psi_df,
                              coeff_outcome_df = coeff_outcome_df,
                              sigma_w = sigma_w )


mod_onestage <- onestage_fun(Tab1 = Tab1,
                             k = 10,
                             N = 2000,
                             pbeta = 21,
                             ppsi = 21)

Tab2 <- Individual_lm_Data_sim_f( Tab1 = Tab1, num_covariate = 20 )

mod_twostage <- twostage_fun(Tab2 = Tab2,
                             k = 10,
                             ppsi = 21,
                             pbeta = 21)





