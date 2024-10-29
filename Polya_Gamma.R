polya_gamma <- function(y, X, w,
                        mu_beta, Sigma_beta,
                        n_mcmc){
  
  ###
  ### Packages
  ###
  
  # library(BayesLogit)
  # library(Boom)
  require(pgdraw)
  
  ###
  ### Loop Variables
  ### 
  
  Sigma_beta_inv=solve(Sigma_beta)
  Sigma_beta_inv_times_mu=Sigma_beta_inv%*%mu_beta
  
  ###
  ### Starting Values 
  ###
  
  ### Cool Start
  beta=mu_beta
  
  ###
  ### Storage
  ###
  
  beta_save=matrix(NA, ncol(X), n_mcmc)
  
  ###
  ### MCMC loop
  ###
  
  kappa=w*(y-1/2)
  # n <- nrow(X)
  
  for(q in 1:n_mcmc){
    
    ### Sample omega
    omega <- pgdraw(w, X%*%beta)
    
    ### Sample beta
    omega_X <- X * omega
    V_omega=solve(crossprod(X, omega_X))
    m_omega=V_omega%*%(crossprod(X, kappa) + Sigma_beta_inv_times_mu)
    beta=t(mvnfast::rmvn(1, m_omega, V_omega))
    
    ### Save Samples
    beta_save[,q]=beta
    
    ###
    ### Timer
    ###
    
    if (q %% 1000 == 0) {cat(q, " ")}
    
  }
  
  list(beta=beta_save)
}