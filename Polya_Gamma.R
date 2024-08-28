polya_gamma <- function(y, X,
                        mu_beta, Sigma_beta,
                        n_mcmc){
  
  ###
  ### Packages
  ###
  
  # library(BayesLogit)
  # library(Boom)
  
  
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
  
  for(q in 1:n_mcmc){
    
    ### Sample omega
    omega=rpg(nrow(X), 1, X%*%beta)
    
    ### Sample beta
    kappa=y-1/2
    V_omega=solve(t(X)%*%diag(omega)%*%X+Sigma_beta_inv)
    m_omega=V_omega%*%(t(X)%*%kappa+Sigma_beta_inv_times_mu)
    beta=t(mvnfast::rmvn(1, m_omega, V_omega))
    
    ### Save Samples
    beta_save[,q]=beta

    ###
    ### Timer
    ###
    
    if (q %% 100 == 0) {cat(q, " ")}
    
  }
  
  list(beta=beta_save)
}