polya_gamma <- function(y, X,
                        mu_beta, Sigma_beta,
                        n_mcmc){
  
  ###
  ### Packages
  ###
  
  require(pgdraw)
  
  ###
  ### Loop Variables
  ### 
  
  P <- solve(Sigma_beta)
  Pmu  <- c(P%*%mu_beta)
  
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
  
  kappa=(y-0.5)
  # n <- nrow(X)
  
  for(q in 1:n_mcmc){
    
    ### Sample omega
    
    omega <- pgdraw(1, tcrossprod(X, t(beta)))
    
    ### Sample beta
    # omega_X <- X * omega
    # V_omega=solve(crossprod(X, omega_X))
    # m_omega=V_omega%*%(crossprod(X, kappa) + Sigma_beta_inv_times_mu)
    
    P_vb <- crossprod(X*omega,X) + P
    Sigma <- solve(P_vb) 
    mu <- Sigma %*% (crossprod(X,kappa) + Pmu)
    beta=t(mvnfast::rmvn(1, mu, Sigma))
    
    ### Save Samples
    beta_save[,q]=beta
    
    ###
    ### Timer
    ###
    
    if (q %% 1000 == 0) {cat(q, " ")}
    
  }
  
  list(beta=beta_save)
}