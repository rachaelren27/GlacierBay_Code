PG_VB <- function(y, X, mu_beta, Sigma_beta, n_mcmc){
  
  # starting values
  xi <- X%*%mu_beta
  
  # storage
  m_save <- matrix(NA, ncol(X), n_mcmc)
  V_save <- list()
  
  # loop variables
  Sigma_beta_inv=solve(Sigma_beta)
  Sigma_beta_inv_times_mu=Sigma_beta_inv%*%mu_beta
  
  kappa= y-1/2
  
  calc_xi <- function(x, m, V){
    sqrt(t(X)%*%V%*%x + (t(x)%*%m)^2)
  }
  
  # CAVI
  for(t in 1:10){
    z <- (1/2)*(xi)^(-1)*tanh((1/2)*xi)
    
    z_X <- X * z
    V <- solve(crossprod(X, z_X) + Sigma_beta_inv)
    m <- V %*% (crossprod(X, kappa) + Sigma_beta_inv_times_mu)
    
    m_save[,t] <- m
    V_save[[t]] <- V
    
    xi <- apply(X, 1, calc_xi, m = m, V = v)
  }
}