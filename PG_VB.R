# PG_VB <- function(y, X, mu_beta, Sigma_beta, n_iter){
#   
#   # starting values
#   xi <- X%*%mu_beta
#   
#   # storage
#   m_save <- matrix(NA, ncol(X), n_iter)
#   V_save <- list()
#   
#   # loop variables
#   Sigma_beta_inv=solve(Sigma_beta)
#   Sigma_beta_inv_times_mu=Sigma_beta_inv%*%mu_beta
#   
#   kappa= y-1/2
#   
#   calc_xi <- function(x, m, V){
#     sqrt(t(x)%*%V%*%x + (t(x)%*%m)^2)
#   }
#   
#   # CAVI
#   for(t in 1:n_iter){
#     z <- c((1/2)*(xi)^(-1)*tanh((1/2)*xi))
#     
#     z_X <- X * z
#     V <- solve(crossprod(X, z_X) + Sigma_beta_inv)
#     m <- V %*% (crossprod(X, kappa) + Sigma_beta_inv_times_mu)
#     
#     m_save[,t] <- m
#     V_save[[t]] <- V
#     
#     xi <- apply(X, 1, calc_xi, m = m, V = V)
#     
#     if (t %% 100 == 0) {cat(t, " ")}
#   }
#   
#   return(list(m_save = m_save, V_save = V_save))
# }

# Durante's function

logit_CAVI <- function(X, y, prior, tol = 1e-16, maxiter=10000){
  
  if (is.null(n <- nrow(X)))
    stop("'X' must be a matrix")
  
  # Compute the log-determinant of a matrix
  ldet <- function(X) {
    if(!is.matrix(X)) return(log(X))
    determinant(X,logarithm = TRUE)$modulus
  }
  
  lowerbound <- numeric(maxiter)
  p          <- ncol(X)
  
  P    <- solve(prior$Sigma)
  mu   <- prior$mu
  Pmu  <- c(P%*%mu)
  Pdet <- ldet(P)
  
  # Initialization for omega equal to 0.25
  P_vb       <- crossprod(X*rep(1/4,n),X) + P
  Sigma_vb   <- solve(P_vb) 
  mu_vb      <- Sigma_vb %*% (crossprod(X,y - 0.5) + Pmu)
  eta        <- c(X%*%mu_vb)
  xi         <- sqrt(eta^2 +  rowSums(X %*% Sigma_vb * X))
  omega      <- tanh(xi/2)/(2*xi); 
  omega[is.nan(omega)] <- 0.25
  
  lowerbound[1]  <- 0.5*p + 0.5*ldet(Sigma_vb) + 0.5*Pdet - 0.5*t(mu_vb - mu)%*%P%*%(mu_vb - mu) + sum((y-0.5)*eta +log(plogis(xi)) - 0.5*xi) - 0.5*sum(diag(P %*% Sigma_vb))
  
  # Iterative procedure
  for(t in 2:maxiter){
    
    P_vb       <- crossprod(X*omega,X) + P
    Sigma_vb   <- solve(P_vb) 
    mu_vb      <- Sigma_vb %*% (crossprod(X,y-0.5) + Pmu)
    
    #Update of xi
    eta        <- c(X%*%mu_vb)
    xi         <- sqrt(eta^2 +  rowSums(X %*% Sigma_vb * X))
    omega      <- tanh(xi/2)/(2*xi); 
    omega[is.nan(omega)] <- 0.25
    
    lowerbound[t]  <- 0.5*p + 0.5*ldet(Sigma_vb) + 0.5*Pdet - 0.5*t(mu_vb - mu)%*%P%*%(mu_vb - mu) + sum((y-0.5)*eta +log(plogis(xi)) - 0.5*xi) - 0.5*sum(diag(P %*% Sigma_vb))
    
    if(abs(lowerbound[t] - lowerbound[t-1]) < tol) return(list(mu = matrix(mu_vb,p,1), Sigma=matrix(Sigma_vb,p,p), 
                                                               Convergence=cbind(Iteration=(1:t)-1, 
                                                                                 Lowerbound=lowerbound[1:t]), xi=xi))
  }
  stop("The algorithm has not reached convergence")
}
