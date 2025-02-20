#### debug Q_tau

Q_tau <- function(tau, phi_k, Ey, W_k) {
  d <- length(phi_k)
  log_Cd_val <- log_Cd(tau, d)  # log(C_d(tau))
  dot_products <- Ey %*% phi_k  # phi_k^T E[y_i | x_i, z_i = k]
  Q <- sum(W_k * (log_Cd_val + tau * dot_products))
  return(Q)
}


Q_tau <- function(tau, phi_k, Ey, W_k) {
  d <- length(phi_k)
  
  # Compute log(C_d(tau))
  log_Cd_vals <- (d/2 - 1) * log(tau) - log(besselI(tau, nu = (d/2 - 1)))
  
  # Compute dot products (phi_k^T E[y_i | x_i, z_i = k])
  dot_products <- Ey %*% phi_k
  
  # Compute Q(tau)
  Q <- sum(W_k * (log_Cd_vals + tau * dot_products))
  return(-Q)
}


tau_init <- tau_list

optim(
  par = tau_init[k], 
  fn = Q_tau, 
  phi_k = phi_hat[k,], 
  Ey = Ey_k, 
  W_k = W_k, 
  method = "L-BFGS-B", 
  lower = 1e-2  # Ensuring tau is positive
)






n <- nrow(x)
d <- ncol(x)

################### initialization
initialPars <- initial_kmeans(x, M)
## double check the original initialization
#initialPars <- initial_kmeans_opt1(x, M)


pi_hat <- initialPars$pi_hat
c_hat <- initialPars$c_hat
r_hat <- initialPars$r_hat
sigma_sq <-initialPars$sigma_sq
phi_hat <- initialPars$phi_hat
tau_hat <- initialPars$tau_hat

### use the true value of c_hat, r_hat, sigma_sq and phi_hat
c_hat <- matrix(unlist(c_list),ncol=2,byrow = T)
r_hat <- r_list
phi_hat <- matrix(unlist(phi_list),ncol=2,byrow = T)
sigma_sq <- 0.5^2

############ EM iterations
iter <- 1
  message("iteration ", iter)
  # E-step
  e_step_rst <- E_step(x, pi_hat, c_hat, r_hat, phi_hat, tau_hat, sigma_sq)
  W <- e_step_rst$W
  y_expect_list <- e_step_rst$y_expect_list
  
for (k in 1:M){
   Ey_k <- y_expect_list[[k]]
   W_k <- W[, k]
   temp <- optim(par = tau_init[k], 
                       fn = Q_tau, 
                       phi_k = phi_hat[k,], 
                       Ey = Ey_k, 
                       W_k = W_k, 
                       method = "L-BFGS-B", 
                       lower = 1e-2 )
   tau_hat[k] <- temp$par
}
   
    
  
    tau_values <- seq(0.1, 100, length.out = 100)
    Q_values <- sapply(tau_values, Q_tau, phi_k = phi_hat[k,], Ey = Ey_k, W_k = W_k)
    plot(tau_values, Q_values, type = "l", xlab = "tau_k", ylab = "Q(tau_k)")
  
  
    
    
    a <- b <- numeric(length(tau_values))
    for (i in 1:length(tau_values)){
      
      a[i] <- (- log(besselI(tau_values[i], nu = (d/2 - 1))))
      b[i] <- (- log(besselI(tau_values[i], nu = (d/2 - 1), expon.scaled = T)))
      
    }
    
