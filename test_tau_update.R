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
  log_Cd_vals <- (d/2 - 1) * log(tau) - log(besselI(tau, nu = (d/2 - 1), expon.scaled = T))
  
  # Compute dot products (phi_k^T E[y_i | x_i, z_i = k])
  dot_products <- Ey %*% phi_k
  
  # Compute Q(tau)
  Q <- sum(W_k * (log_Cd_vals + tau * dot_products))
  return(Q)
}




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

### 
############ EM iterations
iter <- 1
  message("iteration ", iter)
  # E-step
  e_step_rst <- E_step(x, pi_hat, c_hat, r_hat, phi_hat, tau_hat, sigma_sq)
  W <- e_step_rst$W
  y_expect_list <- e_step_rst$y_expect_list
  
 k <-1
    Ey_k <- y_expect_list[[k]]
    W_k <- W[, k]
    tau_values <- seq(0.1, 100, length.out = 100)
    Q_values <- sapply(tau_values, Q_tau, phi_k = phi_hat[k,], Ey = Ey_k, W_k = W_k)
    plot(tau_values, Q_values, type = "l", xlab = "tau_k", ylab = "Q(tau_k)")
  
  
    log_Cd_vals<- function(tau, d){
      temp <- (d/2 - 1) * log(tau) - log(besselI(tau, nu = (d/2 - 1), expon.scaled = TRUE))
      return(temp)
    }
    
    I_vals <- function(tau, d){
      temp <- (- log(besselI(tau, nu = (d/2 - 1), expon.scaled = TRUE)))
      return(temp)
    }
    
    log_Cd_val <- sapply(tau_values, log_Cd_vals, d =2)  # l
    I <- sapply(tau_values, I_vals, d =2)
    plot(tau_values, I, type = "l")
    
    a <- b <- numeric(length(tau_values))
    for (i in 1:length(tau_values)){
      
      a[i] <- (- log(besselI(tau_values[i], nu = (d/2 - 1))))
      b[i] <- (- log(besselI(tau_values[i], nu = (d/2 - 1), expon.scaled = T)))
      
    }
    
  # M-step
  m_step_rst <- M_step(x, W, y_expect_list, c_hat, r_hat, phi_hat, tau_hat)
  pi_hat <- m_step_rst$pi_hat
  c_hat <- m_step_rst$c_hat
  r_hat <- m_step_rst$r_hat
  phi_hat <- m_step_rst$phi_hat
  tau_hat <- m_step_rst$tau_hat
  sigma_sq <- m_step_rst$sigma_sq
  