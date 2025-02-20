library(movMF)
library(mvtnorm)
library(matrixStats)



# L2 norm function
norm_vec <- function(v) {
  norm_val <- sqrt(sum(v^2))
  if (norm_val < 1e-10) return(1e-10)  # Avoid division by zero
  return(norm_val)
}



# Function to calculate log(Cd(tau))
# Cd(tau) = tau^{d/2-1}(2pi)^{-d/2}(I_{d/2-1}(tau))^{-1}
# besselI() function in base R package computes the modified Bessel function of the first kind

log_Cd <- function(tau, d) {
  # Check for invalid inputs
  if (tau < 0) stop("tau must be non-negative")
  
  # Handle the case where tau is zero
  if (tau == 0) {
    # For tau = 0, the vMF distribution is uniform on the sphere
    # The normalizing constant is C_d(0) = 1 / SurfaceArea(S^{d-1})
    # The surface area of S^{d-1} is 2 * pi^(d/2) / gamma(d/2)
    # So log(C_d(0)) = -log(2 * pi^(d/2) / gamma(d/2))
    # lgamma computes the logarithm of the gamma function
    return(-log(2) - (d/2) * log(pi) + lgamma(d/2))
  }
  
  # Compute log_Cd(tau) for tau > 0
  log_tau <- ifelse(tau > 0, log(tau), 0)
  bessel_term <- besselI(tau, nu = (d/2 - 1), expon.scaled = FALSE)
  
  # Handle cases where bessel_term is zero or very small
  if (bessel_term <= 0) {
    return(-Inf)  # log(0) = -Inf
  }
  
  # Compute log_Cd(tau)
  res <- (d/2 - 1) * log_tau - log(bessel_term) - (d/2) * log(2*pi)
  return(res)
}



# Function to calculates FG-kernel density given the parameters cntr(center), rd(radius), sigma.sq, phi(direction), tau(concentration)
# this mainly comes for the FG paper code, but instead of using ((norm_vec(z))^2+rd^2-2*sigma.sq*dinom), I used ((norm_vec(z))^2 + rd^2)

FG_kernel <- function(x_i, c_k, r_k, sigma_sq, phi_k, tau_k) { 
  d <- length(x_i)
  
  # Check for invalid inputs
  if (sigma_sq <= 0) stop("sigma_sq must be positive")
  if (tau_k <= 0) stop("tau_k must be positive")
  
  # Compute z = x_i - c_k
  z <- x_i - c_k
  
  # Compute temp = ||tau_k * phi_k + (r_k * z) / sigma_sq||
  temp <- norm_vec(tau_k * phi_k + (r_k * z) / sigma_sq)
  
  # Compute log_Cd(tau_k, d) and log_Cd(temp, d)
  log_Cd_tau <- log_Cd(tau_k, d)
  log_Cd_temp <- log_Cd(temp, d)
  
  # Compute the log-density
  ln_h_x_theta <- log_Cd_tau - (d * log(2 * pi * sigma_sq) / 2) - 
    log_Cd_temp - ((norm_vec(z))^2 + (r_k)^2) / (2 * sigma_sq)
  
  # Handle numerical underflow by returning exp(ln_h_x_theta)
  return(exp(ln_h_x_theta))
}

 

# Function to compute E(y_i|x_i,z_i, Theta)
compute_y_expect <- function(x, M, phi, tau, c, r, sigma_sq){
  n <- nrow(x)
  d <- ncol(x)
  
  E_y <- list()
  
  for(k in 1: M){
    y_expect_k <- matrix(0, n, d)
    
    for (i in 1:n){
    nu <- tau[k] * phi[k,] + (r[k] * (x[i,] - c[k,]))/sigma_sq
    kappa_k <- norm_vec(nu)
    mu_k <- nu/kappa_k
    Ad_k <- besselI(kappa_k, d/2) / besselI(kappa_k, d/2 - 1)
    y_expect_k[i,] <- Ad_k*mu_k	
    }
    E_y[[k]] <- y_expect_k
  }
  return(E_y)
}


## this version is more efficient
compute_y_expect <- function(x, M, phi, tau, c, r, sigma_sq) {
  n <- nrow(x)
  d <- ncol(x)
  
  E_y <- list()
  
  for (k in 1:M) {
    # Compute nu for all points in the cluster at once
    nu <- t(tau[k] * phi[k, ] + t(r[k] * (x - matrix(c[k, ], n, d, byrow = TRUE)) / sigma_sq))
    
    # Compute kappa_k and mu_k
    kappa_k <- apply(nu, 1, norm_vec)
    # Handle cases where kappa_k is zero
    kappa_k[kappa_k == 0] <- 1e-10  # Small value to avoid division by zero
    mu_k <- nu / kappa_k
    
   
    # Compute A_d(kappa_k) using exponentially scaled Bessel functions for numerical stability
    log_Ad_k <- log(besselI(kappa_k, d/2, expon.scaled = FALSE)) - 
      log(besselI(kappa_k, d/2 - 1, expon.scaled = FALSE))
    Ad_k <- exp(log_Ad_k)
    Ad_k <- pmin(pmax(Ad_k, 1e-10), 1 - 1e-10)
    
    # Compute E[y_i | x_i, z_i = k]
    y_expect_k <- Ad_k * mu_k
    
    E_y[[k]] <- y_expect_k
  }
  
  return(E_y)
}



newton_tau <- function(tau_init, phi_hat, y_expect_list, W){
  for (k in 1: M){
    NUMERATOR <- numeric(n)
    DENUMERATOR <- numeric(n)
    phi_k <- phi_hat[k, ]
    d <- length(phi_k)
    y_expect_k <- y_expect_list[[k]]
    
    for (i in 1:n){
      Ad_tau_k <- besselI(tau_init[k], nu = d/2, expon.scaled = FALSE)/besselI(tau_init[k], d/2 - 1, expon.scaled = FALSE)
      Ad_tau_k_prime <- 1 - Ad_tau_k^2 - (d-1)/tau_init[k]*Ad_tau_k
      dot_product <- sum(phi_k * y_expect_k[i,])
      NUMERATOR[i] <- W[i,k]*(Ad_tau_k - dot_product)
      DENUMERATOR[i] <- W[i,k]*Ad_tau_k_prime
    }
    gradient_tau_k <- sum(NUMERATOR)/sum(DENUMERATOR)
    tau_hat[k] <- tau_init[k] - gradient_tau_k
  }
  return(tau_hat)
  
}


### another version
estimate_tau_k <- function(R_bar, p, max_iter = 100, tol = 1e-6) {
  if (R_bar < 1e-10) {
    return(0)  # Handle the case where R_bar is close to 0
  }
  
  # Initial guess for tau_k
  tau_k <- 1  # Start with a small positive value
  
  for (iter in 1:max_iter) {
    # Compute A_p(tau_k) and its derivative
    A_p_val <- A_p(tau_k, p)
    A_p_deriv <- A_p_prime(tau_k, p)
    
    # Newton update
    delta <- (A_p_val - R_bar) / A_p_deriv
    tau_k_new <- tau_k - delta
    
    # Ensure tau_k remains non-negative
    tau_k_new <- max(0, tau_k_new)
    
    # Check for convergence
    if (abs(tau_k_new - tau_k) < tol) {
      break
    }
    
    tau_k <- tau_k_new
  }
  
  return(tau_k)
}





## update tau based on yi, vMF
update_tau <- function(x, W, c_hat, r_hat){
  n <- nrow(x)
  M <- ncol(W)
  
  # recover y_i and update tau_k
  y <- matrix(0, nrow = n)
  for(i in 1:n){
    cluster_idx <- which.max(W[i,])
    y_i <- (x[i,] - c_hat[cluster_idx,  ])/ r_hat[cluster_idx]
    y[i, ] <- y_i/sqrt(sum(y_i^2)) # normalzied y_i
  }
  
  tau_hat <- numeric(M)
  
  for (k in 1:M){
    weighted_sum <- colSums(W[, k]*y)
    R_bar <- sqrt(sum(weighted_sum^2)) / sum(gamma_ik[, k])
    tau_hat[k] <- estimate_tau_k(R_bar, p)  # Use the function defined earlier
  }
  
  return(tau_hat)
  
}


### alternative way to update tau

Q_tau <- function(tau, phi_k, Ey, W_k) {
  d <- length(phi_k)
  log_Cd_val <- log_Cd(tau, d)  # log(C_d(tau))
  dot_products <- Ey %*% phi_k  # phi_k^T E[y_i | x_i, z_i = k]
  Q <- sum(W_k * (log_Cd_val + tau * dot_products))
  return(Q)
}

optimize_tau <- function(phi_k, Ey, W_k) {
  # Define the negative Q function (since optimize minimizes by default)
  neg_Q_tau <- function(tau) {
    - Q_tau(tau, phi_k, Ey, W_k)
  }
  
  # Find tau_k that maximizes Q(tau)
  result <- optimize(neg_Q_tau, interval = c(1e-6, 100))  # Adjust interval as needed
  return(result$minimum)
}





# Function to calculate the incomplete log-likelihood

log_likelihood <- function(x, pi_hat, c_hat, r_hat, phi_hat, tau_hat, sigma_sq) {
  n <- nrow(x)
  M <- length(pi_hat)
  log_lik <- matrix(0, n, M)
  temp <- numeric(n)
  for (i in 1:n) {	
    for (k in 1:M) {
      # Compute p(x_i | z_i = k, theta), which is FG(theta)
      # we can calculate using the FG kernel
      log_lik[i,k] <- pi_hat[k] * FG_kernel(x[i,], c_hat[k,], r_hat[k],sigma_sq, phi_hat[k, ], tau_hat[k])
    }
	temp[i] <-log(sum(log_lik[i,]))
  }
  rst <- sum(temp)
  
  return(rst)
}

## this version produce the same results as the above function
### logSumExp calculate the log of the sum of exponentials
## eg., originally, we calculate pi_1*FG1, pi_2*FG2, etc, then we calculate log(pi_1*FG1 + pi_2*FG2)
## now we calculate t1= log(pi_1*FG1), t2= log(pi_2*FG2), etc, then logSumExp(t1, t2), which is log(exp(t1)+ exp(t2))
## which is just log(pi_1*FG1 + pi_2*FG2)


log_likelihood <- function(x, pi_hat, c_hat, r_hat, phi_hat, tau_hat, sigma_sq) {
  n <- nrow(x)
  M <- length(pi_hat)
  log_lik <- matrix(0, n, M)
  temp <- numeric(n)
  
  for (i in 1:n) {
    for (k in 1:M) {
      fg_value <- FG_kernel(x[i, ], c_hat[k, ], r_hat[k], sigma_sq, phi_hat[k, ], tau_hat[k])
      if(fg_value <=0){
        cat(sprintf("⚠️ Warning: FG_kernel returned non-positive value at i=%d, k=%d: %f\n", i, k, fg_value))
      }
      
      log_lik[i, k] <- log(pi_hat[k]) + log(FG_kernel(x[i, ], c_hat[k, ], r_hat[k], sigma_sq, phi_hat[k, ], tau_hat[k])+1e-10)
    }
    temp[i] <- logSumExp(log_lik[i, ])  # Use logSumExp for numerical stability
    if (i %% 10 == 0) {  # Print every 10 iterations to avoid too much output
      cat(sprintf("🔹 Iteration %d: log-lik row sum = %f\n", i, temp[i]))
    }
  }
  
  rst <- sum(temp)
  return(rst)
}



## function to compute the posterior probability w_{ik}= P(z_i=k|x_i, Theta)
compute_W <- function(x, pi_hat, c_hat, r_hat, phi_hat, tau_hat, sigma_sq){
  n <- nrow(x)  # Number of data points
  M <- length(pi_hat)  # Number of clusters
  
  W <- matrix(0, n, M)  
  for(i in 1:n){
    for (k in 1: M){
      # compute p(x_i|z_i=k, Theta)
      p_x_given_z <- FG_kernel(x[i,], c_hat[k, ], r_hat[k], sigma_sq, phi_hat[k, ], tau_hat[k])
      W[i, k] <- pi_hat[k] * p_x_given_z
    }
  }
  return(W)
}

### note to calculate (x_i-c_k)'E(y_i|x_i,z_i=k), 
## if denote a=x_i-c_k, b=E(y_i|x_i,z_i=k), a'b, we simply use sum(a*b)



E_step <- function(x, pi_hat, c_hat, r_hat, phi_hat, tau_hat, sigma_sq){
  n <- nrow(x)
  M <- length(pi_hat)
  W <- matrix(0, n, M)
  y_expect_list <- vector('list', M)
  
  W <- compute_W(x, pi_hat, c_hat, r_hat, phi_hat, tau_hat, sigma_sq)
  W <- W/rowSums(W) # normalize 
  
  
  # Compute E(y_i|x_i, z_i=k, Theta)
  y_expect_list <- compute_y_expect(x, M, phi = phi_hat, tau = tau_hat, c= c_hat, r= r_hat, sigma_sq)
    
  return(list(W = W, y_expect_list = y_expect_list))
}


M_step <- function(x, W, y_expect_list, c_hat, r_hat, phi_hat, tau_hat){
  n <- nrow(x)
  d <- ncol(x)
  M <- ncol(W)
  
  # update mixture weight pi_k
  pi_hat <- colSums(W)/n
  
  # update c_k, r_k, phi_k and tau_k for each cluster
  
  
for (k in 1: M){
    y_expect_k <- y_expect_list[[k]]
	
    # update c_k --> checked, correct
    c_hat[k, ] <- colSums(W[, k] * (x - r_hat[k] * y_expect_k)) / sum(W[, k])
    
    # update r_k --- checked, correct
    temp <- numeric(n)
    for(i in 1:n){
      temp[i] <- W[i,k]*sum((x[i,] - c_hat[k,])*y_expect_k[i,])
    }
    r_hat[k] <- sum(temp)/sum(W[,k])
    
    
    # update phi_k --- checked, correct
    phi_hat[k, ] <- colSums(W[, k] * y_expect_k)
    phi_hat[k, ] <- phi_hat[k, ] / norm_vec(phi_hat[k, ])
	}
  
  ##
  #tau_update <- newton_tau(tau_hat, phi_hat, y_expect_list, W)
  #tau_hat <- tau_update
  
	## alternative way to update tau
  for (k in 1:M){
    Ey_k <- y_expect_list[[k]]
    W_k <- W[, k]
    tau_hat[k] <- optimize_tau(phi_hat[k, ], Ey_k, W_k)
  }
  
  
  ### OPT3:
  
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
  
  
  
  # update sigma_sq --- checked, correct
    temp <- matrix(0, n, M)
    for (i in 1:n) {
      for (k in 1:M) {
        y_expect_k <- y_expect_list[[k]]  # Get E[y_i | x_i, z_i = k] for cluster k
        a <- x[i,]-c_hat[k,]
        norm_a <- sum(a^2)
        temp[i, k] <- W[i,k] * (norm_a - 2 * r_hat[k] * sum(a * y_expect_k[i,]) + (r_hat[k])^2)
      }
    }
    sigma_sq <- sum(temp) / (n * d)
 
  
  return(list(pi_hat = pi_hat, c_hat = c_hat, r_hat = r_hat,
              phi_hat = phi_hat, tau_hat = tau_hat, sigma_sq = sigma_sq))
}

####################### The main EM algorithm

FG_EM <- function(x, M, max_iter, tol ){
  n <- nrow(x)
  d <- ncol(x)
  
  ################### initialization
  initialPars <- initial_kmeans(x, M)
  ## double check the original initialization
  #initialPars <- initial_kmeans_opt1(x, M)
  
  log_likelihood_vals <- numeric(max_iter)
  pi_hat <- initialPars$pi_hat
  c_hat <- initialPars$c_hat
  r_hat <- initialPars$r_hat
  sigma_sq <-initialPars$sigma_sq
  phi_hat <- initialPars$phi_hat
  tau_hat <- initialPars$tau_hat
  ############ EM iterations
  for (iter in 1: max_iter){
    message("iteration ", iter)
    # E-step
    e_step_rst <- E_step(x, pi_hat, c_hat, r_hat, phi_hat, tau_hat, sigma_sq)
    W <- e_step_rst$W
    y_expect_list <- e_step_rst$y_expect_list
    
    ## log-likelihood
    log_likelihood_vals[iter] <- log_likelihood(x, pi_hat, c_hat, r_hat, phi_hat, tau_hat, sigma_sq)
    
    
    # M-step
    m_step_rst <- M_step(x, W, y_expect_list, c_hat, r_hat, phi_hat, tau_hat)
    pi_hat <- m_step_rst$pi_hat
    c_hat <- m_step_rst$c_hat
    r_hat <- m_step_rst$r_hat
    phi_hat <- m_step_rst$phi_hat
    tau_hat <- m_step_rst$tau_hat
    sigma_sq <- m_step_rst$sigma_sq
    
    # convergence check
    if (iter > 1 && abs(log_likelihood_vals[iter] - log_likelihood_vals[iter - 1]) < tol) {
      message("Converged at iteration ", iter)
      break
    }
    
  }
  
  return(list(pi = pi_hat, c = c_hat, r = r_hat, phi = phi_hat, 
              tau = tau_hat, sigma_sq = sigma_sq, log_likelihood = log_likelihood_vals, W = W))
  
}
