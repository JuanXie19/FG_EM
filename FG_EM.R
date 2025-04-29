library(movMF)
library(mvtnorm)
library(matrixStats)
library(Directional)



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
    # The normalizing constant is C_d(0) = gamma(d/2)/(2*pi^(d/2))
    # So log(C_d(0)) = log(gamma(d/2)) - log(2) - (d/2)*log(pi)
    # lgamma computes the logarithm of the gamma function
    return(-log(2) - (d/2) * log(pi) + lgamma(d/2))
  }
  
  # use expon.scale =T to avoid overflow
  log_bessel <- log(besselI(tau, nu = (d/2 - 1), expon.scaled = TRUE))
  res <- (d/2 - 1) * log(tau) - log_bessel - (d/2) * log(2*pi)
  return(res)
}



# Function to calculates FG-kernel density given the parameters cntr(center), rd(radius), sigma.sq, phi(direction), tau(concentration)
# this comes from the FG paper code. Note that we 
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
    kappa_k <- pmax(kappa_k, 1e-10)  # Small value to avoid division by zero
    mu_k <- nu / kappa_k
    
   
    # Compute A_d(kappa_k) using exponentially scaled Bessel functions for numerical stability
    log_Ad_k <- log(besselI(kappa_k, d/2, expon.scaled = TRUE)) - 
      log(besselI(kappa_k, d/2 - 1, expon.scaled = TRUE))
    Ad_k <- exp(log_Ad_k)
    Ad_k <- pmin(pmax(Ad_k, 1e-10), 1 - 1e-10)
    
    # Compute E[y_i | x_i, z_i = k]
    y_expect_k <- Ad_k * mu_k
    
    E_y[[k]] <- y_expect_k
  }
  
  return(E_y)
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
  
  ## opt4: based on the formula(2) from the movMF package 
  for (k in 1:M){
    #rbar <- norm_vec(phi_hat[k])/pi_hat[k]  # should not use pi_hat[k], but sum(W[, k]) instead
    rbar <- norm_vec(colSums(W[, k] * y_expect_list[[k]])) / sum(W[, k])  # also not used phi_hat[k], as it is normalized, should use unnormalized phi_hat here
    temp1 <- (rbar * d - rbar^3)/(1 - rbar^2)
    tau_hat[k] <- pmax(temp1, 1e-6)

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
  
  log_likelihood_vals <- rep(NA, max_iter) # use NA to distinguish unconverged iterations
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
  
  converged_loglik <- log_likelihood_vals[1:iter]
  return(list(
    pi = pi_hat, c = c_hat, r = r_hat, 
    phi = phi_hat, tau = tau_hat, sigma_sq = sigma_sq, 
    log_likelihood = converged_loglik,
    n_iter = iter,
    W = W))
}


# Function to compute BIC for a given M
FG_EM_with_criteria <- function(x, M, max_iter, tol) {
  fit <- FG_EM(x, M, max_iter, tol)
  
  # Compute number of parameters (p)
  n <- nrow(x)
  d <- ncol(x)
  p <- (M-1) + M*d + M + M*d + M + 1  # π_k (M-1), c_k (M*d), r_k (M), φ_k (M*d), τ_k (M), σ² (1)
  
  # Use the LAST (converged) log-likelihood value
  final_loglik <- tail(fit$log_likelihood, 1)
  
  # Compute BIC and AIC
  BIC <- -2 * final_loglik + p * log(n)
  AIC <- -2 * final_loglik + 2 * p
  
  # Return model with BIC
  return(list(
    model = fit,
    AIC = AIC,
    BIC = BIC,
    n_parameters = p,
    converged_iter = fit$n_iter
  ))
}

# Updated model selection function

select_best_M <- function(x, M_candidates = 2:5, max_iter = 100, tol = 1e-4) {
  # Initialize storage
  results <- list()
  criteria_df <- data.frame(
    M = M_candidates, 
    AIC = NA_real_,
    BIC = NA_real_,
    converged = NA,
    n_iter = NA
  )
  
  # Fit models for all M candidates
  for (i in seq_along(M_candidates)) {
    M <- M_candidates[i]
    message("\nFitting M = ", M)
    result <- FG_EM_with_criteria(x, M, max_iter, tol)
    
    # Store full results (including model, AIC, BIC, etc.)
    results[[paste0("M_", M)]] <- result
    
    # Update criteria table
    criteria_df$AIC[i] <- result$AIC
    criteria_df$BIC[i] <- result$BIC
    criteria_df$converged[i] <- (result$model$n_iter < max_iter)
    criteria_df$n_iter[i] <- result$model$n_iter
  }
  
  # Find best M by AIC/BIC (optional, since we keep all results)
  best_AIC_idx <- which.min(criteria_df$AIC)
  best_BIC_idx <- which.min(criteria_df$BIC)
  
  # Return ALL models + summary table
  return(list(
    all_models = results,           # Full results for every M
    criteria_table = criteria_df,    # Summary table
    best_M_AIC = M_candidates[best_AIC_idx],
    best_M_BIC = M_candidates[best_BIC_idx],
    best_model_AIC = results[[best_AIC_idx]]$model,  # Optional
    best_model_BIC = results[[best_BIC_idx]]$model   # Optional
  ))
}




### Given the FG-mixture results and a number N, this function generates N samples from estimated density, and also provides scatter plots with and withour error of the same
plot_density_FG_EM <- function(N, rst){
  cntr <- rst$c
  rd <- rst$r
  phi <- rst$phi
  tau <- rst$tau
  Pi <- rst$pi
  sigma.sq <- rst$sigma_sq
  D <- ncol(cntr)
  M <- length(tau)
  
  z <- matrix(data = NA,nrow = N,ncol = D)
  
  
  kr_label <- sapply(1:N,function(u){v <- rmultinom(n = 1,size = 1,prob = Pi); return(which.max(v))})
  for(k in 1:M){
      idx <- which(kr_label == k)
      if(length(idx) > 0){
          yy <- rvmf(length(idx), phi[k, ],tau[k])
          if(length(idx) == 1){
            z[idx,] <- cntr[k, ] + rd[k] * yy
          }
          else{
            z[idx,] <- t(apply(yy, 1, function(u){return(cntr[k, ] + rd[k] * u)})) 
          }
      }
      
  }
  error <- mvrnorm(N, rep(0, D), sigma.sq * diag(D))
  w <- z + error
  if(D == 2){
    par(mfrow=c(1,2))
    plot(z[, 1], z[, 2],main = "Without error", lwd = 3,cex = .3)
    plot(w[, 1], w[, 2],main = "Predictive", lwd = 3, cex = .3)
  }
  if(D == 3){
    library(plot3D)
    scatter3D(z[, 1],z[, 2],z[, 3],main = "Without error",lwd = 3,cex = .3)
    scatter3D(w[, 1],w[, 2],w[, 3],main = "Predictive",lwd = 3,cex = .3)
  }
  return(list(z,w))
}
