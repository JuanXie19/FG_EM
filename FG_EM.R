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
  bessel_term <- besselI(tau, nu = (d/2 - 1), expon.scaled = TRUE)
  
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



newton_tau <- function(tau_init, phi_hat, y_expect_list, W){
  for (k in 1: M){
    NUMERATOR <- numeric(n)
    DENUMERATOR <- numeric(n)
    phi_k <- phi_hat[k, ]
    d <- length(phi_k)
    y_expect_k <- y_expect_list[[k]]
    
    for (i in 1:n){
      Ad_tau_k <- besselI(tau_init[k], nu = d/2, expon.scaled = TRUE)/besselI(tau_init[k], d/2 - 1, expon.scaled = TRUE)
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
  
  tau_update <- newton_tau(tau_hat, phi_hat, y_expect_list, W)
  tau_hat <- tau_update
	  
  
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
  pi_hat <- rep(1/M, M) # mixture weight
  
  set.seed(12345)
  init_kmeans <- kmeans(x, centers = M, nstart = 1000, algorithm = "Lloyd")
  c_hat <- init_kmeans$centers  # Cluster centers
  
  r_hat <- numeric(M)
  
  for (i in 1: M){
    cluster_points <- x[init_kmeans$cluster == i, , drop = FALSE]
    avg_distance <- mean(apply(cluster_points, 1, function(row) sqrt(sum((row - c_hat[i,])^2))))
    r_hat[i] <- avg_distance
  }
  
  # Initialize sigma_sq based on the squared residual
  
  residuals <- numeric(n)
  for (i in 1:n) {
    cluster_center <- init_kmeans$centers[init_kmeans$cluster[i], ]
    residuals[i] <- sum((x[i, ] - cluster_center)^2)
  }
  
  sigma_sq <- mean(residuals)
  
  
# Recover y_i from x_i
	y <- matrix(0, nrow = nrow(x), ncol = ncol(x))
	for (i in 1:nrow(x)) {
		cluster_idx <- init_kmeans$cluster[i]
		y_i <- (x[i, ] - c_hat[cluster_idx, ]) / r_hat[cluster_idx]
		y[i, ] <- y_i / sqrt(sum(y_i^2))  # Normalize y_i
	}

# Initialize phi_k (mean direction) and tau_k (concentration parameter)
	phi_hat <- matrix(0, nrow = M, ncol = ncol(x))
	tau_hat <- numeric(M)

	for (k in 1:M) {
		cluster_indices <- which(init_kmeans$cluster == k)
		if (length(cluster_indices) > 0) {
		# Compute mean direction phi_k
			y_cluster <- y[cluster_indices, , drop = FALSE]
			phi_hat[k, ] <- colSums(y_cluster)
			phi_hat[k, ] <- phi_hat[k, ] / sqrt(sum(phi_hat[k, ]^2))
		}
	}

  # Initialize tau_hat based on the angular spread of points in the cluster
  
	tau_hat <- numeric(M)
	for (k in 1:M) {
		# Get data points in cluster k
		cluster_points <- x[init_kmeans$cluster == k, ]
    
		# Skip if the cluster has no points or only one point
		if (nrow(cluster_points) <= 1) {
			tau_hat[k] <- 1  # Default value for small clusters
			next
		}
    
		phi_k <- phi_hat[k,]
   
		# Compute angular distances for all points in the cluster
		dot_products <- rowSums(cluster_points * phi_k)  # Dot product between each point and phi_k
		norms <- sqrt(rowSums(cluster_points^2))  # Norm of each point
		angular_distances <- acos(pmin(pmax(dot_products / norms, -1), 1))  # Ensure acos input is in [-1, 1]
    
		# Estimate tau_k as the inverse of the average angular distance
		# Use median for robustness to outliers
		tau_hat[k] <- 1 / median(angular_distances)
	}
  
  
  log_likelihood_vals <- numeric(max_iter)
  
  
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
              tau = tau_hat, sigma_sq = sigma_sq, log_likelihood = log_likelihood_vals))
  
}
