## the initialization for the EM algorithm
# use k-kmeans to cluster the data, use cluster centers as the intial values for c_hat, 
initial_kmeans <- function(x, M){
  n <- nrow(x)
  pi_hat <- rep(1/M, M)
  set.seed(12345)
  init_kmeans <- kmeans(x, centers = M, nstart = 1000, algorithm = 'Lloyd')
  c_hat <- init_kmeans$centers
  r_hat <- numeric(M)
  for (k in 1: M){
    cluster_points <- x[init_kmeans$cluster == k, , drop = FALSE]
    avg_distance <- mean(apply(cluster_points, 1, function(row) sqrt(sum((row - c_hat[k,])^2))))
    r_hat[k] <- avg_distance
  }

  # initialize sigma_sq based on the squared residual
  residuals <- numeric(n)
  for (i in 1:n) {
    cluster_center <- init_kmeans$centers[init_kmeans$cluster[i], ]
    residuals[i] <- sum((x[i, ] - cluster_center)^2)
  }
  sigma_sq <- mean(residuals)

  # recover y from x
  y <- matrix(0, nrow = nrow(x), ncol = ncol(x))
	for (i in 1:nrow(x)) {
		cluster_idx <- init_kmeans$cluster[i]
		y_i <- (x[i, ] - c_hat[cluster_idx, ]) / r_hat[cluster_idx]
		y[i, ] <- y_i / sqrt(sum(y_i^2))  # Normalize y_i
	}
                               
  # initialize phi_hat (mean direction)
  phi_hat <- matrix(0, nrow = M, ncol = ncol(x))
	for (k in 1:M) {
		cluster_points <- which(init_kmeans$cluster == k)
		if (length(cluster_points) > 0) {
		# Compute mean direction phi_k
			y_cluster <- y[cluster_points, , drop = FALSE]
			phi_hat[k, ] <- colSums(y_cluster)
			phi_hat[k, ] <- phi_hat[k, ] / sqrt(sum(phi_hat[k, ]^2))
		}
	}
  # initialize tau_hat based on the angular spread of points in the cluster
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
  return(list(pi_hat = pi_hat, c_hat = c_hat, r_hat = r_hat, sigma_sq = sigma_sq, phi_hat = phi_hat, tau_hat = tau_hat))
}
