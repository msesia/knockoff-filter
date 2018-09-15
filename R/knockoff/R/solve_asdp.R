#' Relaxed optimization for fixed-X and Gaussian knockoffs
#'
#' This function solves the optimization problem needed to create fixed-X and Gaussian SDP knockoffs
#' on a block-diagonal approximation of the covariance matrix. This will be less
#' powerful than \code{\link{create.solve_sdp}}, but more computationally efficient.
#' 
#' @param Sigma positive-definite p-by-p covariance matrix.
#' @param max.size size of the largest block in the block-diagonal approximation of Sigma (default: 500). See Details.
#' @param maxit the maximum number of iterations for the solver (default: 1000).
#' @param gaptol tolerance for duality gap as a fraction of the value of the objective functions (default: 1e-6).
#' @param verbose whether to display progress (default: FALSE).
#' @return The solution \eqn{s} to the semidefinite program defined above.
#'
#' @details Solves the following two-step semidefinite program:
#'
#'   (step 1)  \deqn{ \mathrm{maximize}     \; \mathrm{sum}(s) \quad
#'                    \mathrm{subject} \; \mathrm{to:}  \; 0 \leq s \leq 1, \;
#'                                          2 \Sigma_{\mathrm{approx}} - \mathrm{diag}(s) \geq 0}
#'                              
#'   (step 2) \deqn{ \mathrm{maximize}      \; \gamma \quad
#'                   \mathrm{subject} \; \mathrm{to:}    \; \mathrm{diag}(\gamma s) \leq 2 \Sigma}
#'
#' Each smaller SDP is solved using the interior-point method implemented in \code{\link[Rdsdp]{dsdp}}.
#' 
#' The parameter max.size controls the size of the largest semidefinite program that needs to be solved.
#' A larger value of max.size will increase the computation cost, while yielding a solution closer to
#' that of the original semidefinite program.
#'
#' If the matrix Sigma supplied by the user is a non-scaled covariance matrix 
#' (i.e. its diagonal entries are not all equal to 1), then the appropriate scaling is applied before
#' solving the SDP defined above. The result is then scaled back before being returned, as to match 
#' the original scaling of the covariance matrix supplied by the user.
#' 
#' @family optimization
#' 
#' @export
create.solve_asdp <- function(Sigma, max.size=500, gaptol=1e-6, maxit=1000, verbose=FALSE) {
  # Check that covariance matrix is symmetric
  stopifnot(isSymmetric(Sigma))
  
  if(ncol(Sigma) <= max.size) return(create.solve_sdp(Sigma, gaptol=gaptol, maxit=maxit, verbose=verbose))
  
  # Approximate the covariance matrix as block diagonal
  if(verbose) cat(sprintf("Dividing the problem into subproblems of size <= %s ... ", max.size))
  cluster_sol = divide.sdp(Sigma, max.size=max.size)
  n.blocks = max(cluster_sol$clusters)
  if(verbose) cat("done. \n")
  
  # Solve the smaller SDPs corresponding to each block
  if(verbose) cat(sprintf("Solving %s smaller SDPs ... \n", n.blocks))
  s_asdp_list = list()
  if(verbose) pb <- utils::txtProgressBar(min = 0, max = n.blocks, style = 3)
  for(k in 1:n.blocks) {
    s_asdp_list[[k]] = create.solve_sdp(as.matrix(cluster_sol$subSigma[[k]]), gaptol=gaptol, maxit=maxit)
    if(verbose) utils::setTxtProgressBar(pb, k)
  }
  if(verbose) cat("\n")
  
  # Assemble the solutions into one vector of length p
  p = dim(Sigma)[1]
  idx_count = rep(1, n.blocks)
  s_asdp = rep(0,p)
  for( j in 1:p ){
    cluster_j = cluster_sol$clusters[j]
    s_asdp[j] = s_asdp_list[[cluster_j]][idx_count[cluster_j]]
    idx_count[cluster_j] = idx_count[cluster_j]+1
  }
  
  # Maximize the shrinkage factor
  if(verbose) cat(sprintf("Combinining the solutions of the %s smaller SDPs ... ", n.blocks))
  tol = 1e-12
  maxitr=100000
  gamma_range = seq(0,1,len=1000)
  options(warn=-1)
  gamma_opt = gtools::binsearch( function(i) {
    G = 2*Sigma - gamma_range[i]*diag(s_asdp)
    lambda_min = RSpectra::eigs(G, 1, which = "SR", opts = list(retvec = FALSE, maxitr=maxitr, tol=tol))$values
    if (length(lambda_min)==0) {
      lambda_min = 1  # Not converged
    }
    lambda_min
  }, range=c(1,length(gamma_range)) )
  s_asdp_scaled = gamma_range[min(gamma_opt$where)]*s_asdp
  options(warn=0)
  if(verbose) cat("done. \n")
  
  if(verbose) cat("Verifying that the solution is correct ... ")
  # Verify that the solution is correct
  if (!is_posdef(2*Sigma-diag(s_asdp_scaled,length(s_asdp_scaled)))) {
    warning('In creation of approximate SDP knockoffs, procedure failed. Knockoffs will have no power.',immediate.=T)
    s_asdp_scaled = 0*s_asdp_scaled
  }
  if(verbose) cat("done. \n")
  
  # Return result
  s_asdp_scaled
}

#' Merge consecutive clusters of correlated variables while ensuring 
#' that no cluster has size larger than max.size
#'  
#' @rdname merge.clusters
#' @keywords internal
merge.clusters <- function(clusters, max.size) {
  cluster.sizes = table(clusters)
  clusters.new = rep(0, length(clusters))
  g = 1
  g.size = 0
  for(k in 1:max(clusters)) {
    if(g.size + cluster.sizes[k] > max.size) {
      g = g + 1
      g.size = 0
    }
    clusters.new[clusters==k] = g
    g.size = g.size + cluster.sizes[k]
  }
  return(clusters.new)
}

#' Approximate a covariance matrix by a block diagonal matrix with blocks
#' of approximately equal size using Ward's method for hierarchical clustering
#'  
#' @rdname divide.sdp
#' @keywords internal
divide.sdp <- function(Sigma, max.size) {
  # Convert the covariance matrix into a dissimilarity matrix
  # Add a small perturbation to stabilize the clustering in the case of highly symmetrical matrices
  p = ncol(Sigma)
  Eps = matrix(rnorm(p*p),p)*1e-6
  dissimilarity = 1 - abs(cov2cor(Sigma)+Eps)
  distance = as.dist(dissimilarity)
  
  # Hierarchical clustering
  fit = hclust(distance, method="single")
  # Cut tree into clusters of size smaller than a threshold
  n.blocks.min = 1
  n.blocks.max = ncol(Sigma)
  for(it in 1:100) {
    n.blocks = ceiling((n.blocks.min+n.blocks.max)/2)
    clusters = cutree(fit, k=n.blocks)
    size = max(table(clusters))
    if(size <= max.size) {
      n.blocks.max = n.blocks
    }
    if(size >= max.size) {
      n.blocks.min = n.blocks
    }
    if(n.blocks.min == n.blocks.max) {
      break
    }
  }
  
  # Merge small clusters    
  clusters.new = merge.clusters(clusters, max.size)
  while(sum(clusters.new != clusters)>0) {
    clusters = clusters.new
    clusters.new = merge.clusters(clusters, max.size)
  }
  clusters = clusters.new
  
  # Create covariance submatrices for each cluster
  subSigma = vector("list", max(clusters)) 
  for( k in 1:length(subSigma) ) {
    indices_k = clusters==k
    subSigma[[k]] = Sigma[indices_k,indices_k]
  }
  
  # Return the cluster assignments and the cluster covariance submatrices
  structure(list(clusters=clusters, subSigma=subSigma), class='knockoff.clusteredCovariance')
}