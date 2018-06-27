#' Relaxed optimization for fixed-X and Gaussian knockoffs
#'
#' This function solves the optimization problem needed to create fixed-X and Gaussian SDP knockoffs
#' on a block-diagonal approximation of the covariance matrix. This will be less
#' powerful than \code{\link{create.solve_sdp}}, but more computationally efficient.
#' 
#' @param Sigma positive-definite p-by-p covariance matrix.
#' @param nBlocks number of blocks in the block-diagonal approximation of Sigma (default: 10).
#' @param cores number of cores used to solve the smaller SDPs (default: 1).
#' @param maxit the maximum number of iterations for the solver (default: 1000).
#' @param gaptol tolerance for duality gap as a fraction of the value of the objective functions (default: 1e-6).
#' @return The solution \eqn{s} to the semidefinite programming problem defined above.
#'
#' @details Solves the following two-step semidefinite programming problem:
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
#' If the matrix Sigma supplied by the user is a non-scaled covariance matrix 
#' (i.e. its diagonal entries are not all equal to 1), then the appropriate scaling is applied before
#' solving the SDP defined above. The result is then scaled back before being returned, as to match 
#' the original scaling of the covariance matrix supplied by the user.
#' 
#' @family optimization
#' 
#' @export
create.solve_asdp <- function(Sigma, nBlocks=10, cores=1, gaptol=1e-6, maxit=1000) {
  # Check that covariance matrix is symmetric
  stopifnot(isSymmetric(Sigma))
  
  # Approximate the covariance matrix as block diagonal
  cluster_sol = divide_sdp(Sigma, nBlocks=nBlocks)

  # Solve the smaller SDPs corresponding to each block
  s_asdp_list = sapply(1:nBlocks, function(i) {
  create.solve_sdp(as.matrix(cluster_sol$subSigma[[i]]), gaptol=gaptol, maxit=maxit) })

  # Assemble the solutions into one vector of length p
  p = dim(Sigma)[1]
  idx_count = rep(1, nBlocks)
  s_asdp = rep(0,p)
  for( j in 1:p ){
    cluster_j = cluster_sol$clusters[j]
    s_asdp[j] = s_asdp_list[[cluster_j]][idx_count[cluster_j]]
    idx_count[cluster_j] = idx_count[cluster_j]+1
  }
  
  # Maximize the shrinkage factor
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
  
  # Verify that the solution is correct
  if (!is_posdef(2*Sigma-diag(s_asdp_scaled,length(s_asdp_scaled)))) {
    warning('In creation of approximate SDP knockoffs, procedure failed. Knockoffs will have no power.',immediate.=T)
    s_asdp_scaled = 0*s_asdp_scaled
  }
  
  # Return result
  s_asdp_scaled
}

#' Approximate a covariance matrix by a block diagonal matrix with blocks
#' of approximately equal size using Ward's method for hierarchical clustering
#'  
#' @rdname divide_sdp
#' @keywords internal
divide_sdp <- function(Sigma, nBlocks=10) {
  
  # Convert the covariance matrix into a dissimilarity matrix
  # Add a small perturbation to stabilize the clustering in the case of highly symmetrical matrices
  p = ncol(Sigma)
  Eps = matrix(rnorm(p*p),p)*1e-6
  dissimilarity = 1 - cov2cor(Sigma+Eps)
  distance = as.dist(dissimilarity)
  
  # Hierarchical clustering wiht Ward's method for clusters of roughly equal size
  fit = hclust(distance, method="ward.D2")
  # Cut tree into nBlocks clusters
  clusters = cutree(fit, k=nBlocks)

  # Create covariance submatrices for each cluster
  subSigma = vector("list", nBlocks) 
  for( k in 1:nBlocks ) {
    indices_k = clusters==k
    subSigma[[k]] = Sigma[indices_k,indices_k]
  }
  
  # Return the cluster assignments and the cluster covariance submatrices
  structure(list(clusters=clusters, subSigma=subSigma), class='knockoff.clusteredCovariance')
}