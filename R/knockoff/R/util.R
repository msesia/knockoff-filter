# Fast versions of diag(d) %*% X and X %*% diag(d).
`%diag*%` <- function(d, X) d * X
`%*diag%` <- function(X, d) t(t(X) * d)

# Efficient test for matrix positive-definiteness
# 
# Computes the smallest eigenvalue of a matrix A to verify whether
# A is positive-definite
#' @keywords internal
is_posdef = function(A, tol=1e-9) {
  p = nrow(matrix(A))
  
  if (p<500) {
    lambda_min = min(eigen(A)$values)
  }
  else {
    oldw <- getOption("warn")
    options(warn = -1)
    lambda_min = RSpectra::eigs(A, 1, which="SM", opts=list(retvec = FALSE, maxitr=100, tol))$values
    options(warn = oldw)
    if( length(lambda_min)==0 ) {
      # RSpectra::eigs did not converge. Using eigen instead."
      lambda_min = min(eigen(A)$values)
    }
  }
  return (lambda_min>tol*10)
}
  
# Reduced SVD with canonical sign choice.
# 
# Our convention is that the sign of each vector in U is chosen such that the
# coefficient with the largest absolute value is positive.
#' @keywords internal
canonical_svd = function(X) {
  X.svd = tryCatch({
    svd(X)
  }, warning = function(w){}, error = function(e) {
      stop("SVD failed in the creation of fixed-design knockoffs. Try upgrading R to version >= 3.3.0")
  }, finally = {})
  
  for (j in 1:min(dim(X))) {
    i = which.max(abs(X.svd$u[,j]))
    if (X.svd$u[i,j] < 0) {
      X.svd$u[,j] = -X.svd$u[,j]
      X.svd$v[,j] = -X.svd$v[,j]
  }
    }
  return(X.svd)
}

# Scale the columns of a matrix to have unit norm.
#' @keywords internal
normc = function(X,center=T) {
  X.centered = scale(X, center=center, scale=F)
  X.scaled = scale(X.centered, center=F, scale=sqrt(colSums(X.centered^2)))
  X.scaled[,] # No attributes
}

# Generate a random matrix with i.i.d. normal entries.
#' @keywords internal
rnorm_matrix = function(n, p, mean=0, sd=1) {
  matrix(rnorm(n*p, mean, sd), nrow=n, ncol=p)
}

# Generate a random, sparse regression problem.
#' @keywords internal
random_problem = function(n, p, k=NULL, amplitude=3) {
  if (is.null(k)) k = max(1, as.integer(p/5))
  X = normc(rnorm_matrix(n, p))
  nonzero = sample(p, k)
  beta = amplitude * (1:p %in% nonzero)
  y.sample <- function() X %*% beta + rnorm(n)
  list(X = X, beta = beta, y = y.sample(), y.sample = y.sample)
}

# Evaluate an expression with the given random seed, then restore the old seed.
#' @keywords internal
with_seed = function(seed, expr) {
  seed.old = if (exists('.Random.seed')) .Random.seed else NULL
  set.seed(seed)
  on.exit({
    if (is.null(seed.old)) {
      if (exists('.Random.seed')) 
        rm(.Random.seed, envir=.GlobalEnv)
    } else {
      .Random.seed <<- seed.old
    }
  })
  expr
}