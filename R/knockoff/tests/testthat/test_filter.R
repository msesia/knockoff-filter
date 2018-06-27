test_that('knockoff.filter verifies input dimensions', {
  expect_error(knockoff.filter(knockoff:::rnorm_matrix(10, 10), rnorm(10), knockoffs=create.fixed), 'dimensions')
  expect_error(knockoff.filter(knockoff:::rnorm_matrix(20, 10), rnorm(19), knockoffs=create.fixed))
  expect_warning(knockoff.filter(knockoff:::rnorm_matrix(20, 15), rnorm(20), knockoffs=create.fixed), 'dimensions')
})

test_that('knockoff.filter for fixed design is invariant under permutations of the columns of the design matrix.', {
  # Problem parameters
  n = 250          # number of observations
  p = 100          # number of variables
  k = 30           # number of variables with nonzero coefficients
  amplitude = 5    # signal amplitude (for noise level = 1)
  
  # Generate the variables from a multivariate normal distribution
  X = matrix(rnorm(n*p),n)
  
  # Generate the response from a linear model
  nonzero = sample(p, k)
  beta = amplitude * (1:p %in% nonzero) / sqrt(n)
  y.sample <- function(X) X %*% beta + rnorm(n)
  y = y.sample(X)
  
  # Select variables
  set.seed(123)
  S = knockoff.filter(X, y, knockoffs=create.fixed)$selected
  
  # Permute columns
  idx_perm = sample(p)
  X_perm = X[,idx_perm]
  set.seed(123)
  S_perm = idx_perm[knockoff.filter(X_perm, y, knockoffs=create.fixed)$selected]
  
  # Verify consistency
  expect_true(setequal(S,S_perm))
})