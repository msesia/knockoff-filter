test_that('Fixed-design equicorrelated knockoffs have the right correlation structure', {
  n = 20; p = 10
  X = knockoff:::normc(knockoff:::rnorm_matrix(n,p))
  knock_variables_default = create.fixed(X, method='equi', randomize=F)
  knock_variables_randomized = create.fixed(X, method='equi', randomize=T)
  X = knock_variables_default$X
  Xko_default = knock_variables_default$Xk
  Xko_randomized = knock_variables_randomized$Xk
  
  G = t(X) %*% X
  s = min(2*min(eigen(G)$values), 1)
  for (Xko in list(Xko_default, Xko_randomized)) {
    expect_equal(t(Xko) %*% Xko, G)
    expect_equal(t(X) %*% Xko, G - diag(s,p,p))
  }
})

# Test case from Weijie Su.
test_that('equicorrelated knockoffs are created in numerically sensitive case', {
  n = 15; p = 5
  M = matrix(0, p, p)
  diag(M) = 1
  for (i in 1:p) {
    for (j in 1:p) {
      if ((i==j+1) || (j==i+1))
        M[i,j] <- 0.6
      if ((i==j+2) || (j==i+2))
        M[i,j] <- 0.1
    }
  }
  X = knockoff:::with_seed(2, matrix(rnorm(n*p),n) %*% chol(M) )
  k = 4
  
  Z = knockoff:::normc(X[,-k])
  Z_ko = create.fixed(Z, method='equi', randomize=F)$Xk
  expect_false(any(is.nan(Z_ko)))
})

test_that('Fixed-design SDP knockoffs have the right correlation structure', {
  skip_on_cran()
  
  n = 20; p = 10
  X = knockoff:::normc(knockoff:::rnorm_matrix(n,p))
  knock_variables_default = create.fixed(X, method='sdp', randomize=F)
  knock_variables_randomized = create.fixed(X, method='sdp', randomize=T)
  X = knock_variables_default$X
  Xko_default = knock_variables_default$Xk
  Xko_randomized = knock_variables_randomized$Xk
  
  offdiag <- function(A) A - diag(diag(A))
  G = t(X) %*% X
  tol = 1e-4
  for (Xko in list(Xko_default, Xko_randomized)) {
    expect_equal(t(Xko) %*% Xko, G, tolerance=tol)
    expect_equal(offdiag(t(X) %*% Xko), offdiag(G), tolerance=tol)
    expect_true(all(diag(t(X) %*% Xko) < 1+tol))
  }
})

test_that('Gaussian equicorrelated knockoffs have the right correlation structure', {
  # Problem parameters
  n = 10000000   # number of observations
  p = 3          # number of variables
  
  # Generate the variables from a multivariate normal distribution
  mu = c(1,2,3); Sigma = matrix(c(1,0.55,0.2, 0.55,1,0.55, 0.2, 0.55, 1),3)
  
  X = matrix(rep(mu,each=n),n) + matrix(rnorm(n*p),n) %*% chol(Sigma)
  Xk = create.gaussian(X, mu, Sigma, method='equi')
  
  SigmaHat = cov(Xk)
  SigmaHatCross = cov(X, y=Xk)
  muHat = colMeans(Xk)
  
  lambda_min = eigen(Sigma, symmetric=T, only.values = T)$values[p]
  diag_s = diag(rep(1, nrow(Sigma)) * min(2*lambda_min, min(diag(Sigma))))
  
  expect_equal(mu, muHat, tolerance=2e-3)
  expect_equal(Sigma, SigmaHat, tolerance=2e-3)
  expect_equal(Sigma-diag_s, SigmaHatCross, tolerance=2e-3)
})

test_that('Fixed-design SDP knockoffs (with intercept adjustment) have the right correlation structure', {
  skip_on_cran()
  
  n = 21; p = 10
  X = knockoff:::normc(knockoff:::rnorm_matrix(n,p))
  knock_variables_default = create.fixed(X, intercept=T, method='sdp', randomize=F)
  knock_variables_randomized = create.fixed(X, intercept=T, method='sdp', randomize=T)
  X = knock_variables_default$X
  Xko_default = knock_variables_default$Xk
  Xko_randomized = knock_variables_randomized$Xk
  
  offdiag <- function(A) A - diag(diag(A))
  G = t(X) %*% X
  tol = 1e-4
  for (Xko in list(Xko_default, Xko_randomized)) {
    expect_equal(t(Xko) %*% Xko, G, tolerance=tol)
    expect_equal(offdiag(t(X) %*% Xko), offdiag(G), tolerance=tol)
    expect_equal(colSums(X), colSums(Xko), tolerance=tol)
    expect_true(all(diag(t(X) %*% Xko) < 1+tol))
  }
})

test_that('Fixed-design equicorrelated knockoffs (with intercept adjustment) have the right correlation structure', {
  n = 21; p = 10
  X = knockoff:::normc(knockoff:::rnorm_matrix(n,p))
  knock_variables_default = create.fixed(X, intercept=T, method='equi', randomize=F)
  knock_variables_randomized = create.fixed(X, intercept=T, method='equi', randomize=T)
  X = knock_variables_default$X
  Xko_default = knock_variables_default$Xk
  Xko_randomized = knock_variables_randomized$Xk
  
  G = t(X) %*% X
  s = min(2*min(eigen(G)$values), 1)
  for (Xko in list(Xko_default, Xko_randomized)) {
    expect_equal(t(Xko) %*% Xko, G)
    expect_equal(t(X) %*% Xko, G - diag(s,p,p))
    expect_equal(colSums(X), colSums(Xko))
  }
})