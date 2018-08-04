test_that('Statistics obey antisymmetry property', {
  n = 10; p = 5;
  prob = random_problem(n, p)
  knock.variables = create.fixed(prob$X)
  X = knock.variables$X
  Xk = knock.variables$Xk
  G = cbind(X, Xk)
  y = prob$y
    
  i = sort(sample(p, sample(p))) # Indices to swap.
  G_swap = G
  G_swap[,c(i,i+p)] <- G[,c(i+p,i)]
  
  expect_antisymmetric <- function(stat) {
    orig = 1:p; ko = (p+1):(2*p);
    expect_equal(stat(G[,orig],G[,ko],y), 
                 stat(G_swap[,orig],G_swap[,ko],y) * ifelse(1:p %in% i, -1, 1),tolerance = 1e-3)
  }
  expect_antisymmetric(stat.forward_selection)
  stats_fs_omp = function(X,Xk,y) stat.forward_selection(X, Xk, y, omp=FALSE)
  expect_antisymmetric(stats_fs_omp)
  stats_lasso_diff = function(X,Xk,y) stat.lasso_lambdadiff(X, Xk, y, nlambda=100000)
  expect_antisymmetric(stats_lasso_diff)
  stats_lasso_signed_max = function(X,Xk,y) stat.lasso_lambdasmax(X, Xk, y, nlambda=100000)
  expect_antisymmetric(stats_lasso_signed_max)
})

test_that('Finding the max lambda in lasso works for orthonormal design', {
  n = 30; p = 10; amplitude = 3.5;
  X = qr.Q(qr(rnorm_matrix(n,p)))
  beta = amplitude * rnorm(p)
  y = X %*% beta + rnorm(n)
    
  beta_ls = as.vector(t(X) %*% y)
  expect_equal(lasso_max_lambda_glmnet(X, y, nlambda = 1e4, intercept=F, standardize=F), abs(beta_ls),
               tolerance = 1e-3)
})