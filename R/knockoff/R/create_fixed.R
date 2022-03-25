#' Fixed-X knockoffs
#' 
#' Creates fixed-X knockoff variables.
#' 
#' @param X normalized n-by-p matrix of original variables.(\eqn{n \geq p}).
#' @param intercept Will the intercept be fitted (TRUE) or set to zero
#' (default=FALSE).
#' @param method either "equi" or "sdp" (default: "sdp").
#' This determines the method that will be used to minimize the correlation between the original variables and the knockoffs.
#' @param sigma the noise level, used to augment the data with extra rows if necessary (default: NULL).
#' @param y vector of length n, containing the observed responses. 
#' This is needed to estimate the noise level if the parameter \code{sigma} is not provided, 
#' in case \eqn{p \leq n < 2p}.
#' When \code{intercept=TRUE}, it is needed to estimate the intercept in case \eqn{p \leq n < 2p+1} (default: NULL).
#' @param randomize whether the knockoffs are constructed deterministically or randomized (default: F).
#' @return An object of class "knockoff.variables". This is a list 
#'  containing at least the following components:
#'  \item{X}{n-by-p matrix of original variables (possibly augmented or transformed).}
#'  \item{Xk}{n-by-p matrix of knockoff variables.}
#'  \item{y}{vector of observed responses (possibly augmented). }
#' 
#' @family create
#' 
#' @references 
#'   Barber and Candes,
#'   Controlling the false discovery rate via knockoffs. 
#'   Ann. Statist. 43 (2015), no. 5, 2055--2085.
#' 
#' @details
#' Fixed-X knockoffs assume a homoscedastic linear regression model for \eqn{Y|X}. Moreover, they only guarantee
#' FDR control when used in combination with statistics satisfying the "sufficiency" property. 
#' In particular, the default statistics based on the cross-validated lasso does not satisfy this 
#' property and should not be used with fixed-X knockoffs.
#' 
#' @examples
#' set.seed(2022)
#' p=100; n=200; k=15
#' X = matrix(rnorm(n*p),n)
#' nonzero = sample(p, k)
#' beta = 5.5 * (1:p %in% nonzero)
#' y = X %*% beta + rnorm(n)
#' 
#' # Basic usage with default arguments
#' result = knockoff.filter(X, y, knockoffs=create.fixed)
#' print(result$selected)
#' 
#' # Advanced usage with custom arguments
#' knockoffs = function(X) create.fixed(X, method='equi')
#' result = knockoff.filter(X, y, knockoffs=knockoffs)
#' print(result$selected) 
#' 
#' @export
create.fixed <- function(X, intercept=F, method=c('sdp','equi'),
                         sigma=NULL, y=NULL, randomize=F) {
  method = match.arg(method)
  
  # Validate dimensions, if using fixed-X knockoffs
  n = nrow(X); p = ncol(X)
  if (n <= p+intercept){
    stop(paste0('Input X must have dimensions n > p',
                ifelse(intercept, "+1. ", ". ")))
  } else if (n < 2*p+intercept) {
    warning(paste0('Input X has dimensions p',
                   ifelse(intercept, "+1 ", " "),
                   '< n < 2p',
                   ifelse(intercept, "+1. ", ". ")),
            'Augmenting the model with extra rows.', immediate.=T)
    
    if( intercept && is.null(y) ) {
      stop('The response variables "y" must be provided in order to augment the 
           data with extra rows.')
    }
    fitted_intercept <- if(intercept){
      unname(lm(y ~ X + 1)$coefficients[1])
    } else{ 0 }
    
    X.svd = svd(X, nu=n, nv=0)
    u2 = X.svd$u[,(p+1):n]
    X = rbind(X, matrix(0, 2*p+intercept-n, p))
    
    if( is.null(sigma) ) {
      if( is.null(y) ) {
        stop('Either the noise level "sigma" or the response variables "y" must
             be provided in order to augment the data with extra rows.')
      }
      else{
        sigma = sqrt(mean((t(u2) %*% y)^2)) # = sqrt(RSS/(n-p))
      }
    }
    if (randomize)
      y.extra = rnorm(2*p-n+intercept, sd=sigma)
    else
      y.extra = with_seed(0, rnorm(2*p-n+intercept, sd=sigma))
    y = c(y, y.extra + fitted_intercept)
  }
  # Normalize X and y, if using fixed-X knockoffs
  X = normc(X, center = intercept)
  if( !is.null(y) ) {
    y = scale(y, center = intercept, scale = F)
  }
  
  Xk = switch(match.arg(method), 
               "equi" = create_equicorrelated(X,intercept,randomize),
               "sdp"  = create_sdp(X,intercept,randomize)
              )
  structure(list(X=X, Xk=Xk, y=y), class='knockoff.variables')
}

#' Create equicorrelated fixed-X knockoffs.
#'  
#' @rdname create_equicorrelated
#' @keywords internal
create_equicorrelated <- function(X, intercept, randomize) {
  # Compute SVD and U_perp.
  X.svd = decompose(X, intercept, randomize)
  
  # Set s = min(2 * smallest eigenvalue of X'X, 1), so that all the correlations
  # have the same value 1-s.
  if (any(X.svd$d <= 1e-5 * max(X.svd$d)))
    stop(paste('Data matrix is rank deficient.',
               'Equicorrelated knockoffs will have no power.'))
  lambda_min = min(X.svd$d)^2
  s = min(2*lambda_min, 1)
  
  # Construct the knockoff according to Equation 1.4.
  s_diff = pmax(0, 2*s - (s/X.svd$d)^2) # can be negative due to numerical error
  X_ko = (X.svd$u %*diag% (X.svd$d - s / X.svd$d) +
          X.svd$u_perp %*diag% sqrt(s_diff)) %*% t(X.svd$v)
}

#' Create SDP fixed-X knockoffs.
#'  
#' @rdname create_sdp
#' @keywords internal
create_sdp <- function(X, intercept, randomize) {
  # Compute SVD and U_perp.
  X.svd = decompose(X, intercept, randomize)
  
  # Check for rank deficiency.
  tol = 1e-5
  d = X.svd$d
  d_inv = 1 / d
  d_zeros = d <= tol*max(d)
  if (any(d_zeros)) {
    warning(paste('Data matrix is rank deficient.',
                  'Model is not identifiable, but proceeding with SDP knockoffs'),immediate.=T)
    d_inv[d_zeros] = 0
  }
  
  # Compute the Gram matrix and its (pseudo)inverse.
  G = (X.svd$v %*diag% d^2) %*% t(X.svd$v)
  G_inv = (X.svd$v %*diag% d_inv^2) %*% t(X.svd$v)
  
  # Optimize the parameter s of Equation 1.3 using SDP.
  s = create.solve_sdp(G)
  s[s <= tol] = 0
  
  # Construct the knockoff according to Equation 1.4:
  C.svd = canonical_svd(2*diag(s) - (s %diag*% G_inv %*diag% s))
  X_ko = X - (X %*% G_inv %*diag% s) + 
    (X.svd$u_perp %*diag% sqrt(pmax(0, C.svd$d))) %*% t(C.svd$v)
}


#' Compute the SVD of X and construct an orthogonal matrix U_perp such that U_perp * U = 0.
#'  
#' @rdname decompose
#' @keywords internal
decompose <- function(X, intercept, randomize) {
  n = nrow(X); p = ncol(X)
  stopifnot(n >= 2*p+intercept)
  
  result = knockoff:::canonical_svd(X)
  X_basis <- if(intercept){
    cbind(result$u, rep(1, n))
  } else{ result$u }
  Q = qr.Q(qr(cbind(X_basis, matrix(0,n,p))))
  u_perp = Q[,((p+1):(2*p))+intercept]
  if (randomize) {
    Q = qr.Q(qr(knockoff:::rnorm_matrix(p,p)))
    u_perp = u_perp %*% Q
  }
  result$u_perp = u_perp
  result
}