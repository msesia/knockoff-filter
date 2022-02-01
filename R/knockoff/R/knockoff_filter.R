#' @docType package
#' @name knockoff
#' @import stats methods
NULL

#' The Knockoff Filter
#' 
#' This function runs the Knockoffs procedure from start to finish, selecting variables
#' relevant for predicting the outcome of interest.
#' 
#' This function creates the knockoffs, computes the importance statistics, 
#' and selects variables. 
#' It is the main entry point for the knockoff package.
#' 
#' @param X n-by-p matrix or data frame of predictors.
#' @param y response vector of length n.
#' @param knockoffs method used to construct knockoffs for the \eqn{X} variables.
#' It must be a function taking a n-by-p matrix as input and returning a n-by-p matrix of knockoff variables. 
#' By default, approximate model-X Gaussian knockoffs are used.
#' @param statistic statistics used to assess variable importance. By default, 
#' a lasso statistic with cross-validation is used. See the Details section for more information.
#' @param fdr target false discovery rate (default: 0.1).
#' @param offset either 0 or 1 (default: 1). This is the offset used to compute the rejection threshold on the
#' statistics. The value 1 yields a slightly more conservative procedure ("knockoffs+") that
#' controls the false discovery rate (FDR) according to the usual definition, 
#' while an offset of 0 controls a modified FDR.
#'
#' @return An object of class "knockoff.result". This object is a list 
#'  containing at least the following components:
#'  \item{X}{matrix of original variables}
#'  \item{Xk}{matrix of knockoff variables}
#'  \item{statistic}{computed test statistics}
#'  \item{threshold}{computed selection threshold}
#'  \item{selected}{named vector of selected variables}
#'
#' @details
#' 
#' The parameter \code{knockoffs} controls how knockoff variables are created.
#' By default, the model-X scenario is assumed and a multivariate normal distribution 
#' is fitted to the original variables \eqn{X}. The estimated mean vector and the covariance 
#' matrix are used to generate second-order approximate Gaussian knockoffs.
#' In general, the function \code{knockoffs} should take a n-by-p matrix of
#' observed variables \eqn{X} as input and return a n-by-p matrix of knockoffs.
#' Two default functions for creating knockoffs are provided with this package.
#' 
#' In the model-X scenario, under the assumption that the rows of \eqn{X} are distributed 
#' as a multivariate Gaussian with known parameters, then the function 
#' \code{create.gaussian} can be used to generate Gaussian knockoffs, 
#' as shown in the examples below.
#' 
#' In the fixed-X scenario, one can create the knockoffs using the function 
#' \code{create.fixed}. This requires \eqn{n \geq p} and it assumes 
#' that the response \eqn{Y} follows a homoscedastic linear regression model.
#' 
#' For more information about creating knockoffs, type \code{??create}.
#' 
#' The default importance statistic is \link{stat.glmnet_coefdiff}.
#' For a complete list of the statistics provided with this package, 
#' type \code{??stat}.
#' 
#' It is possible to provide custom functions for the knockoff constructions 
#' or the importance statistics. Some examples can be found in the vignette.
#' 
#' @references 
#'   Candes et al., Panning for Gold: Model-free Knockoffs for High-dimensional Controlled Variable Selection,
#'   arXiv:1610.02351 (2016).
#'   \href{https://web.stanford.edu/group/candes/knockoffs/index.html}{https://web.stanford.edu/group/candes/knockoffs/index.html}
#'   
#'   Barber and Candes,
#'   Controlling the false discovery rate via knockoffs. 
#'   Ann. Statist. 43 (2015), no. 5, 2055--2085.
#'   \doi{10.1214/15-AOS1337}
#' 
#' @examples
#' set.seed(2022)
#' p=200; n=100; k=15
#' mu = rep(0,p); Sigma = diag(p)
#' X = matrix(rnorm(n*p),n)
#' nonzero = sample(p, k)
#' beta = 3.5 * (1:p %in% nonzero)
#' y = X %*% beta + rnorm(n)
#' 
#' # Basic usage with default arguments
#' result = knockoff.filter(X, y)
#' print(result$selected)
#' 
#' # Advanced usage with custom arguments
#' knockoffs = function(X) create.gaussian(X, mu, Sigma)
#' k_stat = function(X, Xk, y) stat.glmnet_coefdiff(X, Xk, y, nfolds=5)
#' result = knockoff.filter(X, y, knockoffs=knockoffs, statistic=k_stat)
#' print(result$selected)
#' 
#' 
#' @export
knockoff.filter <- function(X, y,
                              knockoffs=create.second_order,
                              statistic=stat.glmnet_coefdiff, 
                              fdr=0.10,
                              offset=1
                              ) {
  
  # Validate input types.
  if (is.data.frame(X)) {
    X.names = names(X)
    X = as.matrix(X, rownames.force = F)  
  } else if (is.matrix(X)) {
    X.names = colnames(X)
  } else {
    stop('Input X must be a numeric matrix or data frame')
  }
  if (!is.numeric(X)) stop('Input X must be a numeric matrix or data frame')
  
  if (!is.factor(y) && !is.numeric(y)) {
    stop('Input y must be either of numeric or factor type')
  }
  if( is.numeric(y) ) y = as.vector(y)
  
  if(offset!=1 && offset!=0) {
    stop('Input offset must be either 0 or 1')
  }
  
  if (!is.function(knockoffs)) stop('Input knockoffs must be a function')
  if (!is.function(statistic)) stop('Input statistic must be a function')
  
  # Validate input dimensions
  n = nrow(X); p = ncol(X)
  stopifnot(length(y) == n)

  # If fixed-design knockoffs are being used, provive them with the response vector
  # in order to augment the data with new rows if necessary
  if( identical(knockoffs, create.fixed) )
    knockoffs = function(x) create.fixed(x, y=y)
  
  # Create knockoff variables
  knock_variables = knockoffs(X)
  
  # If fixed-design knockoffs are being used, update X and Y with the augmented observations (if present)
  if (is(knock_variables,"knockoff.variables")){
    X  = knock_variables$X
    Xk = knock_variables$Xk
    if(!is.null(knock_variables$y)) y  = knock_variables$y
    rm(knock_variables)
  } else if (is(knock_variables,"matrix")){
    Xk = knock_variables
    rm(knock_variables)
  } else {
    stop('Knockoff variables of incorrect type')
  }
  
  # Compute statistics
  W = statistic(X, Xk, y)
  
  # Run the knockoff filter
  t = knockoff.threshold(W, fdr=fdr, offset=offset)
  selected = sort(which(W >= t))
  if (!is.null(X.names))
    names(selected) = X.names[selected]
  
  # Package up the results.
  structure(list(call = match.call(),
                 X = X,
                 Xk = Xk,
                 y = y,
                 statistic = W,
                 threshold = t,
                 selected = selected),
            class = 'knockoff.result')
}

#' Threshold for the knockoff filter
#' 
#' Computes the threshold for the knockoff filter.
#' 
#' @param W the test statistics
#' @param fdr target false discovery rate (default: 0.1)
#' @param offset either 0 or 1 (default: 1). The offset used to compute the rejection threshold on the
#' statistics. The value 1 yields a slightly more conservative procedure ("knockoffs+") that
#' controls the FDR according to the usual definition, while an offset of 0 controls a modified FDR.
#' @return The threshold for variable selection.
#' 
#' @export
knockoff.threshold <- function(W, fdr=0.10, offset=1) {
  if(offset!=1 && offset!=0) {
    stop('Input offset must be either 0 or 1')
  }
  ts = sort(c(0, abs(W)))
  ratio = sapply(ts, function(t)
    (offset + sum(W <= -t)) / max(1, sum(W >= t)))
  ok = which(ratio <= fdr)
  ifelse(length(ok) > 0, ts[ok[1]], Inf)
}

#' Print results for the knockoff filter
#' 
#' Prints the list of variables selected by the knockoff filter and the corresponding function call.
#' 
#' @param x the output of a call to knockoff.filter
#' @param ... unused
#' 
#' @method print knockoff.result
#' @export
print.knockoff.result <- function(x, ...) {
  cat('Call:\n')
  print(x$call)
  cat('\nSelected variables:\n')
  print(x$selected)
}

#' Verify dependencies for chosen statistics
#' 
#' @param statistic the statistic chosen by the user
#' 
#' @keywords internal
verify_stat_depends <- function(statistic) {
  
}