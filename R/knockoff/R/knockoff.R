#' knockoff: A package for controlled variable selection
#'
#' This package implements the Knockoff Filter, which is a powerful and versatile tool for 
#' controlled variable selection.
#' 
#' @section Outline:
#' The procedure is based on the contruction of artificial 'knockoff copies' of the variables 
#' present in the given statistical model. Then, it selects those variables that are clearly better 
#' than their corresponding knockoffs, based on some measure of variable importance.
#' A wide range of statistics and machine learning tools can be exploited to estimate the 
#' importance of each variable, while guaranteeing finite-sample control of the false
#' discovery rate (FDR).
#' 
#' The Knockoff Filter controls the FDR in either of two statistical scenarios:
#' \itemize{
#'  \item{The "model-X" scenario: }{the response \eqn{Y} can depend on the variables \eqn{X=(X_1,\ldots,X_p)}
#'  in an arbitrary and unknown fashion, but the distribution of \eqn{X} must be known. In thise case
#'  there are no constraints on the dimensions \eqn{n} and \eqn{p} of the problem.}
#'  \item{The "fixed-X" scenario: }{the response \eqn{Y} depends upon \eqn{X} through a 
#'  homoscedastic Gaussian linear model and the problem is low-dimensional (\eqn{n \geq p}). 
#'  In this case, no modeling assumptions on \eqn{X} are required. }
#' }
#' 
#' For more information, see the website below and the accompanying paper.
#' 
#' \url{https://web.stanford.edu/group/candes/knockoffs/index.html}
#' 
#' @docType package
#' @name knockoff
NULL
