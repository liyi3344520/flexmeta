#' Bayesian meta-analysis using a Sinh-arcsinh prior
#'
#' @name sas_stan
#' @rdname sas_stan
#' @param data A list of the data
#' @param par parameters of interest
#' @param chains the number of Markov chains
#' @param iter the number of iterations for each chain including warmup
#' @param warmup the number of warmup iterations per chain
#' @param thin the period for saving samples
#' @param init A list of the initial values
#' @return
#' \itemize{
#' \item \code{f}: An object of S4 class \code{stanfit} representing the fitted results.
#' Slot \code{mode} for this object indicates if the sampling is done or not.
#' }
#' @examples
#' require("rstansas")
#' m1 <- c(15,12,29,42,14,44,14,29,10,17,38,19,21)
#' n1 <- c(16,16,34,56,22,54,17,58,14,26,44,29,38)
#' m2 <- c( 9, 1,18,31, 6,17, 7,23, 3, 6,12,22,19)
#' n2 <- c(16,16,34,56,22,55,15,58,15,27,45,30,38)
#' dat1 <- pimeta::convert_bin(m1, n1, m2, n2, type = "logOR")
#' dat2 <- list(K = dim(dat1)[1], yi = dat1$y, si = dat1$se)
#' init1 <- function() {
#'   list(sigma = runif(1, 0, 2), delta = runif(1, 0, 2),
#'        mu = runif(1, -2, 2), epsilon = runif(1, -2, 2),
#'        theta = runif(dim(dat1)[1], -2, 2), theta_new = runif(1, -2, 2))
#' }
#' sas_stan(data = dat2, init = init1)
#' @export
sas_stan <- function(data, par = NA, chains = 1, iter = 3500, warmup = 500, thin = 1, init = 'random') {

  if (is.na(par)) {
    par <- c("mu", "sigma", "epsilon", "delta", "theta_new", "x_post_mean", "x_post_var")
  }

  return(
    rstan::sampling(stanmodels$sas, data = data, par = par, chains = chains,
                    iter = iter, warmup = warmup, thin = thin, init = init)
  )

}
