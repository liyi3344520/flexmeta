#' Generate Summary of Predictive Distribution
#'
#' @name pred_f
#' @rdname pred_f
#' @param fit An object of S4 class \code{stanfit} representing
#'   the fitted results.
#' @param plot logical; if \code{TRUE}, then this function plots
#'   predictive distribution of \eqn{\theta_{new}}.
#' @importFrom stats quantile sd density
#' @importFrom e1071 skewness kurtosis
#' @importFrom rstan extract
#' @examples
#' require("flexmeta")
#' @export
pred_f <- function(fit, plot = TRUE){

  e1 <- extract(fit)
  y <- e1$theta_new

  X1 <- mean(y)
  X2 <- sd(y)
  X3 <- quantile(y, c(.025, .975))
  X4 <- mean(y < 0)
  X5 <- quantile(y, c(.25, .50, .75))
  X6 <- skewness(y)
  X7 <- kurtosis(y)

  if(plot == TRUE){
    hist(y, xlab = "theta_new", freq = FALSE, br = 100,
         main = "Predictive Distribution")
    ds <- density(y)
    lines(ds$x, ds$y, col = 4, lwd = 2)
  }

  return(list(mean = X1, sd = X2, "95%PI" = X3, "Pr(theta_new<0)" = X4,
              quartiles = X5, skewness = X6, kurtosis = X7))

}


