#' Generate posterior distribution of \eqn{V}
#'
#' @name V_f
#' @rdname V_f
#' @param fit An object of S4 class \code{stanfit} representing
#'   the fitted results.
#' @param plot logical; if \code{TRUE}, then this function plots
#'   posterior distribution of \eqn{V}.
#' @importFrom stats quantile sd density
#' @importFrom e1071 skewness kurtosis
#' @importFrom rstan extract
#' @examples
#' require("flexmeta")
#' @export
V_f <- function(fit, plot = TRUE){

  e1 <- extract(fit)
  v <- e1$V

  X1 <- mean(v)
  X2 <- sd(v)
  X3 <- quantile(v, c(.025, .975))
  X4 <- quantile(v, c(.25, .50, .75))
  X5 <- skewness(v)
  X6 <- kurtosis(v)

  if(plot == TRUE){
    hist(v, xlab = "V", freq = FALSE, br = 100,
         main = "Posterior Distribution of V")
    ds <- density(v)
    lines(ds$x,ds$y,col=4,lwd=2)
  }

  return(list(mean=X1,sd=X2,"95%CrI"=X3,quartiles=X4,skewness=X5,kurtosis=X6))

}
