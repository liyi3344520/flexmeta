#' An Object of initial values for the Sinh-arcsinh distribution
#'
#' @name sas_stan
#' @rdname sas_stan
#' @export
init_sas <- function() {
  list(sigma = runif(1, 0, 2), delta = runif(1, 0, 2),
       mu = runif(1, -2, 2), epsilon = runif(1, -2, 2),
       theta_new = runif(1, -2, 2))
}
