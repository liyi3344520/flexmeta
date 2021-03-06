% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/func-dic.r
\name{dic}
\alias{dic}
\title{Generate DIC}
\usage{
dic(fit, dat)
}
\arguments{
\item{fit}{An object of S4 class \code{stanfit} representing
the fitted results.}

\item{dat}{A list of the study data.
\itemize{
\item \code{y}: A vector of the effect size estimates.
\item \code{se}: A vector of the within studies standard error
  estimates.
}}
}
\value{
A numeric value of \eqn{DIC}.
}
\description{
Generate DIC
}
\examples{
\donttest{
require("flexmeta")
require("rstan")
options(mc.cores = parallel::detectCores())

data("clbp")
dat <- convert_mean(clbp$n1, clbp$m1, clbp$s1, clbp$n2, clbp$m2, clbp$s2)
dat <- list(K = length(dat$y), y = dat$y, se = dat$se)

fit <- sampling(sm_sn, data = dat,
  chains = 1, warmup = 5000, iter = 30000, seed = 2645751)

dic(fit, dat)
}
}
\seealso{
\code{\link[rstan]{sampling}}, \code{\link[rstan]{stanmodel-class}},
  \code{\link[rstan]{stanfit-class}}, \code{\link[rstan]{stan}}.
}
