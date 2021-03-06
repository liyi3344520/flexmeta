% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zzz-stanmodels.r
\docType{data}
\name{sm_norm}
\alias{sm_norm}
\alias{stanmodels}
\alias{sm}
\alias{sm_t}
\alias{sm_sn}
\alias{sm_st}
\alias{sm_as}
\alias{sm_jf}
\alias{sm_sas}
\title{Compiled \code{stanmodel} objects of flexible random-effects distribution models}
\format{An object of class \code{stanmodel}.}
\usage{
sm_norm

sm_t

sm_sn

sm_st

sm_as

sm_jf

sm_sas
}
\description{
Compiled \code{stanmodel} objects of flexible random-effects distribution models
}
\details{
Compiled \code{stanmodel} objects of flexible random
effects distribution models that are defined as
\deqn{Y_i \sim N(\theta_i, \sigma_i^2)}
\deqn{\theta_i \sim F}
for \eqn{i=1,\ldots,K}, where \eqn{\theta_i} is the true effect
size of the ith study, \eqn{\sigma_i^2} is within-studies
variances that is usually assumed to be known and fixed to their
valid estimates, and \eqn{F} corresponds to the random-effects
distribution that expresses the probability distribution how
\eqn{\theta_i} is heterogeneity distributed.
The \code{flexmeta} package supports seven random-effects
distributions, \eqn{F}, as follows:
\itemize{
\item \code{sm_norm}: A \code{stanmodel} object for a normal random
  effects distribution, \eqn{F \sim N(\xi, \omega)},
  where \eqn{\xi} is location (mean), and \eqn{\omega} is scale
  parameter (SD).
\item \code{sm_t}: A \code{stanmodel} object for a non-central
  student-\eqn{t} random-effects distribution,
  \eqn{F \sim t(\nu, \mu, \omega)}, where \eqn{\xi} is location (mean),
  \eqn{\omega} is scale parameter (SD), and \eqn{\nu} is d.f.
\item \code{sm_sn}: A \code{stanmodel} object for a skew normal random
  effects distribution, \eqn{F \sim skew-N(\xi, \omega, \alpha)},
  where \eqn{\xi} is location (mean), \eqn{\omega} is scale parameter (SD),
  and \eqn{\alpha} is shape parameter.
\item \code{sm_t}: A \code{stanmodel} object for a skew \eqn{t} random
  effects distribution, \eqn{F \sim skew-t(\mu, \omega, \nu, \alpha)},
  where \eqn{\xi} is location (mean), \eqn{\omega} is scale parameter (SD),
  and \eqn{\alpha} is shape parameter.
\item \code{sm_as}: A \code{stanmodel} object for an asymmetric Subboton
  (Type II) random-effects distribution,
  \eqn{F \sim AS2(\mu, \omega, \nu, \alpha)},
  where \eqn{\xi} is location (mean), \eqn{\omega} is scale parameter (SD),
  and \eqn{\alpha} is shape parameter.
\item \code{sm_jf}: A \code{stanmodel} object for an Jones-Faddy
  random-effects distribution, \eqn{F \sim JF(\mu, \omega, a, b)},
  where \eqn{\xi} is location (mean), \eqn{\omega} is scale parameter (SD),
  and \eqn{a} and \eqn{b} are distribution parameters regulate
  the skewness and kurtosis.
\item \code{sm_sas}: A \code{stanmodel} object for an Sinh-arcsinh
  random-effects distribution, \eqn{F \sim SAS(\mu, \omega, \delta, \epsilon)},
  where \eqn{\xi} is location (mean), \eqn{\omega} is scale parameter (SD),
  \eqn{\delta} is kurtosis parameter and \eqn{\epsilon} is skewness parameter.
}
The \code{data} should be a list of study data:
\itemize{
\item \code{K}: A numeric value of the number of studies.
\item \code{y}: A numeric vector of the effect size estimates.
\item \code{se}: A numeric vector of the within studies standard error estimates.
}
}
\examples{
\donttest{
require("flexmeta")
require("rstan")
options(mc.cores = parallel::detectCores())

data("clbp")
dat <- convert_mean(clbp$n1, clbp$m1, clbp$s1, clbp$n2, clbp$m2, clbp$s2)
dat <- list(K = length(dat$y), y = dat$y, se = dat$se)

fit1 <- sampling(sm_norm, data = dat,
  chains = 1, warmup = 5000, iter = 30000, seed = 2645751)
fit2 <- sampling(sm_t, data = dat,
  chains = 1, warmup = 5000, iter = 30000, seed = 2645751)
fit3 <- sampling(sm_sn, data = dat,
  chains = 1, warmup = 5000, iter = 30000, seed = 2645751)
fit4 <- sampling(sm_st, data = dat,
  chains = 1, warmup = 5000, iter = 30000, seed = 2645751)
fit5 <- sampling(sm_as, data = dat,
  chains = 1, warmup = 5000, iter = 30000, seed = 2645751)
fit6 <- sampling(sm_jf, data = dat,
  chains = 1, warmup = 5000, iter = 30000, seed = 2645751)
fit7 <- sampling(sm_sas, data = dat, init = init_sas,
  chains = 1, warmup = 5000, iter = 30000, seed = 2645751)

e1 <- extract(fit1)
e2 <- extract(fit2)
e3 <- extract(fit3)
e4 <- extract(fit4)
e5 <- extract(fit5)
e6 <- extract(fit6)
e7 <- extract(fit7)

par(mfrow = c(2, 4))
hist(e1$theta_new, br = 150, xlim = c(-40, 30), ylim = c(0, .1),
  freq = FALSE, main = "(a) Normal distribution",
  xlab = "Mean difference", border = "gray")
ds <- density(e1$theta_new)
lines(ds$x, ds$y, col = 4, lwd = 2)

hist(e2$theta_new, br = 900, xlim = c(-40, 30), ylim = c(0, .1),
  freq = FALSE, main = "(b) t-distribution",
  xlab = "Mean difference", border = "gray")
ds <- density(e2$theta_new)
lines(ds$x, ds$y, col = 4, lwd = 2)

hist(e3$theta_new, br = 100, xlim = c(-40, 30), ylim = c(0, .1),
  freq = FALSE, main = "(c) Skew normal distribution",
  xlab = "Mean difference", border = "gray")
ds <- density(e3$theta_new)
lines(ds$x, ds$y, col = 4, lwd = 2)

hist(e4$theta_new, br = 400, xlim = c(-40, 30), ylim = c(0, .1),
  freq = FALSE, main = "(d) Skew t-distribution",
  xlab = "Mean difference", border = "gray")
ds <- density(e4$theta_new)
lines(ds$x, ds$y, col = 4, lwd = 2)

hist(e5$theta_new, br = 200, xlim = c(-40, 30), ylim = c(0, .1),
  freq = FALSE, main = "(e) Asymmetric Subbotin distribution",
  xlab = "Mean difference", border = "gray")
ds <- density(e5$theta_new)
lines(ds$x, ds$y, col = 4, lwd = 2)

hist(e6$theta_new, br = 400, xlim = c(-40, 30), ylim = c(0, .1),
  freq = FALSE, main = "(f) Jones-Faddy distribution",
  xlab = "Mean difference", border = "gray")
ds <- density(e6$theta_new)
lines(ds$x, ds$y, col = 4, lwd = 2)

hist(e7$theta_new, br = 300, xlim = c(-40, 30), ylim = c(0, .1),
  freq = FALSE, main = "(g) Sinh-arcsinh distribution",
  xlab = "Mean difference", border = "gray")
ds <- density(e7$theta_new)
lines(ds$x, ds$y, col = 4, lwd = 2)

pred_f(fit1)
pred_f(fit2)
pred_f(fit3)
pred_f(fit4)
pred_f(fit5)
pred_f(fit6)
pred_f(fit7)

mu_f(fit1)
mu_f(fit2)
mu_f(fit3)
mu_f(fit4)
mu_f(fit5)
mu_f(fit6)
mu_f(fit7)

dic(fit1, dat)
dic(fit2, dat)
dic(fit3, dat)
dic(fit4, dat)
dic(fit5, dat)
dic(fit6, dat)
dic(fit7, dat)
}
}
\seealso{
\code{\link[rstan]{sampling}}, \code{\link[rstan]{stanmodel-class}},
  \code{\link[rstan]{stanfit-class}}, \code{\link[rstan]{stan}}.
}
\keyword{datasets}
