
functions {
  real besselK(real v, real z);
  real sinh_arcsin_lpdf(real y, real mu, real sigma, real epsilon, real delta) {
    return normal_lpdf(sinh(delta*asinh((y - mu)/sigma) - epsilon) | 0, 1) +
      log(delta) - 0.5*log(1 + pow(y - mu, 2)/pow(sigma, 2)) - log(sigma) +
      log(cosh(delta*asinh((y - mu)/sigma) - epsilon));
  }
  real sinh_arcsin_rng(real mu, real sigma, real epsilon, real delta) {
    return mu + sigma*sinh((asinh(normal_rng(0, 1)) + epsilon)/delta);
  }
}
data {
  int<lower=0> K;
  real yi[K];
  real<lower=0> si[K];
}
parameters {
  real mu;
  real<lower=0> sigma;
	real epsilon;
  real<lower=0, upper=1000> delta;
  real theta[K];
}
model {
  for(i in 1:K) {
    yi[i] ~ normal(theta[i], si[i]);
  }
  for(i in 1:K) {
    theta[i] ~ sinh_arcsin(mu, sigma, epsilon, delta);
  }
  mu ~ normal(0, 100);
  sigma ~ uniform(0, 20);
  epsilon ~ normal(0, 4);
  delta ~ uniform(0, 20);
}
generated quantities {
  real sx_post_mean;
  real sx_post_m2;
  real sx_post_var;
  real x_post_mean;
  real x_post_var;
  real theta_new;
  sx_post_mean = sinh(epsilon/delta)*exp(0.25)/sqrt(8.0*pi())*
      (besselK((1.0/delta + 1.0)*0.5, 0.25) + besselK((1.0/delta - 1.0)*0.5, 0.25));
  sx_post_m2 = 0.5*(cosh(2.0*epsilon/delta)*exp(0.25)/sqrt(8.0*pi())*
      (besselK((2.0/delta + 1.0)*0.5, 0.25) + besselK((2.0/delta - 1.0)*0.5, 0.25)) - 1.0);
  sx_post_var = sx_post_m2 - sx_post_mean^2;
  x_post_mean = sigma*sx_post_mean + mu;
  x_post_var = sigma^2*sx_post_var;
  theta_new = sinh_arcsin_rng(mu, sigma, epsilon, delta);
}
