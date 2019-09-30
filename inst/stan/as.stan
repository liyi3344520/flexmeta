
functions {
  int signnum(real x) {
    return x < 0 ? -1 : x > 0;
  }
  real as_lpdf(real y, real xi, real omega, real nu, real alpha) {
    return log(2) - log(omega) -
      log(2) - (1/nu)*log(nu) - log( tgamma(1 + 1/nu) ) -
      (fabs((y - xi)/omega) )^nu / nu +
      normal_lcdf( signnum(alpha*(y-xi)/omega) * fabs( alpha*(y-xi)/omega )^(nu/2) / sqrt(nu/2)  | 0, 1);
  }
}
data{
  int<lower=0> K;
  real yi[K];
  real<lower=0> si[K];
}
parameters {
  real xi;
  real<lower=0> omega;
  real alpha;
  real<lower=.5, upper=100> nu;
  real theta[K];
  real theta_new;
}
transformed parameters {
  real mu;
  real<lower=0> V;
  mu = xi + signnum(alpha)*omega*(nu^(1/nu))*tgamma(2/nu)*(2*student_t_cdf(sqrt((fabs(alpha)^nu)*4/nu),4/nu,0,1) - 1);
  V = omega^2 * ((nu^(2/nu))*tgamma(3/nu)/tgamma(1/nu) -
    ((nu^(1/nu))*tgamma(2/nu)*(2*student_t_cdf(sqrt((fabs(alpha)^nu )*4/nu),(4/nu),0,1) - 1))^2);
}
model{
  for(i in 1:K) {
    yi[i] ~ normal(theta[i], si[i]);
    theta[i] ~ as(xi, omega, nu, alpha);
  }
  xi ~ normal(0, 100);
  omega ~ uniform(0, 20);
  alpha ~ normal(0, 5);
  nu ~ exponential(0.10);
}
generated quantities{
  vector[K] log_lik;
  real theta_new;
  for (i in 1:K) log_lik[i] = normal_lpdf(yi[i] | theta[i], si[i]);
  theta_new = as(xi, omega, nu, alpha);
}
