
data{
  int<lower=0> K;
  real yi[K];
  real<lower=0> si[K];
}
parameters {
  real xi;
  real<lower=0> omega;
  real alpha;
  real theta[K];
}
transformed parameters {
  real mu;
  real<lower=0> V;
  mu = xi + omega*sqrt(2/pi())*alpha/sqrt(1+alpha^2);
  V = omega^2 * (1 - (sqrt(2/pi())*alpha/sqrt(1+alpha^2))^2);
}
model{
  for(i in 1:K) {
    yi[i] ~ normal(theta[i], si[i]);
    theta[i] ~ skew_normal(xi, omega, alpha);
  }
  xi ~ normal(0, 100);
  omega ~ uniform(0, 20);
  alpha ~ normal(0, 5);
}
generated quantities{
  vector[K] log_lik;
  real theta_new;
  for (i in 1:K) log_lik[i] = normal_lpdf(yi[i] | theta[i], si[i]);
  theta_new = skew_normal_rng(xi, omega, alpha);
}
