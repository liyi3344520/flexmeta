
data{
  int<lower=0> K;
  real yi[K];
  real<lower=0> si[K];
}
parameters {
  real xi;
  real<lower=0> omega;
  real<lower=2.5, upper=1000> nu;
  real theta[K];
}
transformed parameters {
  real mu;
  real<lower=0> V;
  mu = xi;
  V = omega*omega;
}
model{
  for(i in 1:K) {
    yi[i] ~ normal(theta[i], si[i]);
    theta[i] ~ student_t(nu, xi, omega);
  }
  xi ~ normal(0, 100);
  omega ~ uniform(0, 20);
  nu ~ exponential(0.10);
}
generated quantities{
  vector[K] log_lik;
  real theta_new;
  for (i in 1:K) log_lik[i] = normal_lpdf(yi[i] | theta[i], si[i]);
  theta_new = student_t_rng(nu, xi, omega);
}
