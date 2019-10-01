
functions {
  real skew_t_lpdf(real y, real xi, real omega, real nu, real alpha) {
    return log(2) - log(omega) + student_t_lpdf( (y - xi)/omega | nu, 0, 1) +
      student_t_lcdf(alpha*(y-xi)*sqrt((nu+1)/(nu+(y-xi)*(y-xi)/(omega*omega)))/omega | (nu+1), 0, 1);
  }
  real skew_t_rng(real xi, real omega, real nu, real alpha) {
    return 1;
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
  real<lower=2.5, upper=1000> nu;
  real theta[K];
}
transformed parameters {
  real mu;
  real<lower=0> V;
  mu = xi + omega*sqrt(nu)*tgamma(0.5*(nu-1))*alpha/(sqrt(pi())*tgamma(0.5*nu)*sqrt(1+alpha^2));
  V = omega^2 * (nu/(nu-2) - (sqrt(nu)*tgamma(0.5*(nu-1))*alpha/(sqrt(pi())*tgamma(0.5*nu)*sqrt(1+alpha^2)))^2);
}
model{
  for(i in 1:K) {
    yi[i] ~ normal(theta[i], si[i]);
    theta[i] ~ skew_t(xi, omega, nu, alpha);
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
  theta_new = skew_t_rng(xi, omega, nu, alpha);
}
