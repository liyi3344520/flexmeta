
functions {
  real tbeta(real a, real b) {
    return tgamma(a)*tgamma(b)/tgamma(a+b);
  }
  real jf_lpdf(real y, real xi, real omega, real a, real b) {
    return -log(omega) -(a + b - 1)*log(2) - 0.5*log(a + b) - log(tbeta(a, b)) +
      (a + 0.5)*log(1 + (y - xi)/(omega*sqrt(a + b + (y - xi)*(y - xi)/(omega*omega)))) +
      (b + 0.5)*log(1 - (y - xi)/(omega*sqrt(a + b + (y - xi)*(y - xi)/(omega*omega))));
  }
  real jf_rng(real xi, real omega, real a, real b) {
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
  real<lower=1.5> a;
  real<lower=1.5> b;
  real theta[K];
}
transformed parameters {
  real mu;
  real<lower=0> V;
  mu = xi + omega*0.5*(a-b)*sqrt(a+b)*tgamma(a-.5)*tgamma(b-0.5)/(tgamma(a)*tgamma(b));
  V = omega^2 * (((a+b)/4)*(( (a-b)^2 + a + b - 2 )/((a-1)*(b-1))) -
    (0.5*(a-b)*sqrt(a+b)*tgamma(a-0.5)*tgamma(b-0.5)/(tgamma(a)*tgamma(b)))^2);
}
model{
  for(i in 1:K) {
    yi[i] ~ normal(theta[i], si[i]);
    theta[i] ~ jf(xi, omega, a, b);
  }
  xi ~ normal(0, 100);
  omega ~ uniform(0, 20);
  a ~ uniform(0, 200);
  b ~ uniform(0, 200);
}
generated quantities{
  vector[K] log_lik;
  real theta_new;
  for (i in 1:K) log_lik[i] = normal_lpdf(yi[i] | theta[i], si[i]);
  theta_new = jf_rng(xi, omega, a, b);
}
