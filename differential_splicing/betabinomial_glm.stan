// DEPRECIATED. 
// Alternative UNIVARIATE beta binomial model, i.e. y ~ BB(n, conc*sigma(beta*x+mu), conc*sigma(-beta*x-mu))
data {
  int<lower=0> N; 
  int<lower=0> ys[N];
  int<lower=0> ns[N];
  real x[N];
  real<lower=0> concShape; 
  real<lower=0> concRate;  
}
parameters {
  real<lower=0> conc; 
  real beta;
  real mu;
}
model {
  real a[N];
  real b[N];
  real p[N]; 
  for (n in 1:N) {
    p[n] <- inv_logit(beta*x[n]+mu); 
    a[n] <- conc*p[n];
    b[n] <- conc*(1.0-p[n]);
  }
  // beta ~ normal(0,5);
  // mu ~ normal(0,5);
  conc ~ gamma(concShape, concRate);
  ys ~ beta_binomial(ns, a, b);
}