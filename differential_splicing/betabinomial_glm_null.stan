// DEPRECIATED. 
// Null beta binomial model, i.e. y ~ BB(n, conc*sigma(mu), conc*sigma(-mu))
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
  real mu;
}
model {
  real a;
  real b;
  real p; 
  p <- inv_logit(mu); 
  a <- conc*p;
  b <- conc*(1.0-p);
  // mu ~ normal(0,5);
  conc ~ gamma(concShape, concRate);
  ys ~ beta_binomial(ns, a, b);
}