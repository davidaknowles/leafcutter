data {
  int<lower=0> N; 
  int<lower=0> P;
  matrix[N,P] x; 
  int<lower=0> ys[N];
  int<lower=0> ns[N];
  real<lower=0> concShape; 
  real<lower=0> concRate;  
}
parameters {
  real<lower=0> conc; 
  vector[P] beta;
}
model {
  vector[N] xb; 
  real a[N];
  real b[N];
  real p[N]; 
  xb = x * beta;
  for (n in 1:N) {
    p[n] = inv_logit(xb[n]); 
    a[n] = conc*p[n];
    b[n] = conc*(1.0-p[n]);
  }
  // beta ~ normal(0,5);
  conc ~ gamma(concShape, concRate);
  ys ~ beta_binomial(ns, a, b);
}
