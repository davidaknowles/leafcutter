data {
  int<lower=0> N; 
  int<lower=0> P;
  matrix[N,P] x; 
  int<lower=0> y[N];
  real<lower=0> r[N];
  real<lower=0> b[N];
}
parameters {
  vector[P] beta;
}
model {
  vector[N] exb;
  exb <- exp(x * beta); 
  for (n in 1:N) {
    real a; 
    a <- 1.0 + r[n] * b[n] / exb[n];
    // excluding -lgamma(y[n]+1)
    increment_log_prob(lgamma(r[n]+y[n]) + lgamma(a+r[n]) + lgamma(b[n]+y[n]) + lgamma(a+b[n]) - lgamma(r[n]) - lgamma(a+b[n]+y[n]+r[n]) - lgamma(a) - lgamma(b[n]));
  }
}
