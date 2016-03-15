data {
  int<lower=0> N; 
  int<lower=0> P;
  matrix[N,P] x; 
  int<lower=0> y[N];
  real<lower=0> concShape;
  real<lower=0> concRate;
}
parameters {
  vector[P] beta;
  real<lower=0> r;
  real<lower=0> b;
}
model {
  vector[N] exb;
  exb <- exp(x * beta); 
  r ~ gamma( concShape, concRate );
  b ~ gamma( concShape, concRate );
  for (n in 1:N) {
    real a; 
    a <- 1.0 + r * b / exb[n];
    // excluding -lgamma(y[n]+1)
    increment_log_prob(lgamma(r+y[n]) + lgamma(a+r) + lgamma(b+y[n]) + lgamma(a+b) - lgamma(r) - lgamma(a+b+y[n]+r) - lgamma(a) - lgamma(b));
  }
}
