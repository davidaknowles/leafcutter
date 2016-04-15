data {
  int<lower=0> N; 
  vector<lower=0>[N] library_size;
  int<lower=0> y[N];
  real<lower=0> sqrtrbShape;
  real<lower=0> sqrtrbRate;
  real<lower=0> sqrtrOverBShape;
  real<lower=0> sqrtrOverBRate;
}
parameters {
  real<lower=0> mu; 
  real<lower=0> sqrtrb; 
  real<lower=0> sqrtrOverB; 
}
model {
  real r;
  real b;
  r <- sqrtrb * sqrtrOverB;
  b <- sqrtrb / sqrtrOverB;
  sqrtrb ~ gamma( sqrtrbShape, sqrtrbRate );
  sqrtrOverB ~ gamma( sqrtrOverBShape, sqrtrOverBRate );
  for (n in 1:N) {
    real a; 
    a <- 1.0 + r * b / ( library_size[n] * mu ); 
    increment_log_prob(lgamma(r+y[n]) + lgamma(a+r) + lgamma(b+y[n]) + lgamma(a+b) - lgamma(r) - lgamma(a+b+y[n]+r) - lgamma(a) - lgamma(b) - lgamma(y[n]+1));
  }
}
