data {
  int<lower=0> N; 
  int<lower=0> D;
  int<lower=0> y[N,D];
}
parameters {
  simplex[N] library_size; 
  vector<lower=0>[D] intron_usage;  
  real<lower=0> sqrtrbShape;
  real<lower=0> sqrtrbRate;
  real<lower=0> sqrtrOverBShape;
  real<lower=0> sqrtrOverBRate;
  real<lower=0> sqrtrb[D]; 
  real<lower=0> sqrtrOverB[D]; 
}
model {
  sqrtrbShape ~ gamma(1,1); 
  sqrtrbRate ~ gamma(1,100);
  sqrtrOverBShape ~ gamma(1,1); 
  sqrtrOverBRate ~ gamma(1,1); 
  sqrtrb ~ gamma( sqrtrbShape, sqrtrbRate );
  sqrtrOverB ~ gamma( sqrtrOverBShape, sqrtrOverBRate );
  for (d in 1:D) {
    real r;
    real b;
    r <- sqrtrb[d] * sqrtrOverB[d];
    b <- sqrtrb[d] / sqrtrOverB[d];
    for (n in 1:N) {
      real a; 
      real mu; 
      mu <- library_size[n] * N * intron_usage[d];
      a <- 1.0 + r * b / mu;
      increment_log_prob(lgamma(r+y[n,d]) + lgamma(a+r) + lgamma(b+y[n,d]) + lgamma(a+b) - lgamma(r) - lgamma(a+b+y[n,d]+r) - lgamma(a) - lgamma(b) - lgamma(y[n,d]+1));
    }  
  }
}
