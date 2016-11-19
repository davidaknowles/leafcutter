data {
  int<lower=0> D;
  int<lower=0> y[D];
  real<lower=0> r;
  simplex[N] library_size; 
}
parameters {
  simplex[N] library_size; 
  vector<lower=0>[D] intron_usage;  
  real<lower=0> r[N];
  real<lower=0> conc[D];
}
model {
  r ~ gamma(1.01,0.01); 
  conc ~ gamma(1.01,0.01); 
  for (d in 1:D) {
    for (n in 1:N) {
      real a; 
      real b; 
      real p;
      real mu; 
      mu <- library_size[n] * N * intron_usage[d];
      p <- r[n] / (mu + r[n]);
      a <- p * conc[d] + 1.0; 
      b <- (1.0-p) * conc[d] ;
      increment_log_prob(lgamma(r[n]+y[n,d]) + lgamma(a+r[n]) + lgamma(b+y[n,d]) + lgamma(a+b) - lgamma(r[n]) - lgamma(a+b+y[n,d]+r[n]) - lgamma(a) - lgamma(b) - lgamma(y[n,d]+1.0));
    }  
  }
}
