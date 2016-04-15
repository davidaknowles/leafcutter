data {
  int<lower=0> N; 
  int<lower=0> P;
  matrix[N,P] x; 
  int<lower=0> y[N];
  real<lower=0> conc[N];
  real<lower=0> nb_size[N];
  vector<lower=0>[N] library_size;
}
parameters {
  vector[P] beta;
}
model {
  vector[N] exb;
  exb <- exp(x * beta); 
  for (n in 1:N) {
    real a; 
    real b; 
    real p;
    p <- nb_size[n] / (library_size[n] * exb[n] + nb_size[n]);
    a <- p * conc[n] + 1.0; 
    b <- (1.0-p) * conc[n] ;
    increment_log_prob(lgamma(nb_size[n]+y[n]) + lgamma(a+nb_size[n]) + lgamma(b+y[n]) + lgamma(a+b) - lgamma(nb_size[n]) - lgamma(a+b+y[n]+nb_size[n]) - lgamma(a) - lgamma(b));
  }
}
