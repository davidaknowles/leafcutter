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
  real<lower=0> conc;
}
model {
  conc ~ gamma( concShape, concRate );
  y ~ neg_binomial_2_log( x * beta, conc ); 
}
