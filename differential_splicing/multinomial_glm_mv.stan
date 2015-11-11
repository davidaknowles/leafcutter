// Multinomial GLM (no overdispersion, and too many degrees of freedom)

data {
  int<lower=0> N;
  int<lower=0> P;
  int<lower=0> K;
  vector[P] x[N];
  int<lower=0> y[N,K];
  //real<lower=0> concShape;
  //real<lower=0> concRate;
}
parameters {
  // real<lower=0> conc; 
  matrix[K,P] beta;
}

model {
  to_vector( beta ) ~ normal(0,5);
  // conc ~ gamma(concShape, concRate);
  for (n in 1:N)
    y[n] ~ multinomial(softmax(beta * x[n]));
}