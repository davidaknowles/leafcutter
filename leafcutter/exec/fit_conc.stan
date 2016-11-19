data {
  int<lower=0> N; 
  real y[N];
  real mu[N]; 
}
parameters {
  real<lower=0> a;
  real<lower=0> b;
}
model {
  for (n in 1:N)
    y[n] ~ normal( log( b + a/mu[n]  ), 1 );  
}
