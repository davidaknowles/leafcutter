// Multinomial GLM (with no overdispersion)
data {
  int<lower=0> N;
  int<lower=0> P;
  int<lower=0> K;
  vector[P] x[N];
  int<lower=0> y[N,K];
}
parameters {
  simplex[K] beta_raw[P];
  real beta_scale[P]; 
}

model {
  matrix[K,P] beta;
  for (k in 1:K)
    for (p in 1:P)
      beta[k,p] <- beta_scale[p] * (beta_raw[p][k] - 1.0 / K);
  for (n in 1:N)
    y[n] ~ multinomial(softmax(beta * x[n]));
}