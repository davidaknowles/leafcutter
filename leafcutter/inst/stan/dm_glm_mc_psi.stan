data {
  int<lower=0> N; // sample size
  int<lower=0> P; // number of covariates
  int<lower=0> K; // number of classes
  vector[P] x[N]; // covariates
  vector[K] y[N]; // counts 
  matrix[K,P] beta;
  real<lower=0> conc[K]; // concentration parameter
  real residual_sigma; 
}
parameters {
  vector[K] residual[N]; 
}
model {
  for (n in 1:N) {
    vector[K] a; 
    real suma;
    vector[K] aPlusY;
    vector[K] lGaPlusY; 
    vector[K] lGaA;
    vector[K] s; 
    residual[n] ~ normal(0,residual_sigma);
    s = softmax(beta * x[n] + residual[n]); 
    for (k in 1:K)
      a[k] = conc[k] * s[k]; 
    // explicit construction of multinomial dirichlet
    // y ~ multinomial_dirichlet( conc * softmax(beta * x[n]) )
    suma = sum(a);
    aPlusY = a + y[n];
    for (k in 1:K) {
      lGaPlusY[k] = lgamma(aPlusY[k]);
      lGaA[k] = lgamma(a[k]);
    }
    target += lgamma(suma)+sum(lGaPlusY)-lgamma(suma+sum(y[n]))-sum(lGaA); 
  }
}
