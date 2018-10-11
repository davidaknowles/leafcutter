data {
  int<lower=0> N; // sample size
  int<lower=0> P; // number of covariates
  int<lower=0> K; // number of classes
  vector[P] x[N]; // covariates
  vector[K] y[N]; // counts 
  real<lower=0> concShape; // concentration shape
  real<lower=0> concRate; // concentration rate
  real<lower=0> outlier_prior_a; 
  real<lower=0> outlier_prior_b; 
}
parameters {
  simplex[K] beta_raw[P]; // symmetric K-1 dof encoding of beta
  real beta_scale[P]; 
  real<lower=0> conc[K]; // concentration parameter
  real<lower=0,upper=1> outlier_prob;
}

model {
  // beta reparameterization (Section 5.6 in the Stan manual)
  matrix[K,P] beta;
  for (k in 1:K)
    for (p in 1:P)
      beta[k,p] = beta_scale[p] * (beta_raw[p][k] - 1.0 / K);

  outlier_prob ~ beta(outlier_prior_a, outlier_prior_b);

  conc ~ gamma(concShape, concRate);
  for (n in 1:N) {
    vector[K] a; 
    real suma;
    real sumy; 
    vector[K] aPlusY;
    vector[K] lGaPlusY; 
    vector[K] lG1PlusY; 
    vector[K] lGaA ;
    vector[K] s; 
    s = softmax(beta * x[n]); 
    for (k in 1:K)
      a[k] = conc[k] * s[k]; 
    // explicit construction of multinomial dirichlet
    // y ~ multinomial_dirichlet( conc * softmax(beta * x[n]) )
    suma = sum(a);
    sumy = sum(y[n]);
    for (k in 1:K) {
      lGaPlusY[k] = lgamma(a[k]+y[n][k]);
      lGaA[k] = lgamma(a[k]);
      lG1PlusY[k] = lgamma(1.0+y[n][k]);
    }
    target += log_sum_exp(log(1.0-outlier_prob)+lgamma(suma)+sum(lGaPlusY)-lgamma(suma+sumy)-sum(lGaA), log(outlier_prob) + lgamma(K)+sum(lG1PlusY)-lgamma(K+sumy));
  }
}
