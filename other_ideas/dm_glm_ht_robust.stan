data {
  int<lower=0> N; // sample size
  int<lower=0> P; // number of covariates
  int<lower=0> K; // number of classes
  int<lower=0> M; // number of mixtures
  vector[P] x[N]; // covariates
  vector[K] y[N]; // counts 
  real<lower=0> concShape; // concentration shape
  real<lower=0> concRate; // concentration rate
  real<lower=0> outlier_prior_a; 
  real<lower=0> outlier_prior_b; 
  //real<lower=0>[M] alpha; // prior counts for theta, usually 1/M
}
parameters {
  simplex[K] beta_raw[P]; // symmetric K-1 dof encoding of beta
  real beta_scale[P]; 
  real<lower=0> conc[M,K]; // concentration parameter
  simplex[M] theta; // mixing proportions
  real<lower=0,upper=1> outlier_prob;
}

model {
  matrix[K,P] beta;
  for (k in 1:K)
    for (p in 1:P)
      beta[k,p] <- beta_scale[p] * (beta_raw[p][k] - 1.0 / K);

  theta ~ dirichlet(rep_vector(1.0/M, M));

  for (m in 1:M)
    for (k in 1:K)
      conc[m,k] ~ gamma(concShape, concRate);
      
  for (n in 1:N) {
    vector[K] s; 
    real ps[M];
    real sumy; 
    vector[K] lG1PlusY; 
    for (k in 1:K)
      lG1PlusY[k] <- lgamma(1.0+y[n][k]);
    sumy <- sum(y[n]);
    s <- softmax(beta * x[n]);
    for (m in 1:M) {    
      vector[K] a; 
      real suma;
      vector[K] aPlusY;
      vector[K] lGaPlusY; 
      vector[K] lGaA ;
      for (k in 1:K)
        a[k] <- conc[m,k] * s[k]; 
      suma <- sum(a);
      aPlusY <- a + y[n];
      for (k in 1:K) {
        lGaPlusY[k] <- lgamma(aPlusY[k]);
        lGaA[k] <- lgamma(a[k]);
      }
      ps[m] <- log(theta[m])+lgamma(suma)+sum(lGaPlusY)-lgamma(suma+sumy)-sum(lGaA);
    }
    increment_log_prob(log_sum_exp(log(1.0-outlier_prob)+log_sum_exp(ps), log(outlier_prob) + lgamma(K)+sum(lG1PlusY)-lgamma(K+sumy))); 
  }
}
