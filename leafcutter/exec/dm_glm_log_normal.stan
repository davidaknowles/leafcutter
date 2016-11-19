functions {
  real dirichlet_multinomial_log(vector y, vector a) {
    real suma;
    real lp;
    int K; 
    vector[rows(y)] aPlusY;
    vector[rows(y)] lGaPlusY; 
    vector[rows(y)] lGaA ;
    K <- rows(y);
    suma <- sum(a);
    aPlusY <- a + y;
    for (k in 1:K) {
      lGaPlusY[k] <- lgamma(aPlusY[k]);
      lGaA[k] <- lgamma(a[k]);
    }
    return lgamma(suma)+sum(lGaPlusY)-lgamma(suma+sum(y))-sum(lGaA);
  }
}
data {
  int<lower=0> N; // sample size
  int<lower=0> P; // number of covariates
  int<lower=0> K; // number of introns
  vector[P] x[N]; // covariates
  vector[K] y[N]; // counts 
  vector[K] mu; // predicted log_dispersion
  real<lower=0> sigma; // sd in log_dispersion
}
parameters {
  simplex[K] beta_raw[P]; // symmetric K-1 dof encoding of beta
  real beta_scale[P]; 
  real log_dispersion[K]; // log(1/concentration)
}
model {
  // beta reparameterization (Section 5.6 in the Stan manual)
  matrix[K,P] beta;
  for (k in 1:K)
    for (p in 1:P)
      beta[k,p] <- beta_scale[p] * (beta_raw[p][k] - 1.0 / K);
  log_dispersion ~ normal( mu , sigma);
  for (n in 1:N) {
    vector[K] a; 
    vector[K] s; 
    s <- softmax(beta * x[n]); 
    for (k in 1:K)
      a[k] <- exp(-log_dispersion[k]) * s[k]; 
    y[n] ~ dirichlet_multinomial(a); 
  }
}
