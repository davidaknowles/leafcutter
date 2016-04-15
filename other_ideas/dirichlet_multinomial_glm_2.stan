data {
  int<lower=0> N; // sample size
  int<lower=0> P; // number of covariates
  int<lower=0> K; // number of classes
  vector[P] x[N]; // covariates
  vector[K] y[N]; // counts 
}
parameters {
  matrix[K,P] beta; 
}

model {
  
  for (n in 1:N) {
    vector[K] a; 
    real suma;
    vector[K] aPlusY;
    vector[K] lGaPlusY; 
    vector[K] lGaA ;
    a <- exp( beta * x[n] ); 
    // explicit construction of multinomial dirichlet
    // y ~ multinomial_dirichlet( exp(beta * x[n]) )
    suma <- sum(a);
    aPlusY <- a + y[n];
    for (k in 1:K) {
      lGaPlusY[k] <- lgamma(aPlusY[k]);
      lGaA[k] <- lgamma(a[k]);
    }
    increment_log_prob(lgamma(suma)+sum(lGaPlusY)-lgamma(suma+sum(y[n]))-sum(lGaA));
  }
}