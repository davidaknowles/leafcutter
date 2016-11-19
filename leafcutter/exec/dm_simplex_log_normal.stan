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
  int<lower=0> K; // number of introns
  vector[K] y[N]; // counts 
  vector[K] mu; // predicted log_dispersion
  real<lower=0> sigma; // sd in log_dispersion
}
parameters {
  simplex[K] proportion; 
  vector[K] log_dispersion; // log(1/concentration)
}
model {
  vector[K] conc;
  vector[K] a;
  log_dispersion ~ normal( mu , sigma);
  conc <-  exp(-log_dispersion);
  for (k in 1:K) {
    a[k] <- proportion[k] * conc[k]; 
  }
  for (n in 1:N) {
    y[n] ~ dirichlet_multinomial(a); 
  }
}
