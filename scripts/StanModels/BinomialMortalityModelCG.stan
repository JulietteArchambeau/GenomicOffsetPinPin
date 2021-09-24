data {
  int N;
  int NbTot[N];     // Total number of trees in the population
  int NbDead[N];    // Number of dead trees in the population
  vector[N] H;      // Mean tree height of the population
  vector[N] X;      // Genomic offset or climatic transfer distance of the population
}

parameters {
  real beta0;
  real betaH;
  real betaX1;
  //real betaX2;
}

model {
  //NbDead ~ binomial_logit(NbTot,beta0 + betaH * H + betaX1 * X + betaX2 * square(X));
  NbDead ~ binomial_logit(NbTot,beta0 + betaH * H + betaX1 * X);
  
  beta0 ~ normal(0,5);//std_normal();
  betaH ~ normal(0,5);//std_normal();
  betaX1 ~ normal(0,5);//std_normal();
  //betaX2 ~ normal(0,5);//std_normal();
}

generated quantities{
  vector[N] log_lik;
  // log likelihood for loo
  for (n in 1:N) {
    //log_lik[n] = binomial_logit_lpmf( NbDead[n] | NbTot[n] , beta0 + betaH * H[n] + betaX1 * X[n] + betaX2 * square(X[n]));
    log_lik[n] = binomial_logit_lpmf( NbDead[n] | NbTot[n] , beta0 + betaH * H[n] + betaX1 * X[n]);
  }
}
