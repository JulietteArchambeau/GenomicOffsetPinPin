data {
  int N;
  int NbTot[N];
  int NbSurv[N];
  vector[N] heightSC;
  vector[N] covariateSC;
}

parameters {
  real beta0;
  real betaCovariate;
}

model {
  NbSurv ~ binomial_logit(NbTot,beta0 + betaCovariate * covariateSC);
  
  beta0 ~ normal(0,5);
  betaCovariate ~ normal(0,5);
}

