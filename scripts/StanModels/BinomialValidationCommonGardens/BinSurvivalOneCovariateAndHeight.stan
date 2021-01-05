data {
  int N;
  int NbTot[N];
  int NbSurv[N];
  vector[N] heightSC;
  vector[N] covariateSC;
}

parameters {
  real beta0;
  real betaHeight;
  real betaCovariate;
}

model {
  NbSurv ~ binomial_logit(NbTot,beta0 + betaCovariate * covariateSC + betaHeight * heightSC);
  
  beta0 ~ normal(0,5);
  betaHeight ~ normal(0,5);
  betaCovariate ~ normal(0,5);
}

