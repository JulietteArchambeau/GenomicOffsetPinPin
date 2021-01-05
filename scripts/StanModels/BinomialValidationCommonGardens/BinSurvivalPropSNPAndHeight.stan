data {
  int N;
  int NbTot[N];
  int NbSurv[N];
  vector[N] heightSC;
  vector[N] propsnpSC;
}

parameters {
  real beta0;
  real betaHeight;
  real betaPropSNP;
}

model {
  NbSurv ~ binomial_logit(NbTot,beta0 + betaPropSNP * propsnpSC + betaHeight * heightSC);
  
  beta0 ~ normal(0,5);
  betaHeight ~ normal(0,5);
  betaPropSNP ~ normal(0,5);
}

