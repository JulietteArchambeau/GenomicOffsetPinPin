data {
  int N;
  vector[N] H;      // Response variable (mean height of the populations)
  vector[N] X;      // Genomic offset or climatic transfer distance of the population
}

parameters {
  real beta0;
  real betaX1;
  real betaX2;
  real<lower = 0>  sigma_r;
}

transformed parameters {
  vector[N] mu;    // linear predictor
  real R_squared;  // R^2 to evaluate the goodness of fit of the model
  
  mu = beta0 + betaX1 * X + betaX2 * square(X);
  R_squared = 1 - variance(H - mu) / variance(H);
}

model {
  
  H ~ normal(mu, sigma_r); // Likelihood
  
  sigma_r ~ exponential(1);
  beta0 ~ normal(0,1);//std_normal();
  betaX1 ~ normal(0,1);//std_normal();
  betaX2 ~ normal(0,1);//std_normal();
}

generated quantities{
  // log likelihood for loo
  vector[N] log_lik;
  for (n in 1:N) {
    log_lik[n] = normal_lpdf(H[n] |mu[n],sigma_r);
  }
}

