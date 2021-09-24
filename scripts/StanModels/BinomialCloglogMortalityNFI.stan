data {
  int N;
  vector[N] logNbYears;                                                           // Offset to account for different census intervals
  vector[N] C;                                                                    // Proxy of the competition among trees (i.e. basal area)
  vector[N] GO;                                                                   // Genomic offset
  int NbDead[N];                                                                  // Number of dead trees in the plot    
  int NbTot[N];                                                                   // Total number of trees in the plot
  int<lower=0> ncountry;                                                          // Number of countries
  int<lower=0, upper=ncountry> country[N];                                        // Countries
}
parameters {
  vector[ncountry] alpha_country;
  vector[ncountry] betaC_country;
  vector[ncountry] betaGO_country;
}

model {
  NbDead ~ binomial(NbTot,inv_cloglog(alpha_country[country] + betaC_country[country] .* C  + betaGO_country[country] .* GO + logNbYears));
  
  alpha_country ~ normal(0, 1);
  betaC_country ~ normal(0, 1);
  betaGO_country ~ normal(0, 1);
}
