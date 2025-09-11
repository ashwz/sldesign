
data {
  int<lower=0>   num_DL; // number of DL's having data
  int<lower=0>   y[num_DL];
  int<lower=0>   N_n[num_DL];
  vector[num_DL] x;
  real           mu_a;
  real           mu_b;
  real<lower=0>  sigma_a;
  real<lower=0>  sigma_b;
}

parameters {
  real alpha;
  real beta;
}

transformed parameters {
  vector<lower=0, upper=1>[num_DL] q;
  
  q = 1/(1+exp(- alpha - beta*x));
  
}

model {
  alpha ~ normal(mu_a, sigma_a);
  beta  ~ normal(mu_b, sigma_b);
  y     ~ binomial(N_n, q);
  
}


