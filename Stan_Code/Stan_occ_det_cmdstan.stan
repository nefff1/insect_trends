functions {
  // used to parallelise part of likelihood estimation
  
  real partial_sum_psi(int[] y_slice,
                   int start, int end,
                   vector logit_) {
    return bernoulli_logit_lpmf(1 | logit_[start:end]);
  }
  
  real partial_sum_p(int[] y_slice,
                   int start, int end,
                   vector logit_) {
    return bernoulli_logit_lpmf(y_slice | logit_[start:end]);
  }
}


data {
  int<lower=1> P; // number of unique site x year combinations
  int<lower=0> V[P]; // Number of visits per site x year combination (vector)
  int<lower=1> N; // Total number of visits
  int<lower=1> OCC; // number of occupied sites
  int<lower=1> OCC_N; // number of occupied sites
  int<lower=1> NOCC; // number of not occupied sites (but visited)

  int<lower=0> K_fo; // Number of fixed  variables occurrence probability
  int<lower=0> K_ro; // Number of random  variables occurrence probability
  int<lower=0> K_fd; // Number of fixed  variables detection probability
  int<lower=0> K_rd; // Number of random  variables detection probability
  
  int<lower=1> L_year;  // Number of years
  int<lower=1> L_area; // number of areas for which separate random walks are estimated
  int<lower=1> L_ro[K_ro]; // Number of levels of random variables at occurrence probability level
  int<lower=1> L_rd[K_rd]; // Number of levels of random variables at detection probability level

  int x_year[P, 1]; // year of observation
  int x_area[P, 1]; // area of observation
  matrix[P, K_fo] x_fo; // fixed variables occurrence probability
  int x_ro[P, K_ro]; // random variables occurrence probability 
  matrix[N, K_fd] x_fd; // fixed variables detection probability 
  int x_rd[N, K_rd]; // random variables detection probability
  
  int<lower=0> y[N];   // Observation

  int<lower=1> occ[OCC]; // vector of occupied siteyears
  int<lower=1> not_occ[NOCC]; // vector not occupied siteyears
  int<lower=1> occ_visit[OCC_N]; // vector of occupied siteyears
  
  int<lower=0> start_i[P]; // start of vector entries per siteyear
  int<lower=0> end_i[P]; // end of vector entries per siteyear
}

transformed data {
  int<lower=0> max_L_ro; // maximum number of levels for random effects occupancy
  int<lower=0> max_L_rd; // maximum number of levels for random effects detecability

   if (K_ro > 0){
     max_L_ro = max(L_ro);
   } else {
     max_L_ro = 0;
   }
   
  if (K_rd > 0){
     max_L_rd = max(L_rd);
   } else {
     max_L_rd = 0;
   }

}

parameters {
  // global intercept occurrence probability
  real mu_o; // Variablengrössenbeschränkung für schnellere Modelle

  // year effect occurrence probability
  matrix[L_year, L_area] alpha_year;
  vector<lower=0,upper=100>[L_area] sigma_a_year;

  // random effects occurrence probability
  vector[max_L_ro] alpha_ro[K_ro];
  vector<lower=0,upper=100>[K_ro] sigma_a_ro;
  
  // fixed effects occurrence probability
  vector<lower=-10,upper=10>[K_fo] beta_fo; // Beschränkung durch Prior gegeben...
  
  // global intercept detection probability
  real mu_d; // Variablengrössenbeschränkung für schnellere Modelle

  // random effects detection probability
  vector[max_L_rd] alpha_rd[K_rd];
  vector<lower=0,upper=100>[K_rd]  sigma_a_rd;

  // fixed effects detection probability
  vector<lower=-10,upper=10>[K_fd] beta_fd; // Beschränkung durch Prior gegeben...
}



transformed parameters {
  vector[P] logit_psi;  // Logit occupancy probability
  vector[N] logit_p_v; // Logit detection probability in Vektor-Form

  for (i in 1:P){
    logit_psi[i] = mu_o + alpha_year[x_year[i, 1], x_area[i, 1]];
  }
   
  if (K_fo > 0) // check whether there are actually any fixed effects in the occupancy model
  logit_psi += x_fo * beta_fo;

  if (K_ro > 0){ // check whether there are actually any random effects in the occupancy model
      for (i in 1:K_ro){
      logit_psi += alpha_ro[i][x_ro[, i]];
    }
  }


  logit_p_v = rep_vector(mu_d, N);
    
  if (K_fd > 0) { // check whether there are actually any covariates in the detection probability model
    logit_p_v += x_fd * beta_fd;
  } 
  
  if (K_rd > 0){ // check whether there are actually any random effects in the occupancy model
      for (i in 1:K_rd){
      logit_p_v += alpha_rd[i][x_rd[, i]];
    }
  }  
}

model {
  // grainsize for paralleling functions
  int grainsize1 = 1;
  int grainsize2 = 1;
  
  // Priors --------------------------------------------------------------------.
  target += normal_lpdf(mu_o | 0, 1.5);
  
  target += cauchy_lpdf(sigma_a_year | 0, 1); 
  
  for (j in 1:L_area){
      target += normal_lpdf(alpha_year[1,j] | 0, 1.5);
  }
  
  for (i in 2:L_year){
    for (j in 1:L_area){
        target += normal_lpdf(alpha_year[i, j] | alpha_year[i-1, j], sigma_a_year[j]); 
    }
  }

  target += cauchy_lpdf(sigma_a_ro | 0, 1);
  
  for (i in 1:K_ro){
    target += normal_lpdf(alpha_ro[i][1:L_ro[i]] | 0, sigma_a_ro[i]);
  }
  
  target += normal_lpdf(beta_fo | 0, 5);
  
  target += normal_lpdf(mu_d | 0, 1.5);

  target += cauchy_lpdf(sigma_a_rd | 0, 1);

  
  for (i in 1:K_rd){
    target += normal_lpdf(alpha_rd[i][1:L_rd[i]] | 0, sigma_a_rd[i]);
  }
  
  target += normal_lpdf(beta_fd | 0, 5);
  
  
  // Likelihood ----------------------------------------------------------------.
  
  target += reduce_sum(partial_sum_psi, V[occ], // V here just a dummy variable
                       grainsize1,
                       logit_psi[occ]);

  target += reduce_sum(partial_sum_p, y[occ_visit],
                       grainsize2,
                       logit_p_v[occ_visit]);
  

  for (i in not_occ){
    // Occurred and not observed
    target += log_sum_exp(bernoulli_logit_lpmf(1 | logit_psi[i])
    + bernoulli_logit_lpmf(0 | logit_p_v[start_i[i]:end_i[i]]),
    // Not occurred
    bernoulli_logit_lpmf(0 | logit_psi[i]));
  }
}

generated quantities {
  int z[P];         // occupancy indicator, 0/1

  z = bernoulli_logit_rng(logit_psi);
}

