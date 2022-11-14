data {
  // vector lengths ------------------------------------------------------------
  int<lower=0> N_obs; // number of observations
  int<lower=0> N_obs_agrspec; // number of observations of "agricultural"" species only 
  int<lower=0> N_spec; // number of species (random effect)
  int<lower=0> N_agrspec; // number of "agricultural"" species 
  int<lower=0> N_reg; // number of regions (random effect)
  int<lower=0> N_int; // number of time intervals
  int<lower=0> N_group; // number of species groups

  // group ID ------------------------------------------------------------------
  
  int<lower=1,upper=N_group> group[N_obs];
  
  // random variables ----------------------------------------------------------
  int<lower=1,upper=N_spec> spec[N_obs]; // species index
  int<lower=1,upper=N_spec> spec_agrspec[N_obs_agrspec]; // species index only accounting for "agricultural" species
  int<lower=1,upper=N_reg> reg[N_obs];

  // fixed variables -----------------------------------------------------------
  vector[N_obs] T_all;
  vector[N_obs] T_seas;
  vector[N_obs] P;
  vector[N_obs_agrspec] AgrArea_prop;
  vector[N_obs_agrspec] LSUGL;
  vector[N_obs] IAR;

  vector[N_obs] Tavg; 
  vector[N_obs] spec_index; 
  
  vector[N_obs] height;
  matrix[N_obs, N_int-1] year_start;
  
  int index_agrspec[N_obs_agrspec]; // index of "agricultural" species' obserations
  
  // response variable ---------------------------------------------------------
  vector[N_obs] y; 
  real y_mean; // mean response
}
parameters {
  vector[N_group] mu_a; // global intercept per group
  
  // random intercept species
  vector[N_spec] alpha_spec_int;
  // random slopes species
  vector[N_spec] alpha_spec_s_T_all;
  vector[N_spec] alpha_spec_s_T_seas;
  vector[N_spec] alpha_spec_s_P;
  vector[N_agrspec] alpha_spec_s_AgrArea_prop;
  vector[N_agrspec] alpha_spec_s_LSUGL;
  vector[N_spec] alpha_spec_s_IAR;
  vector[N_agrspec] alpha_spec_s_T_all_AgrArea_prop;
  vector[N_agrspec] alpha_spec_s_T_seas_AgrArea_prop;
  vector[N_agrspec] alpha_spec_s_P_AgrArea_prop;
  vector[N_agrspec] alpha_spec_s_T_all_LSUGL;
  vector[N_agrspec] alpha_spec_s_T_seas_LSUGL;
  vector[N_agrspec] alpha_spec_s_P_LSUGL;
  vector[N_spec] alpha_spec_s_T_all_IAR;
  vector[N_spec] alpha_spec_s_T_seas_IAR;
  vector[N_spec] alpha_spec_s_P_IAR;
  
  // random intercept region
  vector[N_reg] alpha_reg;
  
  // slope coefficients
  real b_T_all;
  real b_T_seas;
  real b_P;
  real b_AgrArea_prop;
  real b_LSUGL;
  real b_IAR;
  real b_Tavg;
  real b_spec_index;
  real b_height;
  real b_T_all_Tavg;
  real b_T_seas_Tavg;
  real b_P_Tavg;
  real b_AgrArea_prop_spec_index;
  real b_LSUGL_spec_index;
  real b_IAR_spec_index;
  real b_T_all_AgrArea_prop;
  real b_T_seas_AgrArea_prop;
  real b_P_AgrArea_prop;
  real b_T_all_LSUGL;
  real b_T_seas_LSUGL;
  real b_P_LSUGL;
  real b_T_all_IAR;
  real b_T_seas_IAR;
  real b_P_IAR;
  real b_Tavg_height;
  real b_spec_index_height;
  real b_T_all_height;
  real b_T_seas_height;
  real b_P_height;
  real b_AgrArea_prop_height;
  real b_LSUGL_height;
  vector[N_int-1] b_int;

  // sigma of randoms
  real<lower=0,upper=100> sigma_a_spec_int; 
  real<lower=0,upper=100> sigma_a_spec_slope[15]; 
  real<lower=0,upper=100> sigma_a_reg; 

  // sigma of link
  real<lower=0,upper=100> sigma_y;
}
transformed parameters{
  vector[N_obs] y_hat;
  
    y_hat = 
    // Intercept
    mu_a[group] +
    // random intercepts
    alpha_spec_int[spec] +
    alpha_reg[reg] +
    // climate main effects + random slopes
    (b_T_all + alpha_spec_s_T_all[spec]) .* T_all +
    (b_T_seas + alpha_spec_s_T_seas[spec]) .* T_seas +
    (b_P + alpha_spec_s_P[spec]) .* P +
    // land use main effects + random slopes
    (b_IAR + alpha_spec_s_IAR[spec]) .* IAR +
    // trait main effects
    b_Tavg * Tavg +
    b_spec_index * spec_index +
    // height main effect
    b_height * height +
    // Climate * trait interaction
    b_T_all_Tavg * T_all  .* Tavg +
    b_T_seas_Tavg * T_seas .* Tavg +
    b_P_Tavg * P .* Tavg +
    // Land use * trait interaction
    b_IAR_spec_index * IAR .* spec_index +
    // Climate * land use interactions + random slopes
    (b_T_all_IAR + alpha_spec_s_T_all_IAR[spec]) .* T_all .* IAR +
    (b_T_seas_IAR + alpha_spec_s_T_seas_IAR[spec]) .* T_seas .* IAR +
    (b_P_IAR + alpha_spec_s_P_IAR[spec]) .* P .* IAR +
    // height * trait interaction
    b_Tavg_height * Tavg .* height +
    b_spec_index_height * spec_index .* height +
    // height * climate interaction
    b_T_all_height * T_all .* height +
    b_T_seas_height * T_seas .* height +
    b_P_height * P .* height +
    // year interval effect
    year_start * b_int;
    
    // "Agricultural" species only part ----------------------------------------
    
    y_hat[index_agrspec] +=
    // land use main effects + random slopes
    (b_AgrArea_prop + alpha_spec_s_AgrArea_prop[spec_agrspec]) .* AgrArea_prop +
    (b_LSUGL + alpha_spec_s_LSUGL[spec_agrspec]) .* LSUGL +
    // Land use * trait interaction
    b_AgrArea_prop_spec_index * AgrArea_prop .* spec_index[index_agrspec] +
    b_LSUGL_spec_index * LSUGL .* spec_index[index_agrspec] +
    // Climate * land use interactions + random slopes
    (b_T_all_AgrArea_prop + alpha_spec_s_T_all_AgrArea_prop[spec_agrspec]) .* T_all[index_agrspec] .* AgrArea_prop +
    (b_T_seas_AgrArea_prop + alpha_spec_s_T_seas_AgrArea_prop[spec_agrspec]) .* T_seas[index_agrspec] .* AgrArea_prop +
    (b_P_AgrArea_prop + alpha_spec_s_P_AgrArea_prop[spec_agrspec]) .* P[index_agrspec] .* AgrArea_prop +
    (b_T_all_LSUGL + alpha_spec_s_T_all_LSUGL[spec_agrspec]) .* T_all[index_agrspec] .* LSUGL +
    (b_T_seas_LSUGL + alpha_spec_s_T_seas_LSUGL[spec_agrspec]) .* T_seas[index_agrspec] .* LSUGL +
    (b_P_LSUGL + alpha_spec_s_P_LSUGL[spec_agrspec]) .* P[index_agrspec] .* LSUGL +
    // height * land use interaction
    b_AgrArea_prop_height * AgrArea_prop .* height[index_agrspec] +
    b_LSUGL_height * LSUGL .* height[index_agrspec];
    
}
model {
  // priors --------------------------------------------------------------------
  target += cauchy_lpdf(sigma_a_spec_int | 0, 1);
  target += cauchy_lpdf(sigma_a_spec_slope | 0, 1);
  target += cauchy_lpdf(sigma_a_reg | 0, 1);
  
  target += normal_lpdf(alpha_spec_int | 0, sigma_a_spec_int);
  target += normal_lpdf(alpha_reg | 0, sigma_a_reg);
  
  target += normal_lpdf(alpha_spec_s_T_all | 0, sigma_a_spec_slope[1]);
  target += normal_lpdf(alpha_spec_s_T_seas | 0, sigma_a_spec_slope[2]);
  target += normal_lpdf(alpha_spec_s_P | 0, sigma_a_spec_slope[3]);
  target += normal_lpdf(alpha_spec_s_AgrArea_prop | 0, sigma_a_spec_slope[4]);
  target += normal_lpdf(alpha_spec_s_LSUGL | 0, sigma_a_spec_slope[5]);
  target += normal_lpdf(alpha_spec_s_IAR | 0, sigma_a_spec_slope[6]);
  target += normal_lpdf(alpha_spec_s_T_all_AgrArea_prop | 0, sigma_a_spec_slope[7]);
  target += normal_lpdf(alpha_spec_s_T_seas_AgrArea_prop | 0, sigma_a_spec_slope[8]);
  target += normal_lpdf(alpha_spec_s_P_AgrArea_prop | 0, sigma_a_spec_slope[9]);
  target += normal_lpdf(alpha_spec_s_T_all_LSUGL | 0, sigma_a_spec_slope[10]);
  target += normal_lpdf(alpha_spec_s_T_seas_LSUGL | 0, sigma_a_spec_slope[11]);
  target += normal_lpdf(alpha_spec_s_P_LSUGL | 0, sigma_a_spec_slope[12]);
  target += normal_lpdf(alpha_spec_s_T_all_IAR | 0, sigma_a_spec_slope[13]);
  target += normal_lpdf(alpha_spec_s_T_seas_IAR | 0, sigma_a_spec_slope[14]);
  target += normal_lpdf(alpha_spec_s_P_IAR | 0, sigma_a_spec_slope[15]);

  target += normal_lpdf(mu_a | y_mean, 5);

  target += normal_lpdf(b_T_all | 0, 5);
  target += normal_lpdf(b_T_seas | 0, 5);
  target += normal_lpdf(b_P | 0, 5);
  target += normal_lpdf(b_AgrArea_prop | 0, 5);
  target += normal_lpdf(b_LSUGL | 0, 5);
  target += normal_lpdf(b_IAR | 0, 5);
  target += normal_lpdf(b_Tavg | 0, 5);
  target += normal_lpdf(b_spec_index | 0, 5);
  target += normal_lpdf(b_height | 0, 5);
  target += normal_lpdf(b_T_all_Tavg | 0, 5);
  target += normal_lpdf(b_T_seas_Tavg | 0, 5);
  target += normal_lpdf(b_P_Tavg | 0, 5);
  target += normal_lpdf(b_AgrArea_prop_spec_index | 0, 5);
  target += normal_lpdf(b_LSUGL_spec_index | 0, 5);
  target += normal_lpdf(b_IAR_spec_index | 0, 5);
  target += normal_lpdf(b_T_all_AgrArea_prop | 0, 5);
  target += normal_lpdf(b_T_seas_AgrArea_prop | 0, 5);
  target += normal_lpdf(b_P_AgrArea_prop | 0, 5);
  target += normal_lpdf(b_T_all_LSUGL | 0, 5);
  target += normal_lpdf(b_T_seas_LSUGL | 0, 5);
  target += normal_lpdf(b_P_LSUGL | 0, 5);
  target += normal_lpdf(b_T_all_IAR | 0, 5);
  target += normal_lpdf(b_T_seas_IAR | 0, 5);
  target += normal_lpdf(b_P_IAR | 0, 5);
  target += normal_lpdf(b_Tavg_height | 0, 5);
  target += normal_lpdf(b_spec_index_height | 0, 5);
  target += normal_lpdf(b_T_all_height | 0, 5);
  target += normal_lpdf(b_T_seas_height | 0, 5);
  target += normal_lpdf(b_P_height | 0, 5);
  target += normal_lpdf(b_AgrArea_prop_height | 0, 5);
  target += normal_lpdf(b_LSUGL_height | 0, 5);
  target += normal_lpdf(b_int | 0, 5);

  target += cauchy_lpdf(sigma_y | 0, 25);
  
  // likelihood ----------------------------------------------------------------

  target += normal_lpdf(y | y_hat, sigma_y); 
}

generated quantities {
// simulations from the posterior predictive distribution
  vector[N_obs] y_rep; // vector of same length as the data y
  for (i in 1:N_obs) 
    y_rep[i] = normal_rng(y_hat[i], sigma_y);
}
