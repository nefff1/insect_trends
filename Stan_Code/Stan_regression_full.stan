data {
  // vector lengths ------------------------------------------------------------
  int<lower=0> N_obs; // number of observations
  int<lower=0> N_spec; // number of species (random effect)
  int<lower=0> N_reg; // number of regions (random effect)
  int<lower=0> N_int; // number of time intervals
  int<lower=0> N_group; // number of species groups

  // group ID ------------------------------------------------------------------
  
  int<lower=1,upper=N_group> group[N_obs];
  
  // random variables ----------------------------------------------------------
  int<lower=1,upper=N_spec> spec[N_obs];
  int<lower=1,upper=N_reg> reg[N_obs];

  // fixed variables -----------------------------------------------------------
  vector[N_obs] T_all;
  vector[N_obs] T_seas;
  vector[N_obs] P;
  vector[N_obs] AgrArea_prop;
  vector[N_obs] LSUGL;
  vector[N_obs] IAR;

  vector[N_obs] Tniche; 
  vector[N_obs] spec_index; 
  
  vector[N_obs] elevation;
  matrix[N_obs, N_int-1] year_start;
  
  // response variable ----------------------------------------------------------
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
  vector[N_spec] alpha_spec_s_AgrArea_prop;
  vector[N_spec] alpha_spec_s_LSUGL;
  vector[N_spec] alpha_spec_s_IAR;
  vector[N_spec] alpha_spec_s_T_all_AgrArea_prop;
  vector[N_spec] alpha_spec_s_T_seas_AgrArea_prop;
  vector[N_spec] alpha_spec_s_P_AgrArea_prop;
  vector[N_spec] alpha_spec_s_T_all_LSUGL;
  vector[N_spec] alpha_spec_s_T_seas_LSUGL;
  vector[N_spec] alpha_spec_s_P_LSUGL;
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
  real b_Tniche;
  real b_spec_index;
  real b_elevation;
  real b_T_all_Tniche;
  real b_T_seas_Tniche;
  real b_P_Tniche;
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
  real b_Tniche_elevation;
  real b_spec_index_elevation;
  real b_T_all_elevation;
  real b_T_seas_elevation;
  real b_P_elevation;
  real b_AgrArea_prop_elevation;
  real b_LSUGL_elevation;
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
    (b_AgrArea_prop + alpha_spec_s_AgrArea_prop[spec]) .* AgrArea_prop +
    (b_LSUGL + alpha_spec_s_LSUGL[spec]) .* LSUGL +
    (b_IAR + alpha_spec_s_IAR[spec]) .* IAR +
    // trait main effects
    b_Tniche * Tniche +
    b_spec_index * spec_index +
    // elevation main effect
    b_elevation * elevation +
    // Climate * trait interaction
    b_T_all_Tniche * T_all  .* Tniche +
    b_T_seas_Tniche * T_seas .* Tniche +
    b_P_Tniche * P .* Tniche +
    // Land use * trait interaction
    b_AgrArea_prop_spec_index * AgrArea_prop .* spec_index +
    b_LSUGL_spec_index * LSUGL .* spec_index +
    b_IAR_spec_index * IAR .* spec_index +
    // Climate * land use interactions + random slopes
    (b_T_all_AgrArea_prop + alpha_spec_s_T_all_AgrArea_prop[spec]) .* T_all .* AgrArea_prop +
    (b_T_seas_AgrArea_prop + alpha_spec_s_T_seas_AgrArea_prop[spec]) .* T_seas .* AgrArea_prop +
    (b_P_AgrArea_prop + alpha_spec_s_P_AgrArea_prop[spec]) .* P .* AgrArea_prop +
    (b_T_all_LSUGL + alpha_spec_s_T_all_LSUGL[spec]) .* T_all .* LSUGL +
    (b_T_seas_LSUGL + alpha_spec_s_T_seas_LSUGL[spec]) .* T_seas .* LSUGL +
    (b_P_LSUGL + alpha_spec_s_P_LSUGL[spec]) .* P .* LSUGL +
    (b_T_all_IAR + alpha_spec_s_T_all_IAR[spec]) .* T_all .* IAR +
    (b_T_seas_IAR + alpha_spec_s_T_seas_IAR[spec]) .* T_seas .* IAR +
    (b_P_IAR + alpha_spec_s_P_IAR[spec]) .* P .* IAR +
    // elevation * trait interaction
    b_Tniche_elevation * Tniche .* elevation +
    b_spec_index_elevation * spec_index .* elevation +
    // elevation * climate interaction
    b_T_all_elevation * T_all .* elevation +
    b_T_seas_elevation * T_seas .* elevation +
    b_P_elevation * P .* elevation +
    // elevation * land use interaction
    b_AgrArea_prop_elevation * AgrArea_prop .* elevation +
    b_LSUGL_elevation * LSUGL .* elevation +
    // year interval effect
    year_start * b_int;
  
  
}
model {
  // priors --------------------------------------------------------------------.
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
  target += normal_lpdf(b_Tniche | 0, 5);
  target += normal_lpdf(b_spec_index | 0, 5);
  target += normal_lpdf(b_elevation | 0, 5);
  target += normal_lpdf(b_T_all_Tniche | 0, 5);
  target += normal_lpdf(b_T_seas_Tniche | 0, 5);
  target += normal_lpdf(b_P_Tniche | 0, 5);
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
  target += normal_lpdf(b_Tniche_elevation | 0, 5);
  target += normal_lpdf(b_spec_index_elevation | 0, 5);
  target += normal_lpdf(b_T_all_elevation | 0, 5);
  target += normal_lpdf(b_T_seas_elevation | 0, 5);
  target += normal_lpdf(b_P_elevation | 0, 5);
  target += normal_lpdf(b_AgrArea_prop_elevation | 0, 5);
  target += normal_lpdf(b_LSUGL_elevation | 0, 5);
  target += normal_lpdf(b_int | 0, 5);

  target += cauchy_lpdf(sigma_y | 0, 25);
  
  // likelihood -----------------------------------------------------------------.

  target += normal_lpdf(y | y_hat, sigma_y); 
}

generated quantities {
// simulations from the posterior predictive distribution
  vector[N_obs] y_rep; // vector of same length as the data y
  for (i in 1:N_obs) 
    y_rep[i] = normal_rng(y_hat[i], sigma_y);
}
