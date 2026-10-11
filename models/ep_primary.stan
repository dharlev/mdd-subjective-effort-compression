data {
  int<lower=1> N_subj;
  int<lower=1> N_obs;
  int<lower=1> N_group;
  array[N_obs] int<lower=1, upper=N_subj> subj;
  array[N_subj] int<lower=1, upper=N_group> subj_group;
  vector[N_obs] rating_z;
  vector[N_obs] effort_z;
  vector[N_obs] trial_z;
}

parameters {
  vector[N_group] mu_alpha_group;
  vector[N_group] mu_beta_group;
  vector<lower=0>[2] tau;
  cholesky_factor_corr[2] L_Omega;
  matrix[2, N_subj] z_u;

  real beta_trial;
  real<lower=0> sigma_y;
}

transformed parameters {
  matrix[2, N_subj] u_raw;
  vector[N_subj] alpha;
  vector[N_subj] beta;

  u_raw = diag_pre_multiply(tau, L_Omega) * z_u;
  for (s in 1:N_subj) {
    alpha[s] = mu_alpha_group[subj_group[s]] + u_raw[1, s];
    beta[s] = mu_beta_group[subj_group[s]] + u_raw[2, s];
  }
}

model {
  mu_alpha_group ~ normal(0, 1);
  mu_beta_group ~ normal(0.75, 0.20);

  tau[1] ~ exponential(1);
  tau[2] ~ exponential(8);
  L_Omega ~ lkj_corr_cholesky(2);
  to_vector(z_u) ~ std_normal();

  beta_trial ~ normal(0, 0.5);
  sigma_y ~ exponential(1);

  for (n in 1:N_obs) {
    int s = subj[n];
    rating_z[n] ~ normal(alpha[s] + beta[s] * effort_z[n] + beta_trial * trial_z[n], sigma_y);
  }
}

generated quantities {
  corr_matrix[2] Omega;
  real corr_intercept_slope;

  Omega = multiply_lower_tri_self_transpose(L_Omega);
  corr_intercept_slope = Omega[1, 2];
}
