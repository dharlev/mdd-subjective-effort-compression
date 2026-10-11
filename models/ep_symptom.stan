data {
  int<lower=1> N_subj;
  int<lower=1> N_obs;
  array[N_obs] int<lower=1, upper=N_subj> subj;
  vector[N_obs] rating_z;
  vector[N_obs] effort_z;
  vector[N_obs] trial_z;
  vector[N_obs] phq_z;
  vector[N_obs] age_z;
  vector[N_obs] sex_male;
}
parameters {
  real alpha;
  real beta_effort;
  real beta_phq;
  real beta_phq_effort;
  real beta_age;
  real beta_sex;
  real beta_trial;
  vector<lower=0>[2] tau;
  cholesky_factor_corr[2] L_Omega;
  matrix[2, N_subj] z_u;
  real<lower=0> sigma_y;
}
transformed parameters {
  matrix[2, N_subj] u = diag_pre_multiply(tau, L_Omega) * z_u;
}
model {
  alpha ~ normal(0, 1);
  beta_effort ~ normal(0.75, 0.20);
  beta_phq ~ normal(0, 0.5);
  beta_phq_effort ~ normal(0, 0.5);
  beta_age ~ normal(0, 0.5);
  beta_sex ~ normal(0, 0.5);
  beta_trial ~ normal(0, 0.5);
  tau[1] ~ exponential(1);
  tau[2] ~ exponential(8);
  L_Omega ~ lkj_corr_cholesky(2);
  to_vector(z_u) ~ std_normal();
  sigma_y ~ exponential(1);

  for (n in 1:N_obs) {
    int s = subj[n];
    real mu = alpha + u[1, s]
      + (beta_effort + u[2, s]) * effort_z[n]
      + beta_phq * phq_z[n]
      + beta_phq_effort * phq_z[n] * effort_z[n]
      + beta_age * age_z[n]
      + beta_sex * sex_male[n]
      + beta_trial * trial_z[n];
    rating_z[n] ~ normal(mu, sigma_y);
  }
}
generated quantities {
  corr_matrix[2] Omega = multiply_lower_tri_self_transpose(L_Omega);
}

