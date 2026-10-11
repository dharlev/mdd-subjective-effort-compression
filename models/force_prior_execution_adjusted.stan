data {
  int<lower=1> N;
  int<lower=1> N_subj;
  int<lower=1> G;
  array[N] int<lower=1, upper=N_subj> subj;
  array[N_subj] int<lower=1, upper=G> group;
  vector[N] req_force_z;
  vector[N] prior_force_z;
  vector[N] trial_z;
  vector[N] auc_force_z;
  vector[N] success_z;
  vector[N] age_z;
  vector[N] sex_male;
  vector[N] rating_z;
}

parameters {
  vector[G] alpha_group;
  vector[N_subj] z_alpha_subj;
  real<lower=0> sd_alpha_subj;

  vector[G] beta_effort_group;
  vector[G] beta_prior_group;
  vector[G] beta_interaction_group;
  real beta_trial;
  real beta_auc_force;
  real beta_success;
  real beta_age;
  real beta_sex;
  real<lower=0> sigma_rating;
}

transformed parameters {
  vector[N_subj] alpha_subj;
  vector[N] rating_mu;

  for (s in 1:N_subj) {
    alpha_subj[s] = alpha_group[group[s]] + sd_alpha_subj * z_alpha_subj[s];
  }

  for (n in 1:N) {
    int s = subj[n];
    int g = group[s];
    rating_mu[n] =
      alpha_subj[s] +
      beta_effort_group[g] * req_force_z[n] +
      beta_prior_group[g] * prior_force_z[n] +
      beta_interaction_group[g] * prior_force_z[n] * req_force_z[n] +
      beta_trial * trial_z[n] +
      beta_auc_force * auc_force_z[n] +
      beta_success * success_z[n] +
      beta_age * age_z[n] +
      beta_sex * sex_male[n];
  }
}

model {
  alpha_group ~ normal(0, 0.5);
  z_alpha_subj ~ std_normal();
  sd_alpha_subj ~ normal(0, 0.5);

  beta_effort_group ~ normal(0, 1);
  beta_prior_group ~ normal(0, 1);
  beta_interaction_group ~ normal(0, 1);
  beta_trial ~ normal(0, 0.3);
  beta_auc_force ~ normal(0, 1);
  beta_success ~ normal(0, 1);
  beta_age ~ normal(0, 0.5);
  beta_sex ~ normal(0, 0.5);
  sigma_rating ~ normal(0, 1);

  rating_z ~ normal(rating_mu, sigma_rating);
}
