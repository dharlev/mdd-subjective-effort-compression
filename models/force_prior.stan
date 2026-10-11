data {
  int<lower=1> N;
  int<lower=1> N_subj;
  int<lower=1> G;
  array[N] int<lower=1, upper=N_subj> subj;
  array[N_subj] int<lower=1, upper=G> group;
  vector[N] req_force;
  vector[N] prior_force;
  vector[N] trial_z;
  vector[N] rating_obs;
}

parameters {
  vector[G] alpha_group;
  vector[N_subj] z_alpha_subj;
  real<lower=0> sd_alpha_subj;

  vector[G] beta_effort_group;
  vector[G] beta_prior_group;
  vector[G] beta_interaction_group;
  real beta_trial;
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
      beta_effort_group[g] * req_force[n] +
      beta_prior_group[g] * prior_force[n] +
      beta_interaction_group[g] * prior_force[n] * req_force[n] +
      beta_trial * trial_z[n];
  }
}

model {
  alpha_group ~ normal(0, 20);
  z_alpha_subj ~ std_normal();
  sd_alpha_subj ~ normal(0, 10);

  beta_effort_group ~ normal(100, 40);
  beta_prior_group ~ normal(0, 40);
  beta_interaction_group ~ normal(0, 80);
  beta_trial ~ normal(0, 5);
  sigma_rating ~ normal(0, 15);

  rating_obs ~ normal(rating_mu, sigma_rating);
}

generated quantities {
  vector[N] log_lik;
  vector[N] rating_rep;

  for (n in 1:N) {
    log_lik[n] = normal_lpdf(rating_obs[n] | rating_mu[n], sigma_rating);
    rating_rep[n] = normal_rng(rating_mu[n], sigma_rating);
  }
}
