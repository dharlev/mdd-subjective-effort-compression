data {
  int<lower=1> N;
  int<lower=1> N_subj;
  int<lower=1> G;

  array[N_subj] int<lower=1, upper=G> group;

  // Data must be sorted by subject and trial order.
  array[N_subj] int<lower=1, upper=N> subj_start;
  array[N_subj] int<lower=1, upper=N> subj_end;

  vector[N] effort_pct;
  vector[N] trial_z;
  vector[N] y;

  real beta_trial_fixed;
  real<lower=0> sigma_response_fixed;
}

parameters {
  vector[G] mu_prior0_group;
  vector[N_subj] z_mu_prior0;
  real<lower=0> sd_mu_prior0;

  vector<lower=log(1), upper=log(80)>[G] log_sigma_prior_group;
  vector<lower=log(1), upper=log(80)>[G] log_sigma_sensory_group;

  real logit_lambda;
}

transformed parameters {
  vector[N_subj] mu_prior0;
  vector<lower=0>[G] sigma_prior_group;
  vector<lower=0>[G] sigma_sensory_group;
  vector<lower=0, upper=1>[G] w_prior_group;
  real<lower=0, upper=1> lambda;

  lambda = inv_logit(logit_lambda);

  for (g in 1:G) {
    sigma_prior_group[g] = exp(log_sigma_prior_group[g]);
    sigma_sensory_group[g] = exp(log_sigma_sensory_group[g]);

    // Weight on prior mean: precision_prior / (precision_prior + precision_sensory).
    w_prior_group[g] =
      inv_logit(2 * (log_sigma_sensory_group[g] - log_sigma_prior_group[g]));
  }

  for (s in 1:N_subj) {
    mu_prior0[s] =
      mu_prior0_group[group[s]] +
      sd_mu_prior0 * z_mu_prior0[s];
  }
}

model {
  mu_prior0_group ~ normal(50, 20);
  sd_mu_prior0 ~ normal(0, 15);
  z_mu_prior0 ~ std_normal();

  log_sigma_prior_group ~ normal(log(20), 0.75);
  log_sigma_sensory_group ~ normal(log(20), 0.75);

  // Shared, modest learning from the previous posterior estimate.
  // Median is 0.25, with broad but regularized uncertainty.
  logit_lambda ~ normal(logit(0.25), 0.8);

  for (s in 1:N_subj) {
    int g = group[s];
    real m_prior = mu_prior0[s];
    real w_prior = w_prior_group[g];
    real sigma_y =
      sqrt(
        square((1 - w_prior) * sigma_sensory_group[g]) +
        square(sigma_response_fixed)
      );

    for (n in subj_start[s]:subj_end[s]) {
      real mu_y =
        w_prior * m_prior +
        (1 - w_prior) * effort_pct[n] +
        beta_trial_fixed * trial_z[n];

      y[n] ~ normal(mu_y, sigma_y);

      // Dynamic-lite update: only the prior mean changes, and only partially.
      m_prior = (1 - lambda) * m_prior + lambda * y[n];
    }
  }
}

generated quantities {
  vector[N] log_lik;
  vector[N] y_rep;
  vector[N] prior_mean_t;

  real delta_sigma_prior_g2_g1 = negative_infinity();
  real delta_sigma_sensory_g2_g1 = negative_infinity();
  real delta_w_prior_g2_g1 = negative_infinity();
  real delta_mu_prior0_g2_g1 = negative_infinity();

  if (G == 2) {
    delta_sigma_prior_g2_g1 = sigma_prior_group[2] - sigma_prior_group[1];
    delta_sigma_sensory_g2_g1 = sigma_sensory_group[2] - sigma_sensory_group[1];
    delta_w_prior_g2_g1 = w_prior_group[2] - w_prior_group[1];
    delta_mu_prior0_g2_g1 = mu_prior0_group[2] - mu_prior0_group[1];
  }

  for (s in 1:N_subj) {
    int g = group[s];
    real m_prior = mu_prior0[s];
    real m_prior_rep = mu_prior0[s];
    real w_prior = w_prior_group[g];
    real sigma_y =
      sqrt(
        square((1 - w_prior) * sigma_sensory_group[g]) +
        square(sigma_response_fixed)
      );

    for (n in subj_start[s]:subj_end[s]) {
      real mu_y =
        w_prior * m_prior +
        (1 - w_prior) * effort_pct[n] +
        beta_trial_fixed * trial_z[n];

      real mu_y_rep =
        w_prior * m_prior_rep +
        (1 - w_prior) * effort_pct[n] +
        beta_trial_fixed * trial_z[n];

      prior_mean_t[n] = m_prior;
      log_lik[n] = normal_lpdf(y[n] | mu_y, sigma_y);
      y_rep[n] = normal_rng(mu_y_rep, sigma_y);

      m_prior = (1 - lambda) * m_prior + lambda * y[n];
      m_prior_rep = (1 - lambda) * m_prior_rep + lambda * y_rep[n];
    }
  }
}
