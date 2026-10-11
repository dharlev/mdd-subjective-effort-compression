
data {
  int<lower=1> N_subj;
  int<lower=1> N_obs;
  int<lower=1> G;
  array[N_obs] int<lower=1, upper=N_subj> subj;
  array[N_subj] int<lower=1, upper=G> group_id;
  array[N_obs] int<lower=0, upper=1> y;
  vector[N_obs] rew1;
  vector[N_obs] eff1;
  vector[N_obs] rew2;
  vector[N_obs] eff2;
}
parameters {
  vector[G] mu_log_rew;
  vector[G] mu_log_eff;
  vector<lower=0>[G] sigma_log_rew;
  vector<lower=0>[G] sigma_log_eff;
  vector[N_subj] z_rew;
  vector[N_subj] z_eff;
}
transformed parameters {
  vector[N_subj] log_rew;
  vector[N_subj] log_eff;
  vector<lower=0>[N_subj] rewSens;
  vector<lower=0>[N_subj] effSens;

  for (s in 1:N_subj) {
    int g = group_id[s];
    log_rew[s] = mu_log_rew[g] + sigma_log_rew[g] * z_rew[s];
    log_eff[s] = mu_log_eff[g] + sigma_log_eff[g] * z_eff[s];
    rewSens[s] = exp(log_rew[s]);
    effSens[s] = exp(log_eff[s]);
  }
}
model {
  mu_log_rew ~ normal(0, 1);
  mu_log_eff ~ normal(0, 1);
  sigma_log_rew ~ exponential(1);
  sigma_log_eff ~ exponential(1);
  z_rew ~ std_normal();
  z_eff ~ std_normal();

  for (i in 1:N_obs) {
    int s = subj[i];
    real V1 = rewSens[s] * rew1[i] - effSens[s] * square(eff1[i]);
    real V2 = rewSens[s] * rew2[i] - effSens[s] * square(eff2[i]);
    y[i] ~ bernoulli_logit(V2 - V1);
  }
}
generated quantities {
  vector<lower=0>[G] mu_rew;
  vector<lower=0>[G] mu_eff;
  real ybar_obs = 0;
  real ybar_rep = 0;
  array[N_obs] int y_rep;

  for (g in 1:G) {
    mu_rew[g] = exp(mu_log_rew[g]);
    mu_eff[g] = exp(mu_log_eff[g]);
  }

  for (i in 1:N_obs) {
    int s = subj[i];
    real V1 = rewSens[s] * rew1[i] - effSens[s] * square(eff1[i]);
    real V2 = rewSens[s] * rew2[i] - effSens[s] * square(eff2[i]);
    real p = inv_logit(V2 - V1);
    y_rep[i] = bernoulli_rng(p);
    ybar_obs += y[i];
    ybar_rep += y_rep[i];
  }

  ybar_obs /= N_obs;
  ybar_rep /= N_obs;
}

