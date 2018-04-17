data {
  int<lower=0> N; // number of regions
  int<lower=0> n_ctrl;
  int<lower=0> n_case;
  int<lower=0> case_counts[N];
  int<lower=0> ctrl_counts[N];
}
parameters {
  vector<lower=0, upper=1>[N] case_p;
  vector<lower=0, upper=1>[N] ctrl_p;
  real<lower=0> ctrl_a;
  real<lower=0> ctrl_b;
  real<lower=0> case_a;
  real<lower=0> case_b;
}
model {
  ctrl_a ~ gamma(20, 20);
  ctrl_b ~ gamma(20, 20);
  case_a ~ gamma(20, 20);
  case_b ~ gamma(20, 20);
  ctrl_p ~ beta(ctrl_a, ctrl_b);
  case_p ~ beta(case_a, case_b);
  case_counts ~ binomial(n_case, case_p);
  ctrl_counts ~ binomial(n_ctrl, ctrl_p);
}
generated quantities {
  vector[N] case_ctrl_diff;
  case_ctrl_diff = case_p - ctrl_p;
}
