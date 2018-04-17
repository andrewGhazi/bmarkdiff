data {
  int<lower=0> N; // number of regions
  int<lower=0> case_counts[N];
  int<lower=0> ctrl_counts[N];
}
parameters {
  vector<lower=0, upper=1>[N] case_p;
  vector<lower=0, upper=1>[N] ctrl_p;
  real<lower=0> a;
  real<lower=0> b;
}
model {
  a ~ gamma(20, 20);
  b ~ gamma(20, 20);
  ctrl_p ~ beta(a, b);
  case_p ~ beta(a, b);
  case_counts ~ binomial(12, case_p);
  ctrl_counts ~ binomial(12, ctrl_p);
}
