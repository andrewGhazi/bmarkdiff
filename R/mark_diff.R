#' Estimate the joint for a mark
#'
#' Fits a stan model to estimate the joint distribution on mark probability between two groups.
#'
#' This function takes in a file of marks. The input should be a tsv, with the first column giving a region id and the remainder giving peak calls for each subject. Non-zero integer counts are converted to TRUE, while 0's are converted to FALSE.
#'     The package fits a joint beta distribution on the peak probabilities in the two groups.
#'
#' @param input character - path to an mark input file
#' @param n_case integer - the number of columns (after the id column) that are cases. The remainder are considered controls.
#' @param n_iter integer - number of posterior samples to draw
#' @param use_vb logical - use a variational approximation
#' @param use_subsamp logical - use a random subsample of regions for speed
#' @param n_subsamp logical - number of regions to subsample (if applicable)
#' @param n_cores integer - number of cores to use when using multiple MCMC chains
#' @param n_chains integer - number of MCMC chains
#' @param n_warmup integer - number of MCMC warmup iterations to use
#' @param hdi_content real(0,1) - the posterior probability content to include in the posterior density interval cutoff used to make binary calls of differential marks
#' @return a stanfit object containing the marginal joint prior on case_a, case_b, ctrl_a, and ctrl_b
#' @importFrom magrittr %>%
#' @export
estimate_joint_prior = function(input,
                                n_case = 12,
                                n_iter = 3834,
                                use_vb = TRUE,
                                use_subsamp = TRUE,
                                n_subsamp = 500,
                                n_cores = 1,
                                n_chains = 3,
                                n_warmup = 500,
                                hdi_content = .99) {

  calls = readr::read_tsv(input,
                          col_names = FALSE) %>%
    purrr::set_names(c('region',
                       paste0('case_', 1:n_case),
                       paste0('ctrl_', 1:(ncol(.) - 1 - n_case)))) %>%f
    tidyr::gather(subj, peak, -region) %>%
    dplyr::mutate(peak = peak > 0)

  n_ctrl = calls$subj %>%
    unique %>%
    grepl('ctrl', .) %>%
    sum

  n_peaks = calls %>%
    dplyr::group_by(region) %>%
    dplyr::summarise(case = sum(peak[grepl('case', subj)]),
                     ctrl = sum(peak[grepl('ctrl', subj)]))

  if (use_subsamp) {
    region_samp = dplyr::sample_n(n_peaks, size = n_subsamp)
  } else {
    region_samp = n_peaks
  }


  use_model1 = TRUE
  if (use_vb){
    model_res = vb_approx(region_samp)
  } else if(use_model1){
    model_res = fit_model1(region_samp,
                            n_iter = n_iter,
                            n_cores = n_cores,
                            n_chains = n_chains,
                            n_warmup = n_warmup)
  } else {
    model_res = fit_model2(region_samp,
                            n_iter = n_iter,
                            n_cores = n_cores,
                            n_chains = n_chains,
                            n_warmup = n_warmup)
  }

  proportion_prior = model1_res %>%
    stan_to_tibble()

  # outcome_mat = matrix(rep(FALSE, n_ctrl*n_case), nrow = n_case, ncol = n_ctrl)
  # colnames(outcome_mat) = as.character(1:n_ctrl)
  # rownames(outcome_mat) = as.character(1:n_case)

  outcome_df = expand.grid(0:n_ctrl, 0:n_case) %>%
    dplyr::as_tibble() %>%
    purrr::set_names('ctrl', 'case') %>%
    dplyr::mutate(diff_post = purrr::map2(ctrl, case,
                                         prior_to_post,
                                         prior = proportion_prior,
                                         n_ctrl = n_ctrl,
                                         n_case = n_case),
                  different_by_hdi = purrr::map_lgl(diff_post,
                                            ~get_hdi_call(.x, hdi_content = hdi_content)))

  differentially_marked_regions = n_peaks %>%
    dplyr::left_join(outcome_df %>% dplyr::select(-diff_post), by = c('case', 'ctrl')) %>%
    dplyr::filter(different_by_hdi)

  res = list(outcome_df,
             differentially_marked_regions)

  return(res)
}

get_hdi_call = function(post_tbl,
                        hdi_content) {
  post_tbl %>%
    dplyr::pull(case_ctrl_diff) %>%
    coda::mcmc() %>%
    coda::HPDinterval(prob = hdi_content) %>%
    {!dplyr::between(0, .[1], .[2])}


}

prior_to_post = function(prior,
                         ctrl_count,
                         case_count,
                         n_ctrl,
                         n_case){
  prior %>%
    mutate(case_a = case_a + case_count,
           case_b = case_b + n_case - case_count,
           ctrl_a = ctrl_a + ctrl_count,
           ctrl_b = ctrl_b + n_ctrl - ctrl_count,
           case_p = map2_dbl(case_a, case_b, ~rbeta(1, .x, .y)),
           ctrl_p = map2_dbl(ctrl_a, ctrl_b, ~rbeta(1, .x, .y)),
           case_ctrl_diff = case_p - ctrl_p)
}

#' Variational approximation of mark differences
#'
#' This functional uses \code{rstan::vb} to run a variational approximation on the probability
vb_approx = function(mark_dat) {

  data_list = list(N = nrow(mark_dat),
                   case_counts = mark_dat$case,
                   ctrl_counts = mark_dat$ctrl)

  vb_res = rstan::vb(object = stanmodels$mark_diff,
                     data = data_list,
                     output_samples = 5000,
                     pars = c('case_a', 'case_b', 'ctrl_a', 'ctrl_b'),
                     include = TRUE)

  return(vb_res)
}

fit_model1 = function(mark_dat,
                      n_iter = 3834,
                      n_cores = 1,
                      n_chains = 3,
                      n_warmup = 500) {

  data_list = list(N = nrow(mark_dat),
                   case_counts = mark_dat$case,
                   ctrl_counts = mark_dat$ctrl)

  mcmc_res = sampling(stanmodels$mark_diff,
                      data = data_list,
                      chains = n_chains,
                      iter = n_iter,
                      warmup = n_warmup,
                      cores = n_cores,
                      control = list(max_treedepth = 15),
                      pars = c('case_a', 'case_b', 'ctrl_a', 'ctrl_b'),
                      include = TRUE)

  return(mcmc_res)
}

fit_model2 = function(mark_dat,
                      n_iter = 3834,
                      n_cores = 1,
                      n_chains = 3,
                      n_warmup = 500) {

  data_list = list(N = nrow(mark_dat),
                   case_counts = mark_dat$case,
                   ctrl_counts = mark_dat$ctrl)

  mcmc_res = sampling(stanmodels$mark_diff2,
                      data = data_list,
                      chains = n_chains,
                      iter = n_iter,
                      warmup = n_warmup,
                      cores = n_cores,
                      control = list(max_treedepth = 15),
                      pars = c('a', 'b'),
                      include = TRUE)

  return(mcmc_res)
}

res_to_df = function(stan_param, name){
  if(class(stan_param) == 'array'){
    dplyr::data_frame(x = stan_param) %>%
      purrr::set_names(name)
  } else {
    dplyr::as_tibble(stan_param) %>%
      purrr::set_names(paste(name, 1:ncol(stan_param), sep = '_'))
  }
}

stan_to_tibble = function(stan_res){
  res = rstan::extract(stan_res)

  par_names = names(res)

  purrr::map2(res,
       par_names,
       res_to_df) %>%
    dplyr::bind_cols
}
