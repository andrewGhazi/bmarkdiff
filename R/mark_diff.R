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
#' @param n_warmup integer - number of MCMC warmup iterations to use#'
#' @return a stanfit object containing the marginal joint prior on case_a, case_b, ctrl_a, and ctrl_b
#' @importFrom magrittr %>%
#' @export
estimate_joint_prior = function(input,
                                n_case = 12,
                                n_iter = 3834,
                                use_vb = TRUE,
                                use_subsamp = TRUE,
                                n_subsamp = 1000,
                                n_cores = 1,
                                n_chains = 3,
                                n_warmup = 500,) {

  calls = readr::read_tsv(input,
                          col_names = FALSE) %>%
    purrr::set_names(c('region',
                       paste0('case_', 1:n_case),
                       paste0('ctrl_', 1:(ncol(.) - 1 - n_case)))) %>%
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

  if (use_vb){
    vb_approx(region_samp)
  } else {
    mcmc_fit(region_samp,
             n_iter = n_iter,
             n_cores = n_cores,
             n_chains = n_chains,
             n_warmup = n_warmup)
  }

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

mcmc_fit = function(mark_dat,
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
