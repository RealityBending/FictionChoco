library(brms)
library(cmdstanr)
library(rstan)

combine_and_save <- function(name = "sample1_zoib") {
  files <- list.files('./models/', pattern = paste0(name, '_.*rds$'), full.names = TRUE)
  m <- brms::combine_models(mlist = lapply(files, readRDS))

  # m <- brms::add_criterion(m, criterion="waic", model_name = name)
  # loo1 <- loo::loo_subsample(m, observations = 10, loo_approximation_draws = 100, cores = 4)
  # fit_laplace <- rstan::optimizing(m, data = insight::get_data(m), draws = 1000, importance_resampling = TRUE)
  # parameter_draws_laplace <- fit_laplace$theta_tilde # draws from approximate posterior
  # log_p <- fit_laplace$log_p # log density of the posterior
  # log_g <- fit_laplace$log_g # log density of the approximation


  saveRDS(m, paste0("models/", name, ".rds"))
}

combine_and_save("sample1_zoib")
combine_and_save("sample1_bext")
combine_and_save("sample1_choco")
combine_and_save("sample1_chocoattractiveness")
