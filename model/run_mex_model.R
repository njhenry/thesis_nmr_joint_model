## #######################################################################################
##
## RUN MEXICO NMR MODEL
##
## Author: Nat Henry, github: @njhenry
## Created: 19 August 2021
## Purpose: Run Mexico NMR model
##
## #######################################################################################


load_libs <- c(
  'data.table', 'glue', 'stats', 'tictoc', 'Matrix', 'TMB', 'INLA', 'optimx', 'yaml'
)
invisible(lapply(load_libs, library, character.only=T))


## Input settings
## TODO: Make CLI
run_date <- '20210930_allpois_narrow_sds_2'
data_version <- '20210818'
holdout <- 0
use_covs <- c('intercept', paste0(c('yrschool','lowwage','nohca'), '_norm'))

sim <- FALSE
sim_args <- list(
  mort_coefs = c(
    intercept = -5.5, yrschool_norm = -.25, lowwage_norm = 0.2, nohca_norm = .5,
    electric = -.3, refrig = -0.1
  ),
  add_vr_bias = TRUE, add_bh_bias = TRUE, bh_bias = 1.05
)
# Standard deviations of bias terms in low-exclusion, med-exclusion, and high-exclusion
#  municipalities
# vr_bias_thresholds <- c(0.048, 0.207, 0.821) # Wide SDs
vr_bias_thresholds <- rep(0.0001, 3) # Narrow SDs
starting_log_precisions <- log(1/vr_bias_thresholds**2)

# Get input and output filepaths
repo_dir <- '{REDACTED}'
config <- yaml::read_yaml(file.path(repo_dir, 'mexico/config.yaml'))
source(file.path(repo_dir, 'mexico/model/model_functions.R'))
template_fp <- file.path(repo_dir, 'mexico/model/vr_bh_tmb_model.cpp')

data_dir <- file.path(config$dirs$prepped_data, data_version)
model_dir <- file.path(config$dirs$runs, run_date)
dir.create(model_dir, showWarnings = FALSE)

## SET INPUTS --------------------------------------------------------------------------->

# Load input data
prepped_data <- fread(file.path(data_dir, 'prepped_data_stable.csv'))
# Drop clearly erroneous observations from VR
prepped_data[vr_deaths / vr_births > .1, c('vr_deaths', 'vr_births') := 0]
adjmat <- readRDS(file.path(data_dir, 'adjacency_matrix.RDS'))
use_covs <- intersect(use_covs, colnames(prepped_data))

q_icar <- icar_precision_from_adjacency(adjmat)

# IF SIMULATING: CREATE SIMULATED MORTALITY, VR BIAS, AND OBSERVATION DATA
if(sim){
  # Simulate true mortality
  mort_coefs <- sim_args$mort_coefs
  log_mort <- rep(0., nrow(prepped_data))
  for(cname in names(mort_coefs)){
    log_mort <- log_mort + prepped_data[[cname]] * mort_coefs[cname]
  }
  prepped_data$sim_mort <- exp(log_mort)
  # Simulate VR Bias
  prepped_data[, sim_vr_bias := 0 ]
  # The simulated SDs should be plausible, between 50-100% of the 'high' threshold
  sim_threshold_ratios <- runif(n=3, min=.5)
  if(sim_args$add_vr_bias){
    for(idx_bias in unique(prepped_data$excl_group)){
      which_to_sim <- (prepped_data$excl_group == idx_bias)
      prepped_data[which_to_sim,]$sim_vr_bias <- exp(rnorm(
        n = sum(which_to_sim),
        mean = 0,
        sd = vr_bias_thresholds[idx_bias + 1] * sim_threshold_ratios[idx_bias + 1]
      ))
    }
  }
  knitr::kable(prepped_data[
      , .(hilo = range(sim_vr_bias), qq=quantile(sim_vr_bias, probs=c(.1, .9))), by=excl_group
    ][order(excl_group)]
  )
  # Add BH bias, if specified
  prepped_data[, sim_bh_bias := 1 ]
  if(sim_args$add_bh_bias) prepped_data$sim_bh_bias <- sim_args$bh_bias
  # Simulate observations
  prepped_data$cen_deaths <- rpois(
    n = prepped_data$cen_births,
    lambda = prepped_data$cen_births * prepped_data$sim_mort * prepped_data$sim_bh_bias
  )
  prepped_data$vr_deaths <- rpois(
    n = prepped_data$vr_births,
    lambda = prepped_data$vr_births * prepped_data$sim_mort * prepped_data$sim_vr_bias
  )
}

# Get MAP fixed effects estimates
map_est <- find_glm_map_parameter_estimates(
  in_data = copy(prepped_data), events_field = 'cen_deaths', exposure_field = 'cen_births',
  covar_names = use_covs, distribution_family = 'binomial'
)
saveRDS(map_est, file = file.path(model_dir, 'map_est.RDS'))

tmb_data_stack <- list(
  holdout = holdout,
  vr_births = prepped_data$vr_births, vr_deaths = prepped_data$vr_deaths,
  bh_births = prepped_data$cen_births, bh_deaths = prepped_data$cen_deaths,
  X_ij = as.matrix(prepped_data[, ..use_covs]),
  idx_loc = prepped_data$uid, idx_vr_bias = prepped_data$excl_group,
  idx_holdout_vr = prepped_data$idx_holdout_vr,
  idx_holdout_bh = prepped_data$idx_holdout_bh,
  bias_thresh_vec = vr_bias_thresholds,
  Q_icar = q_icar
)
params_list <- list(
  beta_covs = map_est$fixed_effects_map,
  log_tau_loc = -1., logit_phi_loc = 1.,
  vr_bias_logprec_vec = starting_log_precisions,
  log_bh_bias = 0.,
  Z_loc = rep(0, nrow(adjmat)), log_vr_bias = rep(0, nrow(prepped_data))
)


## RUN TMB MODEL ------------------------------------------------------------------------>

model_fit <- setup_run_tmb(
  tmb_data_stack = tmb_data_stack,
  params_list = params_list,
  tmb_random = c('Z_loc', 'log_vr_bias'),
  tmb_map = list(),
  template_fp = template_fp,
  optimization_method = 'L-BFGS-B',
  model_name = 'Mexico BH-VR NMR model',
  parallel_model = TRUE, verbose = TRUE
)
sdrep <- sdreport(model_fit$obj, bias.correct = TRUE, getJointPrecision = TRUE)
prepped_data[, idx_loc := uid ]

all_draws <- generate_draws(
  tmb_sdreport = sdrep,
  data_template = prepped_data,
  num_draws = 1000,
  covariate_names = use_covs
)

param_summs <- cbind(
  data.table(parameter = all_draws$param_names),
  as.data.table(summarize_draws(all_draws$param_draws))
)

pred_summs <- cbind(
  prepped_data,
  as.data.table(summarize_draws(all_draws$predictive_draws))
)
for(pname in c('Z_loc','log_vr_bias')){
  for(val in c('mean','lower','upper')){
    pred_summs[[paste(pname, val, sep='_')]] <- (
      param_summs[[val]][all_draws$param_names==pname]
    )
  }
}

fe_summs <- cbind(
  data.table(covariate = use_covs),
  as.data.table(param_summs[all_draws$param_names=='beta_covs'])
)


## Write outputs ------------------------------------------------------------------------>

# Summaries
fwrite(all_draws$param_draws, file = file.path(model_dir, 'param_draws.csv'), row.names=T)
fwrite(param_summs, file = file.path(model_dir, 'param_summs.csv'))
fwrite(all_draws$predictive_draws, file = file.path(model_dir,'pred_draws.csv'), row.names=T)
fwrite(pred_summs, file = file.path(model_dir,'pred_summs.csv'))
fwrite(fe_summs, file = file.path(model_dir,'fe_summs.csv'))

# Model objects
saveRDS(model_fit, file = file.path(model_dir, 'model_fit.RDS'))
saveRDS(sdrep, file = file.path(model_dir, 'sdrep.RDS'))
saveRDS(vr_bias_thresholds, file = file.path(model_dir, 'vr_bias_thresholds.RDS'))
if(sim) saveRDS(sim_args, file=file.path(model_dir, 'sim_args.RDS'))
