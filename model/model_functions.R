## MODEL RUNNING AND POSTESTIMATION FUNCTIONS ------------------------------------------->


#' Generate an ICAR precision matrix based on an adjacency matrix
#'
#' @description Generate a precision matrix for the intrinsic correlated autoregressive
#'  (ICAR) model specification, a special case of the correlated autoregressive (CAR)
#'  class of Markov random field models. This precision matrix is usually denoted as "Q".
#'
#' @details The precision matrix is fully specified by the adjacency weights, matrix W,
#'   defined as W = {w_ij} where w_ij is 1 if i and j are neighbors, and 0 otherwise. The
#'   precision matrix Q is defined as Q = D_w - W, where D_w is a diagonal matrix with
#'   each diagonal term d_ii equal to the sum of row i in W.
#'
#'   Note that the ICAR model is improper, in that the conditional distributions
#'   specified by the precision matrix do not determine a full joint distribution that
#'   integrates to 1; in other words, the precision matrix Q is not invertible. The ICAR
#'   precision matrix can still be used as a prior in a hierarchical model.
#'
#'   This function includes optional argument `scale_variance`. If set to `TRUE` (the
#'   default), the function will rescale the precision matrix to have a generalized
#'   variance of 1, which may aid in prior specifications that are comparable across
#'   areal spatial models with different geometries.
#'
#'   For more details, see:
#'   Banerjee, Carlin, and Gelfand (2015). Hierarchical Modeling and Analysis for Spatial
#'     Data, 2nd Edition. Section 6.4.3.3: CAR models and their difficulties.
#'   Riebler et al. (2016). An intuitive Bayesian sptial model for disease mapping that
#'     accounts for scaling. Statistical methods in medical research, 25(4):1145-65.
#'
#' @param W Adjacency matrix, with w_ij = w_ji = 1 if areal units i and j are neighbors,
#'   and zero otherwise. See function details for more information
#' @param scale_variance [default TRUE] Should the precision matrix be rescaled so that
#'  the generalized variance is equal to 1? Setting to TRUE may help with prior
#'  specification.
#'
#' @return Sparse ICAR precision matrix Q. See function details for more information.
#'
#' @import Matrix INLA
#' @export
icar_precision_from_adjacency <- function(W, scale_variance = TRUE){
  # Generate and return sparse precision matrix
  Q <- Matrix::Diagonal(n = nrow(W), x = Matrix::rowSums(W)) - W
  if(scale_variance){
    # Scale model to have generalized variance of 1
    constraint_matrix <- matrix(1, nrow = 1, ncol = ncol(Q))
    Q <- INLA::inla.scale.model(Q, constr = list(A = constraint_matrix, e = 0))
  }
  return(Q)
}


#' Get Maximum a Priori (MAP) parameter estimates
#'
#' @description Find the Maximum a Priori (MAP) estimates for fixed effect
#'   parameter values, including age-specific fixed effects, in a simplified
#'   version of a GLM. Starting the full TMB model with fixed effect starting
#'   values set at the MAP has been shown to improve performance and run time.
#'
#' @param in_data Input data.table, including only the data used to train the
#'    model. Fields must include a numerator, a denominator, fields named for
#'    all covariates, and (optionally) an age ID field
#' @param events_field Field name for the number of events (eg deaths)
#' @param exposure_field Field name for the relative exposure (eg. pop-time)
#' @param covar_names Names of all covariates, including 'intercept' if desired
#' @param distribution_family Name of the distribution to to fit the GLM with,
#'   which also determines the link function to be used. Should be a valid
#'   argument to pass to `stats::glm(family = distribution_family)`
#' @param grouping_field [optional] Field that groups observations by ID
#'
#' @return Named list with two items:
#'     - "glm_fit": Full model fit for the GLM
#'     - "fixed_effects_map": Maximum a priori estimates for covariate fixed
#'          effects, organized as a named numeric vector
#'     - "fixed_effects_grouping": Maximum a priori estimates for grouped fixed
#'          effects, organized as a vector of length(num groups). This list item
#'          is NULL if the argument `grouping_field` was not specified
#'
#' @import data.table glue stats
#' @export
find_glm_map_parameter_estimates <- function(
  in_data, events_field, exposure_field, covar_names, distribution_family,
  grouping_field = NULL
){
  # Ensure that all columns are available in input data
  reqd <- c(events_field, exposure_field, covar_names, grouping_field)
  missing_fields <- setdiff(reqd, colnames(in_data))
  if(length(missing_fields) > 0){
    stop("MAP input data missing fields: ", paste(missing_fields, collapse=', '))
  }
  in_data <- na.omit(in_data, cols = reqd)

  # Add group-based fixed effects, if specified
  grp_cols <- c()
  if(!is.null(grouping_field)){
    if(in_data[, uniqueN(get(grouping_field)) ] > 1){
      grp_vals <- sort(unique(in_data[[grouping_field]]))
      # The first group is set as 'default' - others vary with a fixed effect
      for(grp_val in grp_vals[2:length(grp_vals)]){
        in_data[, paste0('grp',grp_val) := 0 ]
        in_data[ get(grouping_field) == grp_val, paste0('grp',grp_val) := 1 ]
      }
      grp_cols <- paste0('grp', grp_vals[2:length(grp_vals)])
    }
  }
  # Get a field representing successes as a proportion of exposure
  in_data[, rate_success := get(events_field) / get(exposure_field) ]

  # Set up formula with an offset based on link function, then run the GLM
  formula_char <- glue::glue(
    "rate_success ~ 0 + {paste(c(covar_names, grp_cols), collapse = ' + ')}"
  )
  .env <- environment()
  formula_parsed <- as.formula(formula_char, env = .env)
  glm_fit <- stats::glm(
    formula_parsed, data = in_data, family = distribution_family,
    weights = in_data[[exposure_field]]
  )

  # Return the full GLM fit, the covariate fixed effects, and (optionally) the
  #  grouped fixed effects
  covs_map <- glm_fit$coefficients[ covar_names ]
  if(length(grp_cols) > 0){
    grps_map <- c(0, glm_fit$coefficients[ grp_cols ])
    names(grps_map) <- paste0('grp',grp_vals)
  } else {
    grps_map <- NULL
  }
  return(list(
    glm_fit = glm_fit,
    fixed_effects_map = covs_map,
    fixed_effects_grouping = grps_map
  ))
}

#' Set up and run TMB
#'
#' @description Generic TMB model run handler. Sets up the ADFun object, applies
#'   model speedups and fixes as specified, and optimizes using `nlminb`. This
#'   is meant to be a helper function run by more specific model functions.
#'
#' @param tmb_data_stack List containing all data inputs expected in the TMB
#'   CPP file
#' @param params_list List containing all parameters expected in the TMB CPP file
#' @param tmb_random Character vector containing all random effects that will be
#'   optimized in the inner optimizer
#' @param tmb_map Named list containing parameters that will be treated in a
#'   particular way by the optimizer
#' @param template_fp Filepath containing the TMB objective function
#' @param tmb_outer_maxsteps Max number of steps taken by the outer optimizer
#' @param tmb_inner_maxsteps Max number of steps taken by the inner optimizer
#'   in a single outer optimizer step
#' @param parallel_model [bool, default FALSE] Is the model implemented in parallel? If
#'   TRUE, opens multiple OMP threads before fitting
#' @param optimization_method [char, default 'nlminb'] Outer optimization method to use
#'   for fitting, implemented in the optimx library. Recommended options include 'nlminb'
#'   and 'L-BFGS-B'
#' @param model_name [char, default "model"] name of the model
#' @param verbose [boolean, default FALSE] Should this function return logging
#'   information about the stage of model fitting, including the outer optizimer
#'   sampling? This will be somewhat long (>100 lines)
#' @param inner_verbose [boolean, default FALSE] Should this function return
#'   logging information about inner optimizer sampling? This can be useful for
#'   debugging, but will show very verbose (>10k lines) output.
#'
#' @return list of two objects: obj (ADFunction object), and opt (optimized
#'   nlminb object)
#'
#' @import TMB glue tictoc optimx
#' @export
setup_run_tmb <- function(
  tmb_data_stack, params_list, tmb_random, tmb_map, template_fp, tmb_outer_maxsteps = 1E3,
  tmb_inner_maxsteps = 1E3, parallel_model = FALSE, optimization_method = 'nlminb',
  model_name="model", verbose=FALSE, inner_verbose=FALSE
){
  # Helper function to send a message only if verbose
  vbmsg <- function(x) if(verbose) message(x)
  # Setup
  vbmsg(paste0(c("\n",rep("*",nchar(model_name)+14)),collapse=''))
  vbmsg(glue::glue("***  {model_name} RUN  ***"))
  vbmsg(paste0(c(rep("*",nchar(model_name)+14),"\n"),collapse=''))

  if(parallel_model){
    # Set up openmp threads
    threads <- system('echo $OMP_NUM_THREADS', intern = TRUE)
    if(threads != '') {
      vbmsg(sprintf('Detected %s threads in OMP environmental variable.',threads))
      openmp(as.numeric(threads))
    } else {
      vbmsg("Did not detect environmental OMP variable, defaulting to 2 cores. \n
             You can set this using OMP_NUM_THREADS.")
      openmp(2)
    }
  }

  # Compile TMB C++ template
  TMB::compile(template_fp)
  current_dir <- getwd()
  setwd(dirname(template_fp))
  compiled_path <- tools::file_path_sans_ext(basename(template_fp))
  dyn.load(TMB::dynlib(compiled_path))

  # Make Autodiff function
  vbmsg("Constructing ADFunction...")
  tictoc::tic("  Making Model ADFun")
  obj <- TMB::MakeADFun(
    data = tmb_data_stack, parameters = params_list, random = tmb_random,
    map = tmb_map, DLL = compiled_path, silent = !inner_verbose,
    inner.control = list(trace=inner_verbose, tol=1E-11)
  )
  obj$env$tracemgc <- as.integer(verbose)
  tictoc::toc()

  # Optimize using the specified outer optimizer, implemented in optimx
  tictoc::tic("  Optimization")
  vbmsg(glue("\n** OPTIMIZING USING METHOD {optimization_method} **"))
  opt <- optimx::optimx(
    par = obj$par, fn = function(x) as.numeric(obj$fn(x)), gr = obj$gr,
    method = optimization_method,
    itnmax = tmb_outer_maxsteps,
    hessian = FALSE,
    control = list(
      trace = as.integer(verbose), follow.on = TRUE,
      dowarn = as.integer(verbose), maxit = tmb_inner_maxsteps,
      factr = 1E-10
    )
  )
  conv_code <- opt$convcode
  vbmsg(glue::glue(
    "{model_name} optimization finished with convergence code {conv_code}.\n"
  ))
  # Clean up
  tictoc::toc()
  setwd(current_dir)
  vbmsg(glue::glue("*** {model_name} RUN COMPLETE **************\n\n"))
  return(list(obj=obj, opt=opt))
}


#' Take multivariate normal draws given a mean vector and precision matrix
#'
#' @param mu vector of parameter means
#' @param prec joint precision matrix
#' @param n.sims number of draws
#'
#' @return length(mu) by n.sims matrix of parameter draws
#'
#' @import Matrix
#' @export
rmvnorm_prec <- function(mu, prec, n.sims) {
  z = matrix(rnorm(length(mu) * n.sims), ncol=n.sims)
  L_inv = Matrix::Cholesky(prec, super=TRUE)
  return(mu + solve(as(L_inv, 'pMatrix'), solve(t(as.matrix(as(L_inv, 'Matrix'))), z)))
}


#' Create post-estimation predictive draws
#'
#' @description Given the outputs from a fitted TMB model object, create
#'   an object with posterior predictive draws for all groupings specified by a
#'   template data.table
#'
#' @param tmb_sdreport output of `TMB::sdreport()` on the fitted model object.
#'   Should include a joint precision matrix (by specifying
#'   `getJointPrecision = TRUE` in the call to `sdreport()`). This object will be
#'   parsed to check for fixed effects, random effects, and the Fourier time
#'   series terms.
#' @param data_template Prepped data with covariates, in random effects order
#' @param num_draws [int] How many posterior predictive samples to take?
#' @param covariate_names [char] All covariate field names, including 'intercept'
#'
#' @return A named list with three items:
#'    - 'param_names': Vector of parameter names in the order they have been
#'         extracted
#'    - 'param_draws': Matrix of parameter draws
#'    - 'pred_draws': Matrix of mortality predictive draws, taken at the
#'         observation points specified in the `template_dt`
#'
#' @import data.table
#' @export
generate_draws <- function(
  tmb_sdreport, data_template, num_draws, covariate_names
){
  # Copy input data
  templ <- data.table::copy(data_template)

  # Get parameter names
  mu <- c(tmb_sdreport$par.fixed, tmb_sdreport$par.random)
  parnames <- names(mu)

  ## Input data checks
  # Check that joint precision matrix was created
  if(!"jointPrecision" %in% names(tmb_sdreport)) stop("Missing joint precision matrix")
  # Check that all the covariate names are there
  missing_covs <- setdiff(covariate_names, names(templ))
  if(length(missing_covs) > 0){
    stop("Missing covariates: ", paste0(missing_covs, collapse=', '))
  }
  if(length(covariate_names) != sum(parnames == 'beta_covs')){
    stop("Wrong number of covariates in model fit")
  }
  # Check that template data.table has all required columns
  template_req_cols <- c('idx_loc')
  missing_templ_cols <- setdiff(template_req_cols, names(templ))
  if(length(missing_templ_cols) > 0){
    stop("Missing columns: ", paste0(missing_templ_cols, collapse=', '))
  }

  ## Get parameter draws
  message(sprintf(" - Generating %i parameter draws...", num_draws))
  prec_mat <- tmb_sdreport$jointPrecision
  if(any(colnames(prec_mat) != parnames )) stop("Issue with parameter ordering")
  param_draws <- rmvnorm_prec(
    mu = mu,
    prec = prec_mat,
    n.sims = num_draws
  )
  rownames(param_draws) <- parnames

  ## Generate predictive draws from parameter draws
  # Sort by id fields beforehand to add random effects more easily
  templ[, row_id := .I ]
  templ <- templ[order(idx_loc)]
  templ[, sorted_row_id := .I ]

  # Prediction = inverse.logit( Covar FEs + age FEs + REs + sometimes seasonality )
  # Covariate fixed effects
  fes <- as.matrix(templ[, ..covariate_names]) %*% param_draws[parnames=='beta_covs', ]
  res <- param_draws[parnames=='Z_loc', ]
  preds <- plogis(fes + res)

  # Reorder by the original row ordering
  templ <- templ[order(row_id)]
  preds <- preds[templ$sorted_row_id, ]
  # Return parameter draws and predictive draws
  return(list(
    param_names = parnames,
    param_draws = param_draws,
    predictive_draws = preds
  ))
}


#' Summarize predictive draws
#'
#' @description Summarize the mean and select quantiles of a matrix of posterior
#'   draws, where draws are stored in columns
#'
#' @param draws [matrix] matrix of dimensions (num obs) by (num draws)
#'
#' @return data.table with columns 'mean','median','upper','lower' (of 95% UI)
#'
#' @import data.table matrixStats
#' @export
summarize_draws <- function(draws){
  if('data.table' %in% class(draws)) draws <- as.matrix(draws)
  summs <- cbind(
    rowMeans(draws), matrixStats::rowQuantiles(draws, probs=c(0.5, 0.025, 0.975))
  )
  colnames(summs) <- c('mean','median','lower','upper')
  return(as.data.table(summs))
}

