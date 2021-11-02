// ///////////////////////////////////////////////////////////////////////////////////////
//
// Author: Nat Henry
// Created: 19 August 2021
// Purpose: TMB objective function for Mexico two-source neonatal mortality model
//
// ///////////////////////////////////////////////////////////////////////////////////////


#include <TMB.hpp>
using namespace density;
using Eigen::SparseMatrix;


// HELPER FUNCTIONS --------------------------------------------------------------------->

// Function for rescaling a precision matrix to have standard deviation sigma
//
// Parameter Q: Unscaled precision matrix
// Parameter sigma: Standard deviation to scale to
//
template<class Type>
SparseMatrix<Type> scale_precision(SparseMatrix<Type> Q, Type sigma){
  SparseMatrix<Type> Q_scaled = Q / (sigma * sigma);
  return Q_scaled;
}

// Function to create an IID precision matrix (AKA a scaled identity matrix)
//
// Parameter dim: Number of rows (and columns) for the precision matrix
// Parameter sigma: Standard deviation of the iid process
//
template<class Type>
SparseMatrix<Type> iid_precision(int dim, Type sigma = 1.0){
  SparseMatrix<Type> I(dim, dim);
  for(int ii=0; ii < dim; ii++){
    I.insert(ii, ii) = 1.0;
  }
  SparseMatrix<Type> I_scaled = scale_precision(I, sigma);
  return I_scaled;
}

// Function for preparing a precision matrix corresponding to a BYM2 spatial model,
//   based on a scaled ICAR precision matrix. For more details, see:
//   Riebler et al. (2016). An intuitive Bayesian sptial model for disease mapping that
//     accounts for scaling. Statistical methods in medical research, 25(4):1145-65.
//
// Parameter Q_icar: A precision matrix corresponding to an intrinsic correlated
//   autoregressive (ICAR) model in space, scaled to have generalized variance 1
// Parameter phi: A mixing parameter indicating the relative contribution of spatial and
//   IID variation, strictly between 0 and 1.
// Parameter sigma: Standard deviation of the LCAR process
//
template<class Type>
SparseMatrix<Type> bym2_precision(SparseMatrix<Type> Q_icar, Type phi, Type sigma = 1.0){
  SparseMatrix<Type> I = iid_precision(Q_icar.rows(), Type(1.0));
  SparseMatrix<Type> Q_bym2 = phi * Q_icar + (1 - phi) * I;
  SparseMatrix<Type> Q_bym2_scaled = scale_precision(Q_bym2, sigma);
  return Q_bym2_scaled;
}


// FUNCTION: Get log-density of penalized complexity prior based on log-precision
//
// Parameter logtau: Log of the precision term for a distribution. The P.C. prior
//   sets a probability that the associated standard deviation sigma,
//   sigma = exp(logtau * -1/2), exceeds some threshold u. In other words,
//   Pr(sigma > u) = alpha.
// Parameter u: Threshold to test against sigma
// Parameter alpha: Probability that sigma exceeds the threshold
//
// Returns the log-density, which can be subtracted from the joint negative log likelihood
//
template<class Type>
Type pc_dprec(Type logtau, Type u, Type alpha){
  Type lambda = Type(-1.0) * log(alpha) / u;
  Type logdens = log(lambda / Type(2.0)) - (lambda * exp(logtau / Type(-2.0))) - logtau / Type(2.0);
  return logdens;
}


// OBJECTIVE FUNCTION ------------------------------------------------------------------->

template<class Type>
Type objective_function<Type>::operator() () {

  // INPUT DATA ------------------------------------------------------------------------->

  // OPTION: Holdout number
  // Any observation where `idx_holdout_vr` or `idx_holdout_bh` is equal to `holdout`
  //   will have the corresponding data type excluded from this model fit
  // Holdouts are 1-indexed, so if holdout==0 all data will be used
  DATA_INTEGER(holdout);

  // Input mortality data and denominators (all have the same length)
  DATA_VECTOR(vr_births);
  DATA_VECTOR(vr_deaths);
  DATA_VECTOR(bh_births);
  DATA_VECTOR(bh_deaths);

  // Covariate matrix: dimensions (n observations) by (m covariates including intercept)
  DATA_MATRIX(X_ij);

  // Indices
  DATA_IVECTOR(idx_loc);     // Index for the location of each observation
  DATA_IVECTOR(idx_vr_bias); // Index determining shape of VR bias parameter
  DATA_IVECTOR(idx_holdout_vr); // Holdout index for each VR observation
  DATA_IVECTOR(idx_holdout_bh); // Holdout index for each BH observation

  // Hyperpriors - thresholds for VR bias SD terms (used in penalized complexity priors)
  DATA_VECTOR(bias_thresh_vec);

  // Precision matrix for an ICAR spatial model
  DATA_SPARSE_MATRIX(Q_icar);


  // INPUT PARAMETERS ------------------------------------------------------------------->

  // Fixed effects
  PARAMETER_VECTOR(beta_covs); // Vector of fixed effects on covariates

  // Log precision of the spatial random effects
  PARAMETER(log_tau_loc);
  // Mixing parameter controlling spatial vs. nonspatial correlation by municipality
  PARAMETER(logit_phi_loc);

  // Vector of VR bias standard deviations (centered around zero)
  PARAMETER_VECTOR(vr_bias_logprec_vec);

  // Parameter for BH bias (a single scalar term)
  PARAMETER(log_bh_bias);
  // Correlated random effect surface for mortality
  PARAMETER_VECTOR(Z_loc);
  // Bias terms for each VR observation
  PARAMETER_VECTOR(log_vr_bias);

  // VR births, VR deaths, BH births, and BH deaths are all vectors of the same length
  int num_obs = vr_births.size();
  // Vector of VR bias standard deviations (from precision)
  vector<Type> vr_bias_sd_vec = exp(vr_bias_logprec_vec / Type(-2.0));
  // Transformed BH and VR bias terms
  Type bh_bias = exp(log_bh_bias);
  vector<Type> vr_bias = exp(log_vr_bias);


  // Instantiate joint negative log-likelihood (JNLL) ----------------------------------->

  parallel_accumulator<Type> jnll(this);

  // JNLL CONTRIBUTION FROM PRIORS ------------------------------------------------------>

  // Gamma(1, 1E3) hyperprior for tau precision prior
  jnll -= dlgamma(exp(log_tau_loc), Type(1.0), Type(1000.0), true);
  // Spatial effect = BYM2 (scaled CAR) model using municipal neighborhood structure
  SparseMatrix<Type> Q_loc = bym2_precision(
    Q_icar,
    invlogit(logit_phi_loc),
    exp(log_tau_loc * Type(-0.5))
  );
  jnll += GMRF(Q_loc)(Z_loc);

  // Penalized complexity priors on standard deviations for each municipality grouping,
  //  parameterized through the log precision
  for(int i = 0; i < vr_bias_logprec_vec.size(); i++){
    // Arguments: Log-precision, threshold, probability of exceeding
    jnll -= pc_dprec(vr_bias_logprec_vec(i), bias_thresh_vec(i), Type(0.01));
  }

  // N(mean=0, sd=3) prior for fixed effects
  // Skip the intercept (index 0)
  for(int cov_j = 1; cov_j < beta_covs.size(); cov_j++){
    jnll -= dnorm(beta_covs(cov_j), Type(0.0), Type(3.0), true);
  }

  // Apply priors to VR bias terms
  for(int obs_i = 0; obs_i < num_obs; obs_i++){
    jnll -= dnorm(log_vr_bias(obs_i), Type(0.0), vr_bias_sd_vec(idx_vr_bias(obs_i)), true);
  }


  // JNLL CONTRIBUTION FROM DATA -------------------------------------------------------->

  // Logit(true mortality) = fixed effects + random effects
  vector<Type> fix_effs = X_ij * beta_covs.matrix();

  for(int i=0; i < num_obs; i++){
    // Incorporate BH data
    if((idx_holdout_bh(i) != holdout) && (bh_births(i) > 0.)){
      jnll -= dpois(
        bh_deaths(i),
        bh_births(i) * exp(fix_effs(i) + Z_loc(idx_loc(i))) * bh_bias,
        true
      );
    }
    // Incorporate VR data
    if((idx_holdout_vr(i) != holdout) && (vr_births(i) > 0.)){
      jnll -= dpois(
        vr_deaths(i),
        vr_births(i) * exp(fix_effs(i) + Z_loc(idx_loc(i))) * vr_bias(i),
        true
      );
    }
  }

  // RETURN JNLL ------------------------------------------------------------------------>

  return jnll;

} // END objective function
