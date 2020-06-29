//
// scMET: Hierarchical Beta Binomial regression on both the mean (proportions) and overdispersion parameters.
// 1. Here we assume \mu ~ LogitNormal(w_{m}^{T}X, s_{m}); we put normal prior on w_{m} and
//    X are the covariates, e.g. # of CpGs.
// 2. We also assume \gamma ~ LogitNormal(w_0 * 1 + w_{g}^{T} g(\mu)); we put normal prior on c(w_0, w_{g}),
//    g() is a basis function (e.g. RBF), that transforms the \mu to allow us to capture non-linear
//    relationships between mean and overdispersion.
//
// NOTE: This is a general model, which contains linear regression and non-regression models as special cases.
//       1. If we assume X = [1, ..., 1] then we do not have regression on the mean,
//          so the prior is the population average across all features. However we still correct for the
//          mean-overdispersion relationship.
//       2. If we do not include any basis function g(), by setting L = 1, then we do not capture the
//          mean-overdispersion trend, however we still capture the effect of covariates on the mean.
//       3. If we do not include any of the above, it reduces to a simple BB model without regression.
//
functions{
  matrix rbf_H(int L, vector X, real c) {
    matrix[num_elements(X), L] H; // Design matrix
    if (L > 1) {
      vector[L - 1] ms; // Mean locations
      real h; // Scaling parameter
      // Compute rbf centres (assume equally spaced)
      // The equally spaced ms are obtained by ms[i] = i * (rmax - rmin)/M + rmin
      for (l in 1:(L - 1)) {
        ms[l] = l * ((max(X) - min(X)) / L ) + min(X);
      }
      // Compute scaling parameter
      h = (ms[2] - ms[1]) * c;
      // Populate design matrix
      for (l in 1:L) {
        if (l == 1) {
          H[, l] = rep_vector(1.0, num_elements(X));
        } else {
          for (i in 1:num_elements(X)) {
            H[i, l] = exp( -0.5 * pow((X[i] - ms[l - 1]) / h, 2));
          }
        }
      }
    } else {
      H[, 1] = rep_vector(1.0, num_elements(X));
    }
    return(H);
  }

  matrix poly_H(int L, vector X) {
    matrix[num_elements(X), L] H; // Design matrix
    for (l in 1:L) {
      if (l == 1) {
        H[, l] = rep_vector(1.0, num_elements(X));
      } else if (l == 2) {
        H[, l] = X;
      } else {
        for (i in 1:num_elements(X)) {
          H[i, l] = pow(X[i], l - 1);
        }
      }
    }
    return(H);
  }
}

data {
  int<lower=0> N;         // Number of data (all features and cells)
	int<lower=0> J;         // # of features
	int<lower=0> N_X;       // # of covariates for mean parameter
	int<lower=0> L;         // Number of basis functions (including intercept term)
	int<lower=0> n[N];      // total number of trials (long form)
	int<lower=0> y[N];      // number of successes (long form)
	int<lower=0> C[J];      // # of cells with obserations in each feature
	real<lower=0> rbf_c;    // Scaling parameter for spatial variance of RBFs.
	matrix[J, N_X] X;       // Covariates for each feature J
	vector[N_X] m_wmu;      // Normal mean vector prior on w_\mu
	real<lower=0> s_wmu;    // Normal std prior on w_\mu
	real<lower=0> s_mu;     // LogitNormal std prior on proportions \mu
	vector[L] m_wgamma;     // Normal mean vector prior on w_\mu
	real<lower=0> s_wgamma; // Normal std prior on w_\mu
	real<lower=0> a_sgamma; // Inv Gamma shape prior on population dispersion \s_gamma
	real<lower=0> b_sgamma; // Inv Gamma rate prior on population dispersion \s_gamma
}

// Parameters to compute posterior, reparametrization of BB in terms of mean
// and overdispersion.
parameters {
  vector<lower=-10, upper=10>[N_X] w_mu;      // Regression coefficients of \mu
  vector<lower=-10, upper=10>[L] w_gamma;     // Regression coefficients of \gamma
  vector<lower=-20, upper=20>[J] logit_mu;    // Mean of BB
  vector<lower=-20, upper=20>[J] logit_gamma; // Dispersion of BB
  real<lower=0> s_gamma;                      // Population std of LogitNormal
}

// Original parameters of the Beta distribution
transformed parameters {
  vector<lower=1e-15, upper=(1-1e-15)>[J] mu = inv_logit(logit_mu);
  vector<lower=1e-15, upper=(1-1e-15)>[J] gamma = inv_logit(logit_gamma);
  vector<lower=-40, upper=40>[J] f_mu = X * w_mu;
  vector<lower=-40, upper=40>[J] f_gamma;
	// Non-linear regression  for mean-overdispersion relationship
	f_gamma = rbf_H(L, mu, rbf_c) * w_gamma;
}

// Joint distribution for Beta Binomial model
model {
  int pos; // Counter to perform iteration only for cells with coverage
  pos = 1;
  s_gamma ~ inv_gamma(a_sgamma, b_sgamma); // Hyperprior on the population overdisperion std
  w_mu ~ normal(m_wmu, s_wmu);             // Regression coefficients w_{\mu} follow MN distribution
  w_gamma ~ normal(m_wgamma, s_wgamma);    // Regression coefficients w_{\gamma} follow MN distribution
  logit_mu ~ normal(f_mu, s_mu);           // \mu follows LogitNormal distribution
	logit_gamma ~ normal(f_gamma, s_gamma);  // \gamma follows LogitNormal distribution
  // Beta-Binomial formulation, i.e. integrating out theta
  for (j in 1:J) {
    //segment(y, pos, C[j]) ~ beta_binomial(segment(n, pos, C[j]), alpha[j], beta[j]); // likelihood
    segment(y, pos, C[j]) ~ beta_binomial(segment(n, pos, C[j]),
                                          mu[j]/gamma[j] - mu[j],
                                          (1 - mu[j])/gamma[j] + mu[j] - 1); // likelihood
    pos = pos + C[j];     // Increment counter
  }
}

// Generated quantities
generated quantities {
  vector[J] epsilon = logit_gamma - f_gamma;
}
