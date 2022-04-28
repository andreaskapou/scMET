#' @title Simulate methylation data from scMET.
#'
#' @description General function for simulating datasets with diverse proprties.
#'   This for instance include, adding covariates X that explain differences in
#'   mean methylation levels. Or also defining the trend for the mean -
#'   overdispersion relationship.
#'
#' @param N_feat Total number of features (genomics regions).
#' @param N_cells Maximum number of cells.
#' @param N_cpgs Maximum number of CpGs per cell and feature.
#' @param L Total number of radial basis functions (RBFs) to fit the
#'   mean-overdispersion trend. For L = 1, this reduces to a model that does not
#'   correct for the mean-overdispersion relationship.
#' @param X Covariates which might explain variability in mean (methylation). If
#'   X = NULL, a 2-dim matrix will be generated, first column containing
#'   intercept term (all values = 1), and second colunn random generated
#'   covariates.
#' @param w_mu Regression coefficients for covariates X. Should match number of
#'   columns of X.
#' @param s_mu Standard deviation for mean parameter `mu`.
#' @param w_gamma Regression coefficients of the basis functions. Should match
#'   the value of L. If NULL, random coefficients will be generated.
#' @param s_gamma Standard deviation of dispersion parameter `gamma`.
#' @param rbf_c Scale parameter for empirically computing the variance of the
#'   RBFs.
#' @param cells_range Range (betwen 0 and 1) to randomly (sub)sample the number
#'   of cells per feature.
#' @param cpgs_range Range (betwen 0 and 1) to randomly (sub)sample the number
#'   of CpGs per cell and feature.
#'
#' @importFrom logitnorm rlogitnorm
#' @return A simulated dataset and additional information for reproducibility
#'   purposes.
#'
#' @examples
#' sim <- scmet_simulate(N_feat = 150, N_cells = 50, N_cpgs = 15, L = 4)
#'
#' @export
#'
scmet_simulate <- function(N_feat = 100, N_cells = 50, N_cpgs = 15, L = 4,
                           X = NULL, w_mu = c(-0.5, -1.5), s_mu = 1,
                           w_gamma = NULL, s_gamma = 0.3, rbf_c = 1,
                           cells_range = c(0.4, 0.8), cpgs_range = c(0.4, 0.8)) {
  # So RMD check does not complain
  Feature <- NULL

  # Parameter checks
  if (is.null(L)) { L <- 4 }
  if (is.null(X) & is.null(w_mu)) { w_mu <- c(-0.5, -1.5)}
  # Initialize w_gamma
  w_gamma <- .init_w_gamma(w_gamma = w_gamma, L = L)

  # Generate total number of CpGs per feature and cell
  cpgs_list <- .generate_cpgs(N_feat = N_feat, N_cells = N_cells,
                              N_cpgs = N_cpgs, cells_range = cells_range,
                              cpgs_range = cpgs_range)
  # Generate mean methylation levels \mu
  tmp <- .generate_means(N_feat = N_feat, N_cpgs = N_cpgs, X = X, w_mu = w_mu,
                         s_mu = s_mu)
  mu <- tmp$mu
  X <- tmp$X
  # Generate overdispersion parameters \gamma
  gamma <- .generate_overdisp(N_feat = N_feat, L = L, mu = mu,
                              w_gamma = w_gamma, s_gamma = s_gamma, rbf_c = rbf_c)

  # Generate number of methylated CpGs from Beta Binomial
  met_cpgs_list <- lapply(X = seq_len(N_feat), function(n)
    VGAM::rbetabinom(length(cpgs_list[[n]]), cpgs_list[[n]],
                     prob = mu[n], rho = gamma[n]))

  # Create observed data object, which will be used as input to the model
  Y <- .generate_ys(N_feat = N_feat, N_cells = N_cells, cpgs_list = cpgs_list,
                    met_cpgs_list = met_cpgs_list)
  # Set rownames of covariates to feature names
  rownames(X) <- unique(Y$Feature)
  # Store parameters
  theta <- data.frame(mu = mu, gamma = gamma, row.names = unique(Y$Feature))

  ## Order rows by Feature name for all three objects
  Y <- Y[order(Feature), , drop = FALSE]
  X <- as.matrix(X[order(rownames(X)), , drop = FALSE])
  theta <- theta[order(rownames(theta)), , drop = FALSE]

  theta_priors <- list(w_mu = w_mu, w_gamma = w_gamma, s_mu = s_mu,
                       s_gamma = s_gamma)
  opts <- list(N_feat = N_feat, N_cells = N_cells, N_cpgs = N_cpgs, L = L,
               rbf_c = rbf_c, cells_range = cells_range,
               cpgs_range = cpgs_range)
  obj <- structure(list(Y = Y, X = X, theta_true = theta,
                        theta_priors_true = theta_priors, opts = opts),
                   class = "scmet_simulate")
  return(obj)
}


#' @title Simulate differential methylation data from scMET.
#'
#' @description General function for simulating two methylation datasets for
#'   performing differential methylation analysis. Differential analysis can be
#'   either performed in detecting changes in mean or variability of methylation
#'   patterns between the two groups. Similar to \code{\link{scmet_simulate}},
#'   the function allows inclusion of covariates X that explain differences in
#'   mean methylation levels. Or also defining the trend for the mean -
#'   overdispersion relationship.
#'
#' @param diff_feat_prcg_mu Percentage of features (betwen 0 and 1) that show
#'   differential mean methylation between the two groups.
#' @param diff_feat_prcg_gamma Percentage of features (betwen 0 and 1) that show
#'   differential variability between the two groups.
#' @param OR_change_mu Effect size change (in terms of odds ratio) of mean
#'   methylation between the two groups.
#' @param OR_change_gamma Effect size change (in terms of odds ratio) of
#'   methylation variability between the two groups.
#' @inheritParams scmet_simulate
#'
#' @return Methylation data from two cell populations/conditions.
#'
#' @examples
#' sim_diff <- scmet_simulate_diff(N_feat = 150, N_cells = 100, N_cpgs = 15, L = 4)
#'
#' @export
#'
scmet_simulate_diff <- function(N_feat = 100, N_cells = 50, N_cpgs = 15, L = 4,
                                diff_feat_prcg_mu = 0,
                                diff_feat_prcg_gamma = 0.2,
                                OR_change_mu = 3, OR_change_gamma = 3,
                                X = NULL, w_mu = c(-.5, -1.5), s_mu = 1,
                                w_gamma = NULL, s_gamma = 0.3, rbf_c = 1,
                                cells_range = c(0.4, 0.8),
                                cpgs_range = c(0.4, 0.8)) {

  # Parameter checks
  if (is.null(L)) { L <- 4 }
  if (is.null(X) & is.null(w_mu)) { w_mu <- c(-0.5, -1.5)}
  assertthat::assert_that(diff_feat_prcg_mu >= 0 & diff_feat_prcg_mu <= 1)
  assertthat::assert_that(diff_feat_prcg_gamma >= 0 & diff_feat_prcg_gamma <= 1)
  assertthat::assert_that(OR_change_mu > 0 )
  assertthat::assert_that(OR_change_gamma > 0)

  # Initialize w_gamma
  w_gamma <- .init_w_gamma(w_gamma = w_gamma, L = L)

  # Generate data from group A
  cat("Simulating from group A\n")
  sim_dt_A <- scmet_simulate(N_feat = N_feat, N_cells = N_cells, N_cpgs = N_cpgs,
                             L = L, X = X, w_mu = w_mu, s_mu = s_mu,
                             w_gamma = w_gamma, s_gamma = s_gamma, rbf_c = rbf_c,
                             cells_range = cells_range, cpgs_range = cpgs_range)

  cat("Simulating from group B\n")
  # Generate total number of CpGs per feature and cell for group B
  cpgs_list_B <- .generate_cpgs(N_feat = N_feat, N_cells = N_cells, N_cpgs = N_cpgs,
                                cells_range = cells_range, cpgs_range = cpgs_range)
  # Extract mean parameters from group A
  mu_B <- sim_dt_A$theta_true$mu
  # Extract dispersion parameters from group A
  gamma_B <- sim_dt_A$theta_true$gamma

  ###
  # Differential features in terms of variability
  ##
  diff_var_feat = diff_var_feat_up = diff_var_feat_down <- 0
  if (diff_feat_prcg_gamma != 0) {
    if (N_feat * diff_feat_prcg_gamma < 5) {
      stop("Too few differential features!\n",
           "Increase 'diff_feat_prcg_gamma' parameter.")
    }
    # Only define differential features that are are not in the edge cases
    selected_feats <- which(gamma_B > 0.05 & gamma_B < 0.95)
    if (length(selected_feats) == 0) {
      stop("Simulated gamma parameters are all near 0 or 1.\n",
           "Should simulate data across the whole range,\n",
           "since differential analysis in edge cases is problematic.")
    }
    # Randomly sample features that will change across groups
    diff_feat_idx <- sample(x = selected_feats,
                            size = round(N_feat * diff_feat_prcg_gamma))
    # Create overdispersion params that have increase in Odds Ratio
    pivot <- round(length(diff_feat_idx)/1.5)
    diff_up_idx <- diff_feat_idx[seq_len(pivot)]
    diff_down_idx <- diff_feat_idx[(pivot + 1):(length(diff_feat_idx))]
    gamma_B[diff_up_idx] <- stats::plogis(
      stats::qlogis(gamma_B[diff_up_idx]) + log(OR_change_gamma) )
    gamma_B[diff_down_idx] <- stats::plogis(
      stats::qlogis(gamma_B[diff_down_idx]) - log(OR_change_gamma) )

    # Create data.frame with differential features
    diff_var_feat <- data.frame(
      feature_name = rownames(sim_dt_A$X[diff_feat_idx, ]),
      feature_idx = diff_feat_idx)
    diff_var_feat_up <- data.frame(
      feature_name = rownames(sim_dt_A$X[diff_up_idx, ]),
      feature_idx = diff_up_idx)
    diff_var_feat_down <- data.frame(
      feature_name = rownames(sim_dt_A$X[diff_down_idx, ]),
      feature_idx = diff_down_idx)
  }

  ##
  # Differential features in terms of mean methylation
  ##
  diff_mean_feat = diff_mean_feat_up = diff_mean_feat_down <- 0
  if (diff_feat_prcg_mu != 0) {
    if (N_feat * diff_feat_prcg_mu < 5) {
      stop("Too few differential features!\n",
           "Increase 'diff_feat_prcg_gamma' parameter.")
    }

    # Only define differential features that are are not in the edge cases
    selected_feats <- which(mu_B > 0.1 & mu_B < 0.9)
    if (length(selected_feats) == 0) {
      stop("Simulated gamma parameters are all near 0 or 1.\n",
           "Should simulate data across the whole range,\n",
           "since differential analysis in edge cases is problematic.")
    }
    # Randomly sample features that will change across groups
    diff_feat_idx <- sample(x = selected_feats,
                            size = round(N_feat * diff_feat_prcg_mu))
    # Create mean params that have increase in Odds Ratio
    pivot <- round(length(diff_feat_idx)/1.5)
    diff_up_idx <- diff_feat_idx[seq_len(pivot)]
    diff_down_idx <- diff_feat_idx[(pivot + 1):(length(diff_feat_idx))]
    mu_B[diff_up_idx] <- stats::plogis(
      stats::qlogis(mu_B[diff_up_idx]) + log(OR_change_mu) )
    mu_B[diff_down_idx] <- stats::plogis(
      stats::qlogis(mu_B[diff_down_idx]) - log(OR_change_mu) )

    # Create data.frame with differential features
    diff_mean_feat <- data.frame(
      feature_name = rownames(sim_dt_A$X[diff_feat_idx, ]),
      feature_idx = diff_feat_idx)
    diff_mean_feat_up <- data.frame(
      feature_name = rownames(sim_dt_A$X[diff_up_idx, ]),
      feature_idx = diff_up_idx)
    diff_mean_feat_down <- data.frame(
      feature_name = rownames(sim_dt_A$X[diff_down_idx, ]),
      feature_idx = diff_down_idx)
  }

  # Generate number of methylated CpGs from Beta Binomial
  met_cpgs_list_B <- lapply(X = seq_len(N_feat), function(n)
    VGAM::rbetabinom(length(cpgs_list_B[[n]]), cpgs_list_B[[n]],
                     prob = mu_B[n], rho = gamma_B[n]))

  # Create observed data object, which will be used as input to the model
  Y_B <- .generate_ys(N_feat = N_feat, N_cells = N_cells,
                      cpgs_list = cpgs_list_B, met_cpgs_list = met_cpgs_list_B,
                      feature_names = unique(sim_dt_A$Y$Feature))
  # Store parameters
  theta <- data.frame(mu = mu_B, gamma = gamma_B,
                      row.names = unique(Y_B$Feature))
  theta_priors <- list(w_mu = w_mu, w_gamma = w_gamma, s_mu = s_mu,
                       s_gamma = s_gamma)
  # Store data from group B
  sim_dt_B <- list(Y = Y_B, X = sim_dt_A$X, theta_true = theta,
                   theta_priors_true = theta_priors, opts = sim_dt_A$opts)

  # Options for the differential analysis
  opts <- list(diff_feat_prcg_mu = diff_feat_prcg_mu,
               diff_feat_prcg_gamma = diff_feat_prcg_gamma,
               OR_change_mu = OR_change_mu, OR_change_gamma = OR_change_gamma)
  obj <- structure(list(scmet_dt_A = sim_dt_A, scmet_dt_B = sim_dt_B,
                        opts = opts, diff_var_features = diff_var_feat,
                        diff_var_features_up = diff_var_feat_up,
                        diff_var_features_down = diff_var_feat_down,
                        diff_mean_features = diff_mean_feat,
                        diff_mean_features_up = diff_mean_feat_up,
                        diff_mean_features_down = diff_mean_feat_down),
                   class = "scmet_simulate_diff")
  return(obj)
}

# Internal function to generate total number of CpGs for each feature and cell.
# The final object is a list of size N_feat (total number of features). Then for
# each feature we have a vector whose size denotes total number of cells that
# have CpG coverage (<N_cells), and each value denotes the total number of CpGs.
.generate_cpgs <- function(N_feat = 100, N_cells = 50, N_cpgs = 15,
                           cells_range = c(0.4,0.8), cpgs_range = c(0.4,0.8)) {
  # Total number of cells, each feature has a different # of cells
  cell_v <- stats::rbinom(N_feat, N_cells,
                          stats::runif(N_feat, min = cells_range[1],
                                       max = cells_range[2]))
  # If we have featyres with less than 4 cells
  idx <- which(cell_v < 4)
  # ... sample values from [4, 6].
  if (length(idx) > 0) {
    cell_v[idx] <- sample(4:6, length(idx), replace = TRUE)
  }
  # Total number of CpGs for each cell
  cpgs_list <- lapply(cell_v, function(x) {
    n <- stats::rbinom(x, N_cpgs, stats::runif(1, min = cpgs_range[1],
                                               max = cpgs_range[2]))
    # If we have features with less than 3 CpGs, ...
    idx <- which(n < 3)
    # ... sample values from [3, 5].
    if (length(idx) > 0) {
      n[idx] <- sample(3:5, length(idx), replace = TRUE)
    }
    return(n)
  })
  return(cpgs_list)
}

# Internal function to generate mean methylation levels \mu.
# There is choice to make \mu depend on some covariates X.
.generate_means <- function(N_feat, N_cpgs, X = NULL, w_mu, s_mu) {
  # If we are not given covariates X
  if (is.null(X)) {
    # Randomly sample probabilities of success
    p <- stats::rbeta(N_feat, 1, 2)
    p[p < 0.05] <- 0.05
    # Sample # of CpGs and log transform
    tmp <- stats::rbinom(N_feat, 10*N_cpgs, p)
    # In case we have zeros
    tmp[tmp < 1] <- 1
    X <- log(tmp)
    # Perform mean centering
    X <- X - mean(X)
    # Add bias term
    X <- cbind(rep(1, N_feat), X)
    colnames(X) <- c("intercept", "cpg_density")
  } else {
    if (NCOL(X) != length(w_mu)) {
      stop("Number of columns of X should match w_mu length!")
    }
  }
  # Generate means mu
  mu <- logitnorm::rlogitnorm(N_feat, mu = X %*% w_mu, sigma = s_mu)
  return(list(mu = mu, X = X))
}

# Internal function to generate dispersion parameter \gamma.
.generate_overdisp <- function(N_feat, L, mu, w_gamma, s_gamma, rbf_c) {
  # Create design matrix
  H <- create_design_matrix(L = L, X = mu, c = rbf_c)
  # Generate dispersion parameters gamma
  if (NCOL(H) != length(w_gamma)) {
    stop("Number of basis functions should match w_gamma length!")
  }
  gamma <- logitnorm::rlogitnorm(N_feat, H %*% w_gamma, sigma = s_gamma)
  return(gamma)
}

# Create object of observed data Y in the format used by scMET
.generate_ys <- function(N_feat, N_cells, cpgs_list, met_cpgs_list,
                         feature_names = NULL) {
  # FIX: 2020/01/25: Fixed issue of biasing last features to having less cells.
  # FIX: 2020/03/03: Fixed issue when simulating differential data and features
  #         are ordered on the first dataset, now adding an additional
  #         parameter with feature_names.
  if (is.null(feature_names)) {
    Y <- data.table::data.table(
      "Feature" = unlist(lapply(
        X = seq_len(N_feat),
        FUN = function(n) rep(paste0("Feature_", n), length(cpgs_list[[n]])))),
      "Cell" = unlist(lapply(
        X = seq_len(N_feat),
        FUN = function(n) {
          cell_id <- sample(N_cells, length(cpgs_list[[n]]))
          paste0("Cell_", cell_id)
          }) ) )
  }
  else {
    Y <- data.table::data.table(
      "Feature" = unlist(lapply(
        X = seq_len(N_feat),
        FUN = function(n) rep(feature_names[n], length(cpgs_list[[n]])))),
      "Cell" = unlist(lapply(
        X = seq_len(N_feat),
        FUN = function(n) {
          cell_id <- sample(N_cells, length(cpgs_list[[n]]))
          paste0("Cell_", cell_id)
          }) ) )
  }
  Y[, c("total_reads", "met_reads") :=
      list(unlist(cpgs_list), unlist(met_cpgs_list)) ]
  return(Y)
}

# Initialize w_gamma parameter, depending on other input parameters
.init_w_gamma <- function(w_gamma = NULL, L = 4) {
  # Randomly initialize w_gamma parameter
  if (is.null(w_gamma)) {
    if (L == 4) {
      w_gamma <- c(-1.2, -.3, 1.1, -.9)
    } else {
      w_gamma <- stats::rnorm(L, mean = 0, sd = 1)
    }
  } else {
    if (length(w_gamma) != L) {
      stop("Number of basis functions L should match `w_gamma` length!")
    }
  }
  return(w_gamma)
}
