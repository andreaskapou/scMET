#' @name scmet_differential
#' @aliases differential_test, differential_methylation,
#'   differential_variability
#'
#' @title Differential testing using scMET
#'
#' @description Function for performing differential methylation testing to
#'   identify differentially methylted (DM) and differentially variable (DV)
#'   features across two groups of pre-specified cell populations.
#'
#' @param obj_A The scMET posterior object for group A.
#' @param obj_B The scMET posterior object for group B.
#' @param psi_m Minimum log odds ratio tolerance threshold for detecting changes
#'   in overall methylation (positive real number). Default value: \code{psi_m =
#'   log(1.5)} (i.e. 50% increase).
#' @param psi_e Minimum log odds ratio tolerance threshold for detecting changes
#'   in residual over-dispersion (positive real number).
#' @param psi_g Minimum log odds ratio tolerance threshold for detecting changes
#'   in biological over-dispersion (positive real number).
#' @param evidence_thresh_m Optional parameter. Posterior evidence probability
#'   threshold parameter `alpha_{M}` for detecting changes in overall
#'   methylation (between 0.6 and 1). If \code{efdr_m = NULL}, then threshold
#'   will be set to \code{evidence_thresh_m}. If a value for \code{EFDR_M} is
#'   provided, the posterior probability threshold is chosen to achieve an EFDR
#'   equal to \code{efdr_m} and \code{evidence_thresh_m} defines a minimum
#'   probability threshold for this calibration (this avoids low values of
#'   \code{evidence_thresh_m} to be chosen by the EFDR calibration. Default
#'   value \code{evidence_thresh_m = 0.8}.
#' @param evidence_thresh_e Optional parameter. Posterior evidence probability
#'   threshold parameter `alpha_{G}` for detecting changes in cell-to-cell
#'   residual over-dispersion. Same usage as above.
#' @param evidence_thresh_g Optional parameter. Posterior evidence probability
#'   threshold parameter `alpha_{G}` for detecting changes in cell-to-cell
#'   biological over-dispersion. Same usage as above.
#' @param efdr_m Target for expected false discovery rate related to the
#'   comparison of means. If \code{efdr_m = NULL}, no calibration is performed,
#'   and `alpha_{M}` is set to \code{evidence_thresh_m}. Default value:
#'   \code{efdr_m = 0.05}.
#' @param efdr_e Target for expected false discovery rate related to the
#'   comparison of residual over-dispersions If \code{efdr_e = NULL}, no
#'   calibration is performed, and `alpha_{E}`` is set to
#'   \code{evidence_thresh_e}. Default value: \code{efdr_e = 0.05}.
#' @param efdr_g Target for expected false discovery rate related to the
#'   comparison of biological over-dispersions If \code{efdr_g = NULL}, no
#'   calibration is performed, and `alpha_{G}` is set to
#'   \code{evidence_thresh_g}. Default value: \code{efdr_g = 0.05}.
#' @param group_label_A Label assigned to group A.
#' @param group_label_B Label assigned to group B.
#' @param features_selected User defined list of selected features to perform
#'   differential analysis. Should be the same length as the total number of
#'   features, with TRUE for features included in the differential analysis, and
#'   FALSE for those excluded from further analysis.
#' @param filter_outlier_features Logical, whether to filter features that have
#'   either mean methylation levels `mu` or overdispersion `gamma` across both
#'   groups near the range edges, i.e. taking values near 0 or 1. This mostly is
#'   an issue due to taking the logit transformation which effectively makes
#'   small changes in actual space (0, 1) to look really large in transformed
#'   space (-Inf, Inf). In general we expect this will not remove many
#'   interesting features with biological information.
#' @param outlier_m Value of average mean methylation across both groups so a
#'   feature is considered as outlier. I.e. if set to 0.05, then will remove
#'   features with `mu` < 0.05 or `mu` > 1 - 0.05. Only used if
#'   `filter_outlier_features = TRUE`.
#' @param outlier_g Value of average overdispersion `gamma` across groups so a
#'   feature is considered as outlier. Same as `outlier_m` parameter above.
#'
#' @return An `scmet_differential` object which is a list containing the
#'   following elements: \itemize{ \item{ \code{diff_mu_summary}: A data.frame
#'   containing differential mean methylation output information per feature
#'   (rows), including posterior median parameters for each group and `mu_LOR`
#'   containing the log odds-ratio between the groups. The `mu_tail_prob` column
#'   contains the posterior tail probability of a feature being called as DM.
#'   The `mu_diff_test` column informs the outcomes of the test.} \item{
#'   \code{diff_epsilon_summary}: Same as above, but for differential
#'   variability based on residual overdispersion. } \item{
#'   \code{diff_gamma_summary}: The same as above but for DV analysis based on
#'   overdispersion.} \item{ \code{diff_mu_thresh}: Information about optimal
#'   posterior evidence threshold search for mean methylation mu. }
#'   \item{\code{diff_epsilon_thresh}: Same as above but for residual
#'   overdispersion epsilon..} \item{\code{diff_gamma_thresh}: Same as above but
#'   for overdispersion gamma.} \item{\code{opts}: The parameters used for
#'   testing. For reproducibility purposes.} }
#'
#' @seealso \code{\link{scmet}}, \code{\link{scmet_hvf_lvf}}
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#' @examples
#' # Fit scMET for each group
#' fit_A <- scmet(Y = scmet_diff_dt$scmet_dt_A$Y,
#' X = scmet_diff_dt$scmet_dt_A$X, L = 4, iter = 50, seed = 12)
#' fit_B <- scmet(Y = scmet_diff_dt$scmet_dt_B$Y,
#' X = scmet_diff_dt$scmet_dt_B$X, L = 4, iter = 50, seed = 12)
#'
#' # Run differential test
#' diff_obj <- scmet_differential(obj_A = fit_A, obj_B = fit_B)
#'
#' @export
#'
scmet_differential <- function(obj_A, obj_B, psi_m = log(1.5), psi_e = log(1.5),
                               psi_g = log(1.5), evidence_thresh_m = 0.8,
                               evidence_thresh_e = 0.8, evidence_thresh_g = 0.8,
                               efdr_m = 0.05, efdr_e = 0.05, efdr_g = 0.05,
                               group_label_A = "GroupA", group_label_B = "GroupB",
                               features_selected = NULL,
                               filter_outlier_features = FALSE,
                               outlier_m = 0.05, outlier_g = 0.05) {
  # So RMD check does not complain
  Feature <- NULL

  ##
  # Initial value checks
  ##
  if (!inherits(obj_A, c("scmet_mcmc", "scmet_vb"))) {
    stop("Posterior object is not generated from scMET.")
  }
  if (!inherits(obj_B, c("scmet_mcmc", "scmet_vb"))) {
    stop("Posterior object is not generated from scMET.")
  }

  # Mu
  if (!is.null(evidence_thresh_m)) {
    assertthat::assert_that(evidence_thresh_m >= 0.6 & evidence_thresh_m <= 1,
              msg = "Posterior evidence threshold must be between (0.6, 1).\n")
  } else {
    message("Posterior evidence threshold hasn't been supplied.\n",
            "Setting initial value to 0.8.\n")
    evidence_thresh_m <- 0.8
  }
  assertthat::assert_that(psi_m > 0)
  assertthat::assert_that(efdr_m >= 0 & efdr_m <= 1)

  # Epsilon
  if (!is.null(evidence_thresh_e)) {
    assertthat::assert_that(evidence_thresh_e >= 0.6 & evidence_thresh_e <= 1,
              msg = "Posterior evidence threshold must be between (0.6, 1).\n")
  } else {
    message("Posterior evidence threshold hasn't been supplied.\n",
            "Setting initial value to 0.8.\n")
    evidence_thresh_e <- 0.8
  }
  assertthat::assert_that(psi_e > 0)
  assertthat::assert_that(efdr_e >= 0 & efdr_e <= 1)

  # Gamma
  if (!is.null(evidence_thresh_g)) {
    assertthat::assert_that(evidence_thresh_g >= 0.6 & evidence_thresh_g <= 1,
              msg = "Posterior evidence threshold must be between (0.6, 1).\n")
  } else {
    message("Posterior evidence threshold hasn't been supplied.\n",
            "Setting initial value to 0.8.\n")
    evidence_thresh_g <- 0.8
  }
  assertthat::assert_that(psi_g > 0)
  assertthat::assert_that(efdr_g >= 0 & efdr_g <= 1)


  # Check if both groups have the same number of features and are ordered
  if (length(obj_A$feature_names) != length(obj_B$feature_names)) {
    message("Number of features does not match between groups.\n",
            "Keeping only matching features!\n")
    joint_feat <- intersect(obj_A$feature_names, obj_B$feature_names)
    idx_A <- which(obj_A$feature_names %in% joint_feat)
    idx_B <- which(obj_B$feature_names %in% joint_feat)

    obj_A$feature_names <- obj_A$feature_names[idx_A]
    obj_B$feature_names <- obj_B$feature_names[idx_B]

    obj_A$Y <- obj_A$Y[Feature %in% joint_feat]
    obj_B$Y <- obj_B$Y[Feature %in% joint_feat]

    obj_A$posterior$mu <- obj_A$posterior$mu[, idx_A]
    obj_A$posterior$gamma <- obj_A$posterior$gamma[, idx_A]
    obj_A$posterior$epsilon <- obj_A$posterior$epsilon[, idx_A]

    obj_B$posterior$mu <- obj_B$posterior$mu[, idx_B]
    obj_B$posterior$gamma <- obj_B$posterior$gamma[, idx_B]
    obj_B$posterior$epsilon <- obj_B$posterior$epsilon[, idx_B]

    # Removing the features selected mode option
    message("Ignoring the `features_selected` option as well.\n")
    features_selected <- NULL
  }

  if (!is.null(features_selected)) {
    assertthat::assert_that(all(features_selected %in% c(TRUE, FALSE)))
  }
  # If all features are to be included
  if (is.null(features_selected)) {
    features_selected <- rep(TRUE, times = length(obj_A$feature_names))
  }

  ##----------------
  ## Extract summary information
  ##----------------

  # Number of cells for group A
  c_A <- length(unique(obj_A$Y$Cell))
  # Number of cells for group B
  c_B <- length(unique(obj_B$Y$Cell))
  # Total number of cells
  c <- c_A + c_B

  factor_levels <- c(paste0(group_label_A, "+"), paste0(group_label_B, "+"),
                     "NoDiff", "ExcludedByUser", "ExcludedFromTesting")

  # Extract posterior median for mean methylation, check for outliers
  mu_A <- matrixStats::colMedians(obj_A$posterior$mu)
  mu_B <- matrixStats::colMedians(obj_B$posterior$mu)
  mu_A <- .fix_outliers(x = mu_A, xmin = 1e-02, xmax = 1 - 1e-2)
  mu_B <- .fix_outliers(x = mu_B, xmin = 1e-02, xmax = 1 - 1e-2)
  mu_overall <- (mu_A * c_A + mu_B * c_B) / c

  # Extract posterior median for biological over-dispersion, check for outliers
  gamma_A <- matrixStats::colMedians(obj_A$posterior$gamma)
  gamma_B <- matrixStats::colMedians(obj_B$posterior$gamma)
  gamma_A <- .fix_outliers(x = gamma_A, xmin = 1e-02, xmax = 1 - 1e-2)
  gamma_B <- .fix_outliers(x = gamma_B, xmin = 1e-02, xmax = 1 - 1e-2)
  gamma_overall <- (gamma_A + gamma_B) / 2

  # Extract posterior median for residual over-dispersion, check for outliers
  epsilon_A <- matrixStats::colMedians(obj_A$posterior$epsilon)
  epsilon_B <- matrixStats::colMedians(obj_B$posterior$epsilon)
  epsilon_A <- .fix_outliers(x = epsilon_A, xmin = -7, xmax = 7)
  epsilon_B <- .fix_outliers(x = epsilon_B, xmin = -7, xmax = 7)
  epsilon_overall <- (epsilon_A + epsilon_B) / 2

  ##----------------
  ## Changes in mean methylation
  ##----------------

  # Chain with log(Odds Ratio)
  chain_lor_m <- stats::qlogis(obj_A$posterior$mu) -
    stats::qlogis(obj_B$posterior$mu)
  # Compute posterior tail probabilities for difference in mean methylation
  prob_m <- .tail_prob(chain = chain_lor_m, tolerance_thresh = psi_m)
  # Compute posterior median of LORs
  median_lor_m <- matrixStats::colMedians(chain_lor_m)
  median_lor_m <- .fix_outliers(x = median_lor_m, xmin = -7, xmax = 7)

  # Remove features from DM/DV testing that are considered as outliers,
  # due to the logit link issue at the edge cases near 0 or 1.
  features_selected_mu <- features_selected
  if (filter_outlier_features) {
    idx <- which(mu_overall < outlier_m | mu_overall > (1 - outlier_m))
    features_selected_mu[idx] <- FALSE
  }

  # Search optimal threshold to identify changes in mean between groups
  mu_thresh <- .thresh_search(evidence_thresh = evidence_thresh_m,
                              prob = prob_m[features_selected_mu],
                              efdr = efdr_m, task = "Differential mean",
                              suffix = "m")
  opt_evidence_thresh_m <- mu_thresh$optimal_evidence_thresh

  # Obtain differential mean test results
  res_diff_mu <- .diff_test_results(prob = prob_m,
        evidence_thresh = opt_evidence_thresh_m[1], estimate = median_lor_m,
        group_label_A = group_label_A, group_label_B = group_label_B,
        features_selected = features_selected_mu, excluded = NULL)
  # Output table
  tbl_mu <- cbind.data.frame(feature_name = obj_A$feature_names,
                             mu_overall = as.numeric(mu_overall),
                             mu_A = mu_A,
                             mu_B = mu_B,
                             mu_LOR = as.numeric(median_lor_m),
                             mu_OR = as.numeric(exp(median_lor_m)),
                             mu_tail_prob = prob_m,
                             mu_diff_test = factor(res_diff_mu,
                                                   levels = factor_levels),
                             stringsAsFactors = FALSE)
  tbl_mu <- tbl_mu[order(tbl_mu$mu_tail_prob,
                         decreasing = TRUE, na.last = TRUE), ]


  ##----------------
  ## Changes in over-dispersion
  ##----------------
  # Chain with log(Odds Ratio)
  chain_lor_g <- stats::qlogis(obj_A$posterior$gamma) -
    stats::qlogis(obj_B$posterior$gamma)
  # Compute posterior tail probabilities for difference in over-dispersion
  prob_g <- .tail_prob(chain = chain_lor_g, tolerance_thresh = psi_g)
  # Compute posterior median of LORs
  median_lor_g <- matrixStats::colMedians(chain_lor_g)
  median_lor_g <- .fix_outliers(x = median_lor_g, xmin = -7, xmax = 7)

  # Remove features from DM/DV testing that are considered as outliers
  features_selected_gamma <- features_selected
  if (filter_outlier_features) {
    idx <- which(gamma_overall < outlier_g | gamma_overall > (1 - outlier_g))
    features_selected_gamma[idx] <- FALSE
  }

  # Features with no change in mean methylation
  not_dm <- res_diff_mu == "NoDiff"
  # Features to calibrate EFDR
  feat_mu_diff <- not_dm & features_selected_gamma

  # Search optimal threshold to identify changes in mean between groups
  gamma_thresh <- .thresh_search(evidence_thresh = evidence_thresh_g,
                                 prob = prob_g[feat_mu_diff], efdr = efdr_g,
                                 task = "Differential overdispersion",
                                 suffix = "g")
  opt_evidence_thresh_g <- gamma_thresh$optimal_evidence_thresh

  # Obtain differential over-dispersion test results
  res_diff_gamma <- .diff_test_results(prob = prob_g,
          evidence_thresh = opt_evidence_thresh_g[1], estimate = median_lor_g,
          group_label_A = group_label_A, group_label_B = group_label_B,
          features_selected = features_selected_gamma, excluded = !not_dm)
  # Output table
  tbl_gamma <- cbind.data.frame(feature_name = obj_A$feature_names,
                                gamma_overall = as.numeric(gamma_overall),
                                gamma_A = gamma_A,
                                gamma_B = gamma_B,
                                gamma_LOR = as.numeric(median_lor_g),
                                gamma_OR = as.numeric(exp(median_lor_g)),
                                gamma_tail_prob = prob_g,
                                gamma_diff_test = factor(res_diff_gamma,
                                                         levels = factor_levels),
                                mu_overall = as.numeric(mu_overall),
                                mu_A = mu_A,
                                mu_B = mu_B,
                                stringsAsFactors = FALSE)
  tbl_gamma <- tbl_gamma[order(tbl_gamma$gamma_tail_prob,
                               decreasing = TRUE, na.last = TRUE), ]

  ##----------------
  ## Changes in residual over-dispersion
  ##----------------

  # Chain with log(Odds Ratio)
  chain_lor_e <- obj_A$posterior$epsilon - obj_B$posterior$epsilon
  # Compute tail probabilities for difference in residual overdispersion
  prob_e <- .tail_prob(chain = chain_lor_e, tolerance_thresh = psi_e)
  # Compute posterior median of LORs
  median_lor_e <- matrixStats::colMedians(chain_lor_e)
  median_lor_e <- .fix_outliers(x = median_lor_e, xmin = -7, xmax = 7)

  # Search optimal evidence threshold to identify changes in epsilon
  epsilon_thresh <- .thresh_search(evidence_thresh = evidence_thresh_e,
                                   prob = prob_e[features_selected_gamma],
                                   efdr = efdr_e,
                                   task = "Differential residual overdispersion",
                                   suffix = "e")
  opt_evidence_thresh_e <- epsilon_thresh$optimal_evidence_thresh

  # Obtain differential over-dispersion test results
  res_diff_epsilon <- .diff_test_results(prob = prob_e,
        evidence_thresh = opt_evidence_thresh_e[1], estimate = median_lor_e,
        group_label_A = group_label_A, group_label_B = group_label_B,
        features_selected = features_selected_gamma, excluded = NULL)
  # Output table
  tbl_epsilon <- cbind.data.frame(feature_name = obj_A$feature_names,
                                  epsilon_overall = as.numeric(epsilon_overall),
                                  epsilon_A = epsilon_A,
                                  epsilon_B = epsilon_B,
                                  epsilon_change = as.numeric(median_lor_e),
                                  epsilon_tail_prob = prob_e,
                                  epsilon_diff_test = factor(res_diff_epsilon,
                                                             levels = factor_levels),
                                  mu_overall = as.numeric(mu_overall),
                                  mu_A = mu_A,
                                  mu_B = mu_B,
                                  stringsAsFactors = FALSE)
  tbl_epsilon <- tbl_epsilon[order(tbl_epsilon$epsilon_tail_prob,
                                   decreasing = TRUE, na.last = TRUE), ]

  if (!is.null(features_selected)) {
    message("--------------------------------------- \n",
            "The user excluded ", sum(!features_selected), " features.\n",
            "These features are marked as 'ExcludedByUser' \n",
            "and are excluded from EFDR calibration.\n")
  }

  if (filter_outlier_features) {
    message("--------------------------------------- \n",
            "Additionally, the user excluded ", sum(!features_selected_mu),
            " due to mean methylation outlier features.\n",
            "Additionally, the user excluded ", sum(!features_selected_gamma),
            " due to overdispersion outlier features.\n",
            "These features are marked as 'ExcludedByUser' \n",
            "and are excluded from EFDR calibration.\n")
  }

  # Summary of total hits
  n_mu_high_A <- sum(res_diff_mu == paste0(group_label_A, "+"))
  n_mu_high_B <- sum(res_diff_mu == paste0(group_label_B, "+"))
  n_gamma_high_A <- sum(res_diff_gamma == paste0(group_label_A, "+"))
  n_gamma_high_B <- sum(res_diff_gamma == paste0(group_label_B, "+"))
  n_epsilon_high_A <- sum(res_diff_epsilon == paste0(group_label_A, "+"))
  n_epsilon_high_B <- sum(res_diff_epsilon == paste0(group_label_B, "+"))


  message("-------------------------------------------------------------\n",
          n_mu_high_A + n_mu_high_B,
          " features with a change in mean methylation:\n",
          "- Higher methylation in ", group_label_A,
          " samples: ", n_mu_high_A, "\n",
          "- Higher methylation in ", group_label_B,
          " samples: ", n_mu_high_B, "\n",
          "- Odds ratio tolerance = ", exp(psi_m), "\n",
          "- Probability evidence threshold=", opt_evidence_thresh_m[1], "\n",
          "- EFDR = ", round(100 * opt_evidence_thresh_m[2], 2), "% \n",
          "- EFNR = ", round(100 * opt_evidence_thresh_m[3], 2), "%. \n",
          "-------------------------------------------------------------\n\n",
          "-------------------------------------------------------------\n",
          n_gamma_high_A + n_gamma_high_B,
          " features with a change in over-dispersion:\n",
          "- Higher dispersion in ", group_label_A,
          " samples: ", n_gamma_high_A,"\n",
          "- Higher dispersion in ", group_label_B,
          " samples: ", n_gamma_high_B,"\n",
          "- Odds ratio tolerance = ", exp(psi_g), "\n",
          "- Probability evidence threshold=", opt_evidence_thresh_g[1], "\n",
          "- EFDR = ", round(100 * opt_evidence_thresh_g[2], 2), "% \n",
          "- EFNR = ", round(100 * opt_evidence_thresh_g[3], 2), "% \n",
          "NOTE: differential dispersion assessment only applied to the \n",
          sum(not_dm), " features for which the mean did not change \n",
          "and that were included for testing. \n",
          "--------------------------------------------------------------\n",
          "-------------------------------------------------------------\n",
          n_epsilon_high_A + n_epsilon_high_B,
          " features with a change in residual over dispersion:\n",
          "- Higher residual dispersion in ", group_label_A,
          " samples: ", n_epsilon_high_A,"\n",
          "- Higher residual dispersion in ", group_label_B,
          " samples: ", n_epsilon_high_B,"\n",
          "- Odds ratio tolerance = ", exp(psi_e), "\n",
          "- Probability evidence threshold=", opt_evidence_thresh_e[1], "\n",
          "- EFDR = ", round(100 * opt_evidence_thresh_e[2], 2), "% \n",
          "- EFNR = ", round(100 * opt_evidence_thresh_e[3], 2), "%. \n",
          "--------------------------------------------------------------\n")

  obj <- structure(list(diff_mu_summary = tbl_mu,
                        diff_epsilon_summary = tbl_epsilon,
                        diff_gamma_summary = tbl_gamma,
    diff_mu_thresh = list(evidence_thresh = opt_evidence_thresh_m[1],
                          efdr = opt_evidence_thresh_m[2],
                          efnr = opt_evidence_thresh_m[3],
                          efdr_grid = mu_thresh$efdr_grid,
                          efnr_grid = mu_thresh$efnr_grid,
                          evidence_thresh_grid = mu_thresh$evidence_thresh_grid,
                          target_efdr = efdr_m),
    diff_epsilon_thresh = list(evidence_thresh = opt_evidence_thresh_e[1],
                               efdr = opt_evidence_thresh_e[2],
                               efnr = opt_evidence_thresh_e[3],
                               efdr_grid = epsilon_thresh$efdr_grid,
                               efnr_grid = epsilon_thresh$efnr_grid,
                               evidence_thresh_grid = epsilon_thresh$evidence_thresh_grid,
                               target_efdr = efdr_e),
    diff_gamma_thresh = list(evidence_thresh = opt_evidence_thresh_g[1],
                            efdr = opt_evidence_thresh_g[2],
                            efnr = opt_evidence_thresh_g[3],
                            efdr_grid = gamma_thresh$efdr_grid,
                            efnr_grid = gamma_thresh$efnr_grid,
                            evidence_thresh_grid = gamma_thresh$evidence_thresh_grid,
                            target_efdr = efdr_g),
    opts = list(psi_m = psi_m,
                psi_g = psi_g,
                psi_e = psi_e,
                evidence_thresh_m = evidence_thresh_m,
                evidence_thresh_g = evidence_thresh_g,
                evidence_thresh_e = evidence_thresh_e,
                efdr_m = efdr_m,
                efdr_g = efdr_g,
                efdr_e = efdr_e,
                group_label_A = group_label_A,
                group_label_B = group_label_B)),
    class = "scmet_differential")
  return(obj)

}
