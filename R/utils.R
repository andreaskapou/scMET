#' @title Create design matrix
#'
#' @description Generic function for crating a radial basis function (RBF)
#'   design matrix for input vector X.
#'
#' @param L Total number of basis functions, including the bias term.
#' @param X Vector of covariates
#' @param c Scaling parameter for variance of RBFs
#'
#' @return A design matrix object H.
#'
#' @seealso \code{\link{scmet}}, \code{\link{scmet_differential}},
#'   \code{\link{scmet_hvf_lvf}}
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#' @export
create_design_matrix <- function(L, X, c = 1.2) {
  H <- .rbf_design_matrix(L = L, X = X, c = c)
  return(H)
}


# RBF evaluation
.rbf_basis <- function(X, mus, h = 1){
  return(exp( -0.5 * sum(((X - mus) / h)^2) ))
}


# @title Create RBF design matrix
#
# @param L Total number of basis functions, including the bias term
# @param X Vector of covariates
# @param c Scaling parameter for variance of RBFs
.rbf_design_matrix <- function(L, X, c = 1){
  N <- length(X)  # Length of the dataset
  if (L > 1) {
    # Compute mean locations
    ms <- rep(0, L - 1)
    for (l in 1:(L - 1)) {
      ms[l] <- l * ((max(X) - min(X)) / L ) + min(X)
    }
    # Compute scaling parameter
    h <- (ms[2] - ms[1]) * c
    H <- matrix(1, nrow = N, ncol = L)
    for (l in 1:(L - 1)) {
      H[, l + 1] <- apply(as.matrix(X), 1, .rbf_basis, mus = ms[l], h = h)
    }
  } else {
    H <- matrix(1, nrow = N, ncol = L)
  }
  return(H)
}


# @title Create polynomial design matrix
#
# @param L The degree of the polynomial basis that will be applied to input X.
# @param X Vector of covariates
.poly_design_matrix <- function(L, X) {
  H <- matrix(1, nrow = length(X), ncol = L)
  if (L > 1) {
    for (l in 2:L) {
      H[, l] <- X ^ (l - 1)
    }
  }
  return(H)
}


# Infer penalized linear regression model
.lm_mle_penalized <- function(y, H, lambda = 0.5){
  if (lambda == 0) {
    qx <- qr(H)             # Compute QR decomposition of H
    W <- c(solve.qr(qx, y)) # Compute (H'H)^(-1)H'y
  }else{
    I <- diag(1, NCOL(H))   # Identity matrix
    # TODO: Should we do this or not for the bias term??
    # I[1,1]  <- 1e-10  # Do not change the intercept coefficient
    qx <- qr(lambda * I + t(H) %*% H) # Compute QR decomposition
    W <- c(solve.qr(qx, t(H) %*% y))  # Comp (lambda*I+H'H)^(-1)H'y
  }
  return(c(W))
}


# Evaluate EFDR for given evidence threshold and posterior tail probabilities
# Adapted from BASiCS package.
.eval_efdr <- function(evidence_thresh, prob) {
  return(sum((1 - prob) * I(prob > evidence_thresh)) /
           sum(I(prob > evidence_thresh)))
}


# Evaluate EFNR for given evidence threshold and posterior tail probabilities
# Adapted from BASiCS package.
.eval_efnr <- function(evidence_thresh, prob) {
  return(sum(prob * I(evidence_thresh >= prob)) /
           sum(I(evidence_thresh >= prob)))
}


# Compute posterior tail probabilities for the differential analysis task.
# Adapted from BASiCS. See also Bochina and Richardson (2007)
.tail_prob <- function(chain, tolerance_thresh) {
  if (tolerance_thresh > 0) {
    prob <- matrixStats::colMeans2(ifelse(abs(chain) > tolerance_thresh, 1, 0))
  } else {
    tmp <- matrixStats::colMeans2(ifelse(abs(chain) > 0, 1, 0))
    prob <- 2 * pmax(tmp, 1 - tmp) - 1
  }
  return(prob)
}


# Search function for optimal posterior evidence threshold \alpha
# Adapted from BASiCS
.thresh_search <- function(evidence_thresh, prob, efdr, task, suffix = "") {
  # Summary of cases
  # 1. If EFDR is provided - run calibration
  #   1.1. If the calibration doesn't completely fail - search \alpha
  #     1.1.1. If optimal \alpha is not too low - set \alpha to optimal
  #     1.1.2. If optimal \alpha is too low - fix to input probs
  #   1.2 If calibration completely fails - default \alpha=0.9 (conservative)
  # 2. If EFDR is not provided - fix to input probs


  # 1. If EFDR is provided - run calibration
  if (!is.null(efdr)) {
    # If threshold is not set a priori (search)
    evidence_thresh_grid <- seq(0.6, 0.9995 , by = 0.00025)
    # Evaluate EFDR and EFNR on this grid
    efdr_grid <- vapply(evidence_thresh_grid, FUN = .eval_efdr,
                        FUN.VALUE = 1, prob = prob)
    efnr_grid <- vapply(evidence_thresh_grid, FUN = .eval_efnr,
                        FUN.VALUE = 1, prob = prob)

    # Compute absolute difference between supplied EFDR and grid search
    abs_diff <- abs(efdr_grid - efdr)
    # If we can estimate EFDR
    if (sum(!is.na(abs_diff)) > 0) {
      # Search EFDR closest to the desired value
      efdr_optimal <- efdr_grid[abs_diff == min(abs_diff, na.rm = TRUE) &
                                  !is.na(abs_diff)]
      # If multiple threholds lead to same EFDR, choose the lowest EFNR
      efnr_optimal <- efnr_grid[efdr_grid == mean(efdr_optimal) &
                                  !is.na(efdr_grid)]
      if (sum(!is.na(efnr_optimal)) > 0) {
        optimal <- which(efdr_grid == mean(efdr_optimal) &
                           efnr_grid == mean(efnr_optimal))
      } else {
        optimal <- which(efdr_grid == mean(efdr_optimal))
      }
      # Quick fix for EFDR/EFNR ties; possibly not an issue in real datasets
      optimal <- stats::median(round(stats::median(optimal)))

      # If calibrated threshold is above the minimum required probability
      if (evidence_thresh_grid[optimal] > evidence_thresh) {
        # 1.1.1. If optimal prob is not too low - set prob to optimal
        optimal_evidence_thresh <- c(evidence_thresh_grid[optimal],
                                     efdr_grid[optimal], efnr_grid[optimal])
        if (abs(optimal_evidence_thresh[2] - efdr) > 0.025) {
          # Message when different to desired EFDR is large
          message("For ", task, " task:\n",
                  "Not possible to find evidence probability threshold (>0.6)",
                  "\n that achieves desired EFDR level (tolerance +- 0.025)\n",
                  "Output based on the closest possible value. \n")
        }
      } else {
        # 1.1.2. If optimal prob is too low - fix to input probs
        efdr_grid <- .eval_efdr(evidence_thresh = evidence_thresh, prob = prob)
        efnr_grid <- .eval_efnr(evidence_thresh = evidence_thresh, prob = prob)
        optimal_evidence_thresh <- c(evidence_thresh, efdr_grid[1],
                                     efnr_grid[1])

        # Only required for differential test function
        if (suffix != "") { suffix <- paste0("_", suffix) }
        message("For ", task, " task:\n",
                "Evidence probability threshold chosen via EFDR valibration",
                " is too low. \n", "Probability threshold set automatically",
                " equal to 'evidence_thresh", suffix, "'.\n")
      }
    }
    else {
      # 1.2 If calibration completely fails - default prob = 0.9 (conservative)
      message("EFDR calibration failed for ", task, " task. \n",
              "Evidence probability threshold set equal to 0.9. \n")
      optimal_evidence_thresh <- c(0.90, NA, NA)
    }
  } else {
    # 2. If EFDR is not provided - fix to given probs
    efdr_grid <- .eval_efdr(evidence_thresh = evidence_thresh, prob = prob)
    efnr_grid <- .eval_efnr(evidence_thresh = evidence_thresh, prob = prob)
    optimal_evidence_thresh <- c(evidence_thresh, efdr_grid[1], efnr_grid[1])
    evidence_thresh_grid <- NULL
  }

  return(list("optimal_evidence_thresh" = optimal_evidence_thresh,
              "evidence_thresh_grid" = evidence_thresh_grid,
              "efdr_grid" = efdr_grid, "efnr_grid" = efnr_grid))
}


# Internal function to collect differential test results
.diff_test_results <- function(prob, evidence_thresh, estimate, group_label_A,
                               group_label_B, features_selected,
                               excluded = NULL) {

  # Which features are + in each group
  high_A <- which(prob > evidence_thresh & estimate > 0)
  high_B <- which(prob > evidence_thresh & estimate < 0)
  res_diff <- rep("NoDiff", length(estimate))
  res_diff[high_A] <- paste0(group_label_A, "+")
  res_diff[high_B] <- paste0(group_label_B, "+")
  if (!is.null(excluded)) { res_diff[excluded] <- "ExcludedFromTesting" }
  res_diff[!features_selected] <- "ExcludedByUser"
  return(res_diff)
}


# Create confusion matrix. Used only for assessing simulated data where ground
# truth information is present.
# TODO: Rewrite this function!
.confusion_matrix <- function(dt, diff_analysis, N_cells, N_feat, b = 1) {
  fp = fn = tp = tn <- matrix(0, ncol = 3, nrow = length(N_cells))
  colnames(fp) = colnames(fn) = colnames(tn) = colnames(tp) <- c("mu", "gamma", "epsilon")
  rownames(fp) = rownames(fn) = rownames(tn) = rownames(tp) <- paste0("Cells", N_cells)
  for (i in 1:length(N_cells)) {
    if (is.numeric(dt[[i]]$sim_dt$diff_var_features)) {
      diff_var_feat <- 0
    } else {
      diff_var_feat <- dt[[i]]$sim_dt$diff_var_features$feature_idx
    }
    if (is.numeric(dt[[i]]$sim_dt$diff_mean_features) ) {
      diff_mean_feat <- 0
    } else {
      diff_mean_feat <- dt[[i]]$sim_dt$diff_mean_features$feature_idx
    }

    hits <- list(which(diff_analysis[[i]]$mean_summary$mean_diff_test %in%
                         c( paste0(diff_analysis[[i]]$opts$group_label_B, "+"),
                            paste0(diff_analysis[[i]]$opts$group_label_A, "+") )),
              which(diff_analysis[[i]]$disp_summary$disp_diff_test %in%
                      c( paste0(diff_analysis[[i]]$opts$group_label_B, "+"),
                         paste0(diff_analysis[[i]]$opts$group_label_A, "+") ) ),
              which(diff_analysis[[i]]$res_disp_summary$res_disp_diff_test %in%
                      c( paste0(diff_analysis[[i]]$opts$group_label_B, "+"),
                         paste0(diff_analysis[[i]]$opts$group_label_A, "+") ) ))
    tp[i, ] <- c( sum(diff_mean_feat %in% hits[[1]]),
                  sum(diff_var_feat %in% hits[[2]]),
                  sum(diff_var_feat %in% hits[[3]]) )
    fp[i, ] <- c( sum(setdiff(seq(1, N_feat), diff_mean_feat) %in% hits[[1]]),
                  sum(setdiff(seq(1, N_feat), diff_var_feat) %in% hits[[2]]),
                  sum(setdiff(seq(1, N_feat), diff_var_feat) %in% hits[[3]]) )

    tn[i, ] <- c( sum(setdiff(seq(1, N_feat), diff_mean_feat) %in%
                        setdiff(seq(1, N_feat), hits[[1]])),
                  sum(setdiff(seq(1, N_feat), diff_var_feat) %in%
                        setdiff(seq(1, N_feat), hits[[2]])),
                  sum(setdiff(seq(1, N_feat), diff_var_feat) %in%
                        setdiff(seq(1, N_feat), hits[[3]])) )

    fn[i, ] <- c( sum(diff_mean_feat %in% setdiff(seq(1, N_feat), hits[[1]])),
                  sum(diff_var_feat %in% setdiff(seq(1, N_feat), hits[[2]])),
                  sum(diff_var_feat %in% setdiff(seq(1, N_feat), hits[[3]])) )
  }

  # FPR: # of incorrect predicted positives divided by
  #       total number of negatives (1 - specificity)
  fpr <- fp / (fp + tn)
  # FNR: # of incorrect predicted negatives divided by
  #       total number of positives (1 - recall)
  fnr <- fn / (fn + tp)
  # FDR: # of incorrect predicted positives divided by
  #       total number of discoveries
  fdr <- fp / (fp + tp)
  precision <- tp / (tp + fp)
  recall <- tp / (tp + fn) # True positive rate or sensitivity
  specificity <- tn / (tn + fp)
  f1_measure <- 2 * ((precision * recall) / (precision + recall))
  fb_measure <- (1 + b^2) * ((precision * recall) / ((b^2*precision) + recall))

  obj <- list(fpr = fpr, fnr = fnr, fdr = fdr, tpr = recall,
              sensitivity = recall, precision = precision, recall = recall,
              specificity = specificity, f1_measure = f1_measure,
              fb_measure = fb_measure)
  return(obj)
}


# Odds Ratio function
.compute_odds_ratio <- function(p1, p2) {
  return((p1/(1 - p1) ) / (p2 / (1 - p2)))
}

# Log odds Ratio function
.compute_log_odds_ratio <- function(p1, p2) {
  return(log(.compute_odds_ratio(p1, p2)))
}
