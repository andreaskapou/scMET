# Define ggplot2 theme for scatter plots
.scatter_theme <- function(legend_pos = "top") {
  p <- ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, size = rel(1.1),
                              margin = ggplot2::margin(0,0,2,0), color = "black"),
    legend.position = legend_pos,
    legend.title = ggplot2::element_blank(),
    axis.line = ggplot2::element_line(),
    #panel.border = element_blank(),
    #panel.grid.major = element_line(colour = "gainsboro"),
    panel.background = ggplot2::element_blank(),
    axis.text = ggplot2::element_text(color = "black", size = rel(0.8)),
    axis.title = ggplot2::element_text(color = "black", size = rel(1.2))
  )
  return(p)
}


# Get density of points in 2 dimensions.
# @param x A numeric vector.
# @param y A numeric vector.
# @param n Create a square n by n grid to compute density.
# @return The density within each square.
.get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}


#' @title Plot EFDR/EFNR grid
#'
#' @description Function for plotting the grid search performed to obtain the
#'   optimal posterior evidence threshold to achieve a specific EFDR.
#'
#' @param obj Either the scMET object after calling the
#'   \code{\link{scmet_hvf_lvf}} functions or the object from calling the
#'   \code{\link{scmet_differential}} function.
#' @param task String. When calling variable features, i.e. output of
#'   \code{\link{scmet_hvf_lvf}}, it can be either "hvf" or "lvf". For
#'   differential analysis, i.e. output of \code{\link{scmet_differential}}, it
#'   can be either: (1) "diff_mu" for diff mean methylation, (2) "diff_epsilon"
#'   for residual overdispersion, or (3) "diff_gamma" for overdispersion
#'   analysis.
#'
#' @return A ggplot2 object.
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @seealso \code{\link{scmet}}, \code{\link{scmet_differential}},
#'   \code{\link{scmet_hvf_lvf}}, \code{\link{scmet_plot_mean_var}},
#'   \code{\link{scmet_plot_vf_tail_prob}}, \code{\link{scmet_plot_volcano}},
#'   \code{\link{scmet_plot_ma}}
#'
#' @examples
#' # Fit scMET
#' obj <- scmet(Y = scmet_dt$Y, X = scmet_dt$X, L = 4, iter = 200)
#' obj <- scmet_hvf(scmet_obj = obj, delta_e = 0.7)
#' scmet_plot_vf_tail_prob(obj = obj, task = "hvf")
#'
#' @import ggplot2
#'
#' @export
scmet_plot_efdr_efnr_grid <- function(obj, task = "hvf") {
  if (!(methods::is(obj, "scmet_mcmc") ||
        methods::is(obj, "scmet_vb") ||
        methods::is(obj, "scmet_differential")) ) {
    stop("Object is not generated from scMET.")
  }
  # Object for HVF/LVF analysis
  if (methods::is(obj, "scmet_mcmc") || methods::is(obj, "scmet_vb")) {
    if (!(tolower(task) %in% c("hvf", "lvf"))) {
      stop("Task can be either 'hvf' or 'lvf'.\n")
    }
    evi_thresh_grid <- obj[[task]]$evidence_thresh_grid
    efdr_grid <- obj[[task]]$efdr_grid
    efnr_grid <- obj[[task]]$efnr_grid
    target_efdr <- obj[[task]]$target_efdr
    evi_thresh <- obj[[task]]$evidence_thresh
    title <- toupper(task)
  } else {
    if (task == "diff_mu") {
      mode <- "diff_mu_thresh"
      title <- "Differential mean"
    } else if (task == "diff_epsilon") {
      mode <- "diff_epsilon_thresh"
      title <- "Differential residual overdispersion"
    } else if (task == "diff_gamma") {
      mode <- "diff_gamma_thresh"
      title <- "Differential overdispersion"
    } else {
      stop("Task can be one of 'diff_mu', 'diff_epsilon' or 'diff_gamma'.\n")
    }
    evi_thresh_grid <- obj[[mode]]$evidence_thresh_grid
    efdr_grid <- obj[[mode]]$efdr_grid
    efnr_grid <- obj[[mode]]$efnr_grid
    target_efdr <- obj[[mode]]$target_efdr
    evi_thresh <- obj[[mode]]$evidence_thresh
  }
  gg <- ggplot() +
    geom_line(aes(evi_thresh_grid, efdr_grid, color = "EFDR")) +
    geom_line(aes(evi_thresh_grid, efnr_grid, color = "EFNR")) +
    geom_hline(aes(color = "Target EFDR", yintercept = target_efdr)) +
    geom_vline(aes(color = "Threshold", xintercept = evi_thresh)) +
    labs(x = "Posterior evidence threshold", y = "Error rate", title = title) +
    ylim(c(0,1)) + theme_bw() +
    .scatter_theme(legend_pos = "right")
  return(gg)
}



#' @title Plot tail probabilities for variable feature analysis
#'
#' @description Function for plotting the tail probabilities associated with the
#'   HVF/LVF analysis. The tail probabilities are plotted on the y-axis, and the
#'   user can choose which parameter can be plotted on the x-axis, using the `x`
#'   parameter.
#'
#' @param obj The scMET object after calling the \code{\link{scmet_hvf_lvf}}
#'   function.
#' @param x The parameter to plot on the x-axis. Values can be `mu` (default),
#'   `epsilon` or `gamma`.
#' @param task The task for identifying variable, either "hvf" or "lvf".
#' @param title Optional title, default NULL.
#' @param nfeatures Optional parameter, denoting a subset of number of features
#'   to plot (only for non HVF/LVF features). Mostly to reduce over-plotting.
#'
#' @return A ggplot2 object.
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @seealso \code{\link{scmet}}, \code{\link{scmet_differential}},
#'   \code{\link{scmet_hvf_lvf}}, \code{\link{scmet_plot_mean_var}},
#'   \code{\link{scmet_plot_efdr_efnr_grid}}, \code{\link{scmet_plot_volcano}},
#'   \code{\link{scmet_plot_ma}}
#'
#' @examples
#' # Fit scMET
#' obj <- scmet(Y = scmet_dt$Y, X = scmet_dt$X, L = 4, iter = 200)
#' obj <- scmet_hvf(scmet_obj = obj, delta_e = 0.7)
#' scmet_plot_vf_tail_prob(obj = obj, x = "mu")
#'
#' @export
#'
scmet_plot_vf_tail_prob <- function(obj, x = "mu", task = "hvf", title = NULL,
                                    nfeatures = NULL){
  if (!(methods::is(obj, "scmet_mcmc") ||
        methods::is(obj, "scmet_vb")) ) {
    stop("Object is not generated from scMET.")
  }
  is_variable <- NULL
  if (!(tolower(task) %in% c("hvf", "lvf"))) {
    stop("Task can be either 'hvf' or 'lvf'.\n")
  }
  if (length(obj[[tolower(task)]]) == 0) {
    stop("You should run HVF/LVF analysis prior to calling this function.\n")
  }
  ylab <- paste(toupper(task), "tail probability")
  if (x == "mu") {
    xlab <- expression(paste("Mean DNA methylation ", mu))
  } else if (x == "epsilon") {
    xlab <- expression(paste("Residual overdispersion ", epsilon))
  } else if (x == "gamma") {
    xlab <- expression(paste("Overdispersion ", gamma))
  } else {
    stop("x parameter can be one of 'mu', 'epsilon' or 'gamma'.\n")
  }

  up_task <- toupper(task)
  size <- c(1.6, 0.6)
  names(size) <- c(up_task, "Other")
  fill <- c("red", "gray80")
  names(fill) <- c(up_task, "Other")
  alpha <- c(0.65, 0.45)
  names(alpha) <- c(task, "Other")

  df <- obj[[tolower(task)]]$summary
  # Subset number of features
  if (!is.null(nfeatures)) {
    assertthat::assert_that(nfeatures > 0)
    # Features that are non HVF/LVF
    tmp <- df[df$is_variable == FALSE, ]
    if (NROW(tmp) > nfeatures) {
      tmp <- tmp[sample(NROW(tmp), nfeatures), ]
      df <- rbind(df[df$is_variable == TRUE, ], tmp)
    }
  }

  gg <- ggplot(df,aes_string(x = x, y = "tail_prob")) +
    geom_point(aes(fill = ifelse(is_variable, up_task, "Other"),
                   size = ifelse(is_variable, up_task, "Other"),
                   alpha = ifelse(is_variable, up_task, "Other")),
               colour = "black", shape = 21, stroke = 0.03) +
    scale_fill_manual(values = fill) +
    scale_size_manual(values = size) +
    scale_alpha_manual(values = alpha) +
    geom_hline(yintercept = obj[[tolower(task)]]$evidence_thresh[1], lty = 2,
               col = "black") +
    labs(x = xlab, y = ylab, title = title) + theme_classic() +
    .scatter_theme(legend_pos = "right")
  return(gg)
}


#' @title Plotting mean-variability relationship
#'
#' @description Function for plotting mean methylation on x-axis and variability
#'   on y-axis (either overdispersion or residual overdispersion). If HVF/LVF
#'   analysis is performed, points will be also coloured accordingly.
#'
#' @param obj The scMET object after calling the \code{\link{scmet_hvf_lvf}}
#'   function.
#' @param y The parameter to plot on the y-axis. Values can be `gamma` (default)
#'   or `epsilon`.
#' @param task If NULL (default) the mean-variability relationship is plotted.
#'   If set to "hvf" or "lvf", points are coloured according the HVF/LVF
#'   analysis task.
#' @param show_fit Logical, whether to show the fitted mean-overdispersion
#'   trend. Applicable only when `y = gamma` and `task = NULL`.
#' @param title Optional title, default NULL.
#' @param nfeatures Optional parameter, denoting a subset of number of features
#'   to plot. Mostly to reduce over-plotting. When `task = hvf or lvf`, the
#'   subsampling is performed on the features that are not called as HVF or LVF
#'   (i.e. not interesting features).
#' @param n Optional integer denoting the number of grid points to colour them
#'   by density. Used by \code{\link[MASS]{kde2d}} function. Used only when
#'   `task = NULL`.
#'
#' @return A ggplot2 object.
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @seealso \code{\link{scmet}}, \code{\link{scmet_differential}},
#'   \code{\link{scmet_hvf_lvf}}, \code{\link{scmet_plot_vf_tail_prob}},
#'   \code{\link{scmet_plot_efdr_efnr_grid}}, \code{\link{scmet_plot_volcano}},
#'   \code{\link{scmet_plot_ma}}
#'
#' @examples
#' # Fit scMET
#' obj <- scmet(Y = scmet_dt$Y, X = scmet_dt$X, L = 4, iter = 200)
#' scmet_plot_mean_var(obj = obj, y = "gamma")
#'
#' @export
#'
scmet_plot_mean_var <- function(obj, y = "gamma", task = NULL, show_fit = TRUE,
                                title = NULL, nfeatures = NULL, n = 80) {
  if (!(methods::is(obj, "scmet_mcmc") ||
        methods::is(obj, "scmet_vb")) ) {
    stop("Object is not generated from scMET.")
  }
  xlab <- expression(paste("Mean DNA methylation ", mu))
  if (y == "epsilon") {
    ylab <- expression(paste("Residual overdispersion ", epsilon))
  } else if (y == "gamma") {
    ylab <- expression(paste("Overdispersion ", gamma))
  } else {
    stop("y parameter can be one of 'epsilon' or 'gamma'.\n")
  }

  is_variable = mu = gamma = epsilon = density_gamma = density_epsilon <- NULL

  # Just plot mean-variability relationship, no HVF/LVF analysis
  if (is.null(task)) {
    # Obtain posterior median parameter estimates
    df <- data.table(mu = matrixStats::colMedians(obj$posterior$mu),
                     epsilon = matrixStats::colMedians(obj$posterior$epsilon),
                     gamma = matrixStats::colMedians(obj$posterior$gamma))
    df <- df %>%
      .[, density_gamma := .get_density(mu, gamma, n = n)] %>%
      .[, density_epsilon := .get_density(mu, epsilon, n = n)]
    col_grad <- ifelse(y == "gamma", "density_gamma", "density_epsilon")

    # Subset number of features
    if (!is.null(nfeatures)) {
      assertthat::assert_that(nfeatures > 0)
      if (NROW(df) > nfeatures) {
        df <- df[sample(NROW(df), nfeatures), ]
      }
    }

    gg <- ggplot(df, aes_string(x = "mu", y = y, color = col_grad)) +
      geom_point(size = 1) +
      viridis::scale_fill_viridis() +
      viridis::scale_color_viridis() + theme_classic() +
      labs(x = xlab, y = ylab, title = title) +
      .scatter_theme(legend_pos = "none")

    # Show RBF fit
    if (show_fit) {
      xs <- seq(min(df$mu), max(df$mu), length.out = 200)
      hs <- create_design_matrix(L = obj$opts$L, X = xs, c = obj$opts$rbf_c)
      ys <- stats::plogis(c(hs %*% matrixStats::colMedians(obj$posterior$w_gamma)))
      if (y == "epsilon") {
        gg <- gg +
          geom_hline(yintercept = 0, linetype = "dashed", color = "gray20",
                     size = 1.25)
      } else {
        gg <- gg +
          geom_line(data = data.frame(xs = xs, ys = ys), aes(x = xs, y = ys),
                    color = "gray20", size = 1.35, alpha = 1)
      }
    }
  } else if (!(tolower(task) %in% c("hvf", "lvf"))) {
    stop("Task can be either 'hvf' or 'lvf'.\n")
  } else {
    if (length(obj[[tolower(task)]]) == 0) {
      stop("HVF/LVF analysis is not performed.\n",
           "Either identify variable features or set `task = NULL` ",
           "to just plot mean-variability relationship.\n")
    }

    task <- toupper(task)
    size <- c(1.6, 0.65)
    names(size) <- c(task, "Other")
    fill <- c("red", "gray80")
    names(fill) <- c(task, "Other")
    alpha <- c(0.65, 0.45)
    names(alpha) <- c(task, "Other")

    df <- obj[[tolower(task)]]$summary
    # Subset number of features
    if (!is.null(nfeatures)) {
      assertthat::assert_that(nfeatures > 0)
      # Features that are non HVF/LVF
      tmp <- df[df$is_variable == FALSE, ]
      if (NROW(tmp) > nfeatures) {
        tmp <- tmp[sample(NROW(tmp), nfeatures), ]
        df <- rbind(df[df$is_variable == TRUE, ], tmp)
      }
    }

    gg <- ggplot(df, aes_string(x = "mu", y = y)) +
      geom_point(aes(fill = ifelse(is_variable, task, "Other"),
                     size = ifelse(is_variable, task, "Other"),
                     alpha = ifelse(is_variable, task, "Other")),
                 colour = "black", shape = 21, stroke = 0.03) +
      scale_fill_manual(values = fill) +
      scale_size_manual(values = size) +
      scale_alpha_manual(values = alpha) +
      labs(x = xlab, y = ylab, title = title) + theme_classic() +
      .scatter_theme(legend_pos = "right")
  }
  return(gg)
}


#' @title Plot true versus inferred parameter estimated.
#'
#' @description Function for plotting true on x-axis and inferred parameter
#'   estimates on y-axis (either mean methylation or overdispersion). Along with
#'   posterior medians, the 80 high posterior density is shown as error bars.
#'   Wehn MLE estimates are provided, a plot showing the shrinkage introduced by
#'   scMET is shown as arrows.
#'
#' @param obj The scMET object after calling the \code{\link{scmet}} function.
#' @param sim_dt The simulated data object. E.g. after calling the
#'   \code{\link{scmet_simulate}} function.
#' @param param The parameter to plot posterior estimates, either "mu" or
#'   "gamma".
#' @param mle_fit A three column matrix of beta-binomial maximum likelihood
#'   estimates. First column feature name, second column mean methylation and
#'   third column overdispersion estimates. Number of features should match the
#'   ones used by scMET.
#' @param diff_feat_idx Vector with locations of features that were simulated to
#'   be differentially variable or methylated. This is stored in the object
#'   after calling the \code{\link{scmet_simulate_diff}} function.
#' @param hpd_thresh The high posterior density threshold, as computed by the
#'   \code{\link[coda]{HPDinterval}} function.
#' @param title Optional title, default NULL.
#' @param nfeatures Optional parameter, denoting a subset of number of features
#'   to plot. Mostly to reduce over-plotting.
#'
#' @return A ggplot2 object.
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @seealso \code{\link{scmet}}, \code{\link{scmet_simulate_diff}},
#'   \code{\link{scmet_simulate}}, \code{\link{scmet_plot_mean_var}},
#'   \code{\link{scmet_plot_vf_tail_prob}},
#'   \code{\link{scmet_plot_efdr_efnr_grid}}, \code{\link{scmet_plot_volcano}},
#'   \code{\link{scmet_plot_ma}}
#'
#' @examples
#' # Fit scMET
#' obj <- scmet(Y = scmet_dt$Y, X = scmet_dt$X, L = 4, iter = 200)
#' scmet_plot_estimated_vs_true(obj = obj, sim_dt = scmet_dt, param = "mu")
#'
#' # BB MLE fit to compare with scMET
#' mle_fit <- scmet_dt$Y[, bb_mle(cbind(total_reads, met_reads))
#' [c("mu", "gamma")], by = c("Feature")]
#' scmet_plot_estimated_vs_true(obj = obj, sim_dt = scmet_dt, param = "mu",
#' mle_fit = mle_fit)
#'
#' @export
#'
scmet_plot_estimated_vs_true <- function(obj, sim_dt, param = "mu",
                                         mle_fit = NULL, diff_feat_idx = NULL,
                                         hpd_thresh = 0.8, title = NULL,
                                         nfeatures = NULL) {
  if (!(methods::is(obj, "scmet_mcmc") ||
        methods::is(obj, "scmet_vb")) ) {
    stop("Object is not generated from scMET.")
  }
  if (!(tolower(param) %in% c("mu", "gamma"))) {
    stop("Param can be either 'mu' or 'gamma'.\n")
  }
  suffix <- ifelse(tolower(param) == "mu", "mean DNA methylation", "overdispersion")
  true = scmet = hpdlow = hpdhigh = bbmle <- NULL
  res <- data.frame(true = sim_dt$theta_true[[param]],
                    scmet = matrixStats::colMedians(obj$posterior[[param]]),
                    hpdlow = coda::HPDinterval(coda::mcmc(obj$posterior[[param]]),
                                               prob = hpd_thresh)[, 1],
                    hpdhigh = coda::HPDinterval(coda::mcmc(obj$posterior[[param]]),
                                              prob = hpd_thresh)[, 2])

  if (!is.null(mle_fit)) {
    if (NCOL(mle_fit) != 3) {
      stop("The 'mle_fit' parameter should have three columns:",
           "'Feature', 'mu' and 'gamma'.\n")
    }
    if (NROW(mle_fit) != length(obj$feature_names)) {
      stop("Number of features between scMET and BB-MLE does not match!\n")
    }
    if (!is.data.frame(mle_fit)) {
      stop("The 'mle_fit' should be a data.frame or data.table object.\n")
    }
    mle_fit <- as.data.frame(mle_fit)
    mle_fit <- mle_fit[match(obj$feature_names, mle_fit[, 1]), ]
    # Extract the correct column
    idx <- ifelse(tolower(param) == "mu", 2, 3)
    res$bbmle <- mle_fit[, idx]
    # Subset number of features
    if (!is.null(nfeatures)) {
      assertthat::assert_that(nfeatures > 0)
      if (NROW(res) > nfeatures) {
        res <- res[sample(NROW(res), nfeatures), ]
      }
    }
    colors <- c("BB MLE" = "#999999", "scMET" = "#E69F00")
    gg <-  ggplot(res) +
      geom_point(aes(y = bbmle, x = true, fill = "BB MLE"), shape = 21,
                 size = 1.7, alpha = 0.7, stroke = 0.2) +
      geom_point(aes(y = scmet, x = true, fill = "scMET"), shape = 21,
                 size = 2.2, stroke = 0.2) +
      geom_abline(intercept = 0, slope = 1, color = "black",
                  linetype = "dashed", alpha = 0.7) +
      geom_segment(aes(x = true, y = bbmle, xend = true, yend = scmet),
                   arrow = arrow(angle = 10, length = unit(0.06, "inches")),
                   size = 0.14, alpha = 0.5) +
      labs(x = paste("True", suffix), y = paste("Estimated", suffix),
           color = "Legend") +
      scale_fill_manual(values = colors) + theme_classic() +
      .scatter_theme(legend_pos = "right")
  } else {
    # Colour differential features
    if (!is.null(diff_feat_idx)) {
      res$diff <- FALSE
      res$diff[diff_feat_idx] <- TRUE
    }
    # Subset number of features
    if (!is.null(nfeatures) && is.null(diff_feat_idx)) {
      assertthat::assert_that(nfeatures > 0)
      if (NROW(res) > nfeatures) {
        res <- res[sample(NROW(res), nfeatures), ]
      }
    }
    # Create plot
    gg <- ggplot(res, aes(x = true, y = scmet)) +
      geom_errorbar(aes(ymin = hpdlow, ymax = hpdhigh),
                    width = .003, position = position_dodge(0.05), size = 0.25) +
      xlab(paste("True", suffix)) + ylab(paste("Estimated", suffix)) +
      ggtitle(title)

    if (!is.null(diff_feat_idx)) {
      gg <- gg +
        geom_point(aes(color = diff), size = 1) +
        scale_color_manual(values = c('#595959', 'red'))
    } else {
      gg <- gg + geom_point(size = 1)
    }
    # Add theme
    gg <- gg +
      geom_abline(intercept = 0, slope = 1, color = "orange", linetype = "dashed") +
      theme_classic() + .scatter_theme(legend_pos = "none")
  }
  return(gg)
}


#' @title Volcano plot for differential analysis
#'
#' @description Function showing volcano plots for differential analysis. The
#'   posterior tail probabilities are ploteted on the y-axis, and depending on
#'   the differential test to plot the effect size will be plotted on the
#'   x-axis. For differential variability (DV) analysis we recommend using the
#'   `epsilon` parameter.
#'
#' @param diff_obj The differential scMET object after calling the
#'   \code{\link{scmet_differential}} function.
#' @param task The differential test to plot. For differential mean methylation:
#'   `diff_mu` that plots the LOR(mu_A, mu_B) on x-axis. For differential
#'   variability: either (1) `diff_epsilon` that plots the change (epsilon_A -
#'   epsilon_B), or (2) `diff_gamma` that plots the LOR(gamma_A, gamma_B) on
#'   x-axis.
#' @param xlab Optional x-axis label.
#' @param ylab Optional y-axis label.
#' @param title Optional title, default NULL.
#' @param nfeatures Optional parameter, denoting a subset of number of features
#'   to plot (only for non-differential features). Mostly to reduce
#'   over-plotting.
#'
#' @return A ggplot2 object.
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @seealso \code{\link{scmet}}, \code{\link{scmet_differential}},
#'   \code{\link{scmet_hvf_lvf}}, \code{\link{scmet_plot_mean_var}},
#'   \code{\link{scmet_plot_vf_tail_prob}},
#'   \code{\link{scmet_plot_efdr_efnr_grid}}, \code{\link{scmet_plot_ma}}
#'
#' @examples
#' # Fit scMET for each group
#' fit_A <- scmet(Y = scmet_diff_dt$scmet_dt_A$Y,
#' X = scmet_diff_dt$scmet_dt_A$X, L = 4, iter = 200, seed = 12)
#' fit_B <- scmet(Y = scmet_diff_dt$scmet_dt_B$Y,
#' X = scmet_diff_dt$scmet_dt_B$X, L = 4, iter = 200, seed = 12)
#'
#' # Run differential test
#' diff_obj <- scmet_differential(obj_A = fit_A, obj_B = fit_B)
#' # Create volcano plot
#' scmet_plot_volcano(diff_obj, task = "diff_epsilon")
#'
#' @export
#'
scmet_plot_volcano <- function(diff_obj, task = "diff_epsilon", xlab = NULL,
                               ylab = "Posterior tail probability",
                               title = NULL, nfeatures = NULL) {
  if (!(methods::is(diff_obj, "scmet_differential") )) {
    stop("Object is not generated from scMET differential function.")
  }
  task <- tolower(task)
  if (!(task %in% c("diff_mu", "diff_epsilon", "diff_gamma"))) {
    stop("The `mode` param can be either 'diff_mu', 'diff_epsilon' ",
         "or 'diff_gamma'.\n")
  }

  # Extract corresponding matrices
  sum_obj <- diff_obj[[paste0(task, "_summary")]]
  thresh_obj <- diff_obj[[paste0(task, "_thresh")]]
  # Define the corresponding info to extract
  if (task == "diff_mu") {
    x <- "mu_LOR"
    y <- "mu_tail_prob"
    psi <- diff_obj$opts$psi_m
    test <- "mu_diff_test"
    if (is.null(xlab)) {
      xlab <- expression(paste("Change in mean: LOR(", mu[A], ", ", mu[B], ")"))
    }

  } else if (task == "diff_epsilon") {
    x <- "epsilon_change"
    y <- "epsilon_tail_prob"
    psi <- diff_obj$opts$psi_e
    test <- "epsilon_diff_test"
    if (is.null(xlab)) {
      xlab <- expression(paste("Change in variability: ( ",
                               epsilon[A] - epsilon[B], " )"))
    }
  } else {
    x <- "gamma_LOR"
    y <- "gamma_tail_prob"
    psi <- diff_obj$opts$psi_g
    test <- "gamma_diff_test"
    if (is.null(xlab)) {
      xlab <- expression(paste("Change in variability: LOR(",
                               gamma[A], ", ", gamma[B], ")"))
    }
  }

  # Store the names of the two groups
  diff_names <- c(diff_obj$opts$group_label_A, diff_obj$opts$group_label_B)
  diff_names <- paste0(diff_names, "+")
  # Set all other outcomes to NoDiff
  idx <- which(!sum_obj[[test]] %in% diff_names)
  sum_obj[[test]][idx] <- "NoDiff"

  size <- c(1.4, 1.4, 0.7)
  names(size) <- c(diff_names, "NoDiff")
  alpha <- c(0.6, 0.6, 0.4)
  names(alpha) <- c(diff_names, "NoDiff")

  # Subset number of features
  if (!is.null(nfeatures)) {
    assertthat::assert_that(nfeatures > 0)
    # Features that are non HVF/LVF
    tmp <- sum_obj[sum_obj[[test]] == "NoDiff", ]
    if (NROW(tmp) > nfeatures) {
      tmp <- tmp[sample(NROW(tmp), nfeatures), ]
      sum_obj <- rbind(sum_obj[sum_obj[[test]] != "NoDiff", ], tmp)
    }
  }

  # Create ggplot
  gg <- ggplot(data = sum_obj, aes_string(x = x, y = y)) +
    geom_point(aes_string(fill = test, size = test, alpha = test),
               shape = 21, stroke = 0.02) +
    scale_size_manual(name = "Hits", values = size) +
    scale_alpha_manual(name = "Hits", values = alpha) +
    scale_fill_manual(name = "Hits",
                      values = c("lightpink3", "darkolivegreen3", "#595959")) +
    geom_vline(xintercept = c(-psi, psi), color = "dodgerblue4",
               linetype = "dashed", alpha = 0.8) +
    geom_hline(yintercept = thresh_obj$evidence_thresh, color = "orange",
               linetype = "dashed", alpha = 0.8) +
    labs(x = xlab, y = ylab, title = title) + theme_classic() +
    .scatter_theme(legend_pos = "right")
  return(gg)
}



#' @title MA plot for differential analysis
#'
#' @description Function showing MA plots for differential analysis. The y-axis
#'   shows difference between measurements across two groups and the x-axis
#'   shows the average measurements across the two groups.
#'
#' @param diff_obj The differential scMET object after calling the
#'   \code{\link{scmet_differential}} function.
#' @param task The differential test to plot. For differential mean methylation:
#'   `diff_mu` that plots the LOR(mu_A, mu_B) on y-axis. For differential
#'   variability: either (1) `diff_epsilon` that plots the change (epsilon_A -
#'   epsilon_B), or (2) `diff_gamma` that plots the LOR(gamma_A, gamma_B) on
#'   y-axis.
#' @param x The average parameter across the two populations to plot on the
#'   x-axis. Can be either `mu`, `epsilon` or `gamma`. When `task = epsilon`, x
#'   can be either `mu` or `epsilon`. When `task = gamma`, x can be either `mu`
#'   or `gamma`. When `task = mu`, x can be only `mu`.
#' @param xlab Optional x-axis label.
#' @param ylab Optional y-axis label.
#' @param title Optional title, default NULL.
#' @param nfeatures Optional parameter, denoting a subset of number of features
#'   to plot (only for non-differential features). Mostly to reduce
#'   over-plotting.
#'
#' @return A ggplot2 object.
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @seealso \code{\link{scmet}}, \code{\link{scmet_differential}},
#'   \code{\link{scmet_hvf_lvf}}, \code{\link{scmet_plot_mean_var}},
#'   \code{\link{scmet_plot_vf_tail_prob}},
#'   \code{\link{scmet_plot_efdr_efnr_grid}}, \code{\link{scmet_plot_volcano}}
#'
#' @examples
#' # Fit scMET for each group
#' fit_A <- scmet(Y = scmet_diff_dt$scmet_dt_A$Y,
#' X = scmet_diff_dt$scmet_dt_A$X, L = 4, iter = 200, seed = 12)
#' fit_B <- scmet(Y = scmet_diff_dt$scmet_dt_B$Y,
#' X = scmet_diff_dt$scmet_dt_B$X, L = 4, iter = 200, seed = 12)
#'
#' # Run differential test
#' diff_obj <- scmet_differential(obj_A = fit_A, obj_B = fit_B)
#' # Create volcano plot
#' scmet_plot_ma(diff_obj, task = "diff_epsilon")
#'
#' @export
#'
scmet_plot_ma <- function(diff_obj, task = "diff_epsilon", x = "mu",
                          xlab = NULL, ylab = NULL, title = NULL,
                          nfeatures = NULL) {
  if (!(methods::is(diff_obj, "scmet_differential") )) {
    stop("Object is not generated from scMET differential function.")
  }
  task <- tolower(task)
  if (!(task %in% c("diff_mu", "diff_epsilon", "diff_gamma"))) {
    stop("The `mode` param can be either 'diff_mu', 'diff_epsilon' ",
         "or 'diff_gamma'.\n")
  }

  # Extract corresponding matrices
  sum_obj <- diff_obj[[paste0(task, "_summary")]]
  # Define the corresponding info to extract
  if (task == "diff_mu") {
    x_param <- "mu_overall"
    y <- "mu_LOR"
    psi <- diff_obj$opts$psi_m
    test <- "mu_diff_test"
    if (is.null(ylab)) {
      ylab <- expression(paste("Change in mean: LOR(", mu[A], ", ", mu[B], ")"))
    }
    if (is.null(xlab)) {
      xlab <- expression(paste("Overall mean methylation ", mu))
    }
  } else if (task == "diff_epsilon") {
    y <- "epsilon_change"
    if (x == "epsilon") {
      x_param <- "epsilon_overall"
      if (is.null(xlab)) {
        xlab <- expression(paste("Overall residual overdispersion ", epsilon))
      }
    } else if (x == "mu") {
      x_param <- "mu_overall"
      if (is.null(xlab)) {
        xlab <- expression(paste("Overall mean methylation ", mu))
      }
    } else {
      stop("For `task = diff_epsilon`, `x` can either be `mu` or `epsilon`\n")
    }
    psi <- diff_obj$opts$psi_e
    test <- "epsilon_diff_test"
    if (is.null(ylab)) {
      ylab <- expression(paste("Change in variability: ( ",
                               epsilon[A] - epsilon[B], " )"))
    }
  } else {
    y <- "gamma_LOR"
    if (x == "gamma") {
      x_param <- "gamma_overall"
      if (is.null(xlab)) {
        xlab <- expression(paste("Overall overdispersion ", gamma))
      }
    } else if (x == "mu") {
      x_param <- "mu_overall"
      if (is.null(xlab)) {
        xlab <- expression(paste("Overall mean methylation ", mu))
      }
    } else {
      stop("For `task = diff_gamma`, `x` can either be `mu` or `gamma`\n")
    }
    psi <- diff_obj$opts$psi_g
    test <- "gamma_diff_test"
    if (is.null(ylab)) {
      ylab <- expression(paste("Change in variability: LOR(",
                               gamma[A], ", ", gamma[B], ")"))
    }
  }

  # Store the names of the two groups
  diff_names <- c(diff_obj$opts$group_label_A, diff_obj$opts$group_label_B)
  diff_names <- paste0(diff_names, "+")
  # Set all other outcomes to NoDiff
  idx <- which(!sum_obj[[test]] %in% diff_names)
  sum_obj[[test]][idx] <- "NoDiff"

  size <- c(1.4, 1.4, 0.7)
  names(size) <- c(diff_names, "NoDiff")
  alpha <- c(0.6, 0.6, 0.4)
  names(alpha) <- c(diff_names, "NoDiff")

  # Subset number of features
  if (!is.null(nfeatures)) {
    assertthat::assert_that(nfeatures > 0)
    # Features that are non HVF/LVF
    tmp <- sum_obj[sum_obj[[test]] == "NoDiff", ]
    if (NROW(tmp) > nfeatures) {
      tmp <- tmp[sample(NROW(tmp), nfeatures), ]
      sum_obj <- rbind(sum_obj[sum_obj[[test]] != "NoDiff", ], tmp)
    }
  }

  # Create ggplot
  gg <- ggplot(data = sum_obj, aes_string(x = x_param, y = y)) +
    geom_point(aes_string(fill = test, size = test, alpha = test),
               shape = 21, stroke = 0.02) +
    scale_size_manual(name = "Hits", values = size) +
    scale_alpha_manual(name = "Hits", values = alpha) +
    scale_fill_manual(name = "Hits",
                      values = c("lightpink3", "darkolivegreen3", "#595959")) +
    geom_hline(yintercept = c(-psi, psi), color = "dodgerblue4",
               linetype = "dashed", alpha = 0.8) +
    labs(x = xlab, y = ylab, title = title) + theme_classic() +
    .scatter_theme(legend_pos = "right")
  return(gg)
}
