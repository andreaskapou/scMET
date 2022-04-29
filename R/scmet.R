#' @name scmet
#' @aliases scMET scmet_vb scmet_mcmc
#'
#' @title Perform inference with scMET
#'
#' @description Compute posterior of scMET model. This is the main function
#'   which infers model parameters and corrects for the mean-overdispersion
#'   relationship. The most important parameters the user should focus are `X`,
#'   `L`, `user_mcmc` and `iter`. Advanced users may want to optimise the model
#'   by changing the prior parameters. For small datasets, we recommend using
#'   MCMC implementation of scMET since it is more stable.
#'
#' @param Y Observed data (methylated reads and total reads) for each feature
#'   and cell, in a long format \code{\link{data.table}}. That is it should have
#'   4 named columns: (Feature, Cell, total_reads, met_reads).
#' @param X Covariates which might explain variability in mean (methylation). If
#'   X = NULL, then we do not perform any correction on the mean estimates. NOTE
#'   that if X is provided, `rownames` of X should be the unique feature names
#'   in Y. If the dimensions or all feature names do not match, an error will
#'   be thrown.
#' @param L Total number of basis function to fit the mean-overdispersion trend.
#'   For L = 1, this reduces to a model that does not correct for the
#'   mean-overdispersion relationship.
#' @param use_mcmc Logical, whether to use the MCMC implementation for posterior
#'   inference. If FALSE, we run the VB implementation (default). For small
#'   datasets, we recommend using MCMC implementation since it is more stable.
#' @param use_eb Logical, whether to use 'Empirical Bayes' for parameter
#'   initialization. If `TRUE` (default), it will intialise the `m_wmu` and
#'   `m_wgamma` parameters below.
#' @param iter Total number of iterations, either MCMC or VB algorithm.
#' @param algorithm Stan algorithm to be used by Stan. If MCMC: Possible values
#'   are: "NUTS", "HMC". If VB: Possible values are: "meanfield" and "fullrank".
#' @param output_samples If VB algorithm, the number of posterior samples to
#'   draw and save.
#' @param chains Total number of chains.
#' @param m_wmu Prior mean of regression coefficients for covariates X.
#' @param s_wmu Prior standard deviation of regression coefficients for
#'   covariates X.
#' @param s_mu Prior standard deviation for mean parameter `mu`.
#' @param m_wgamma Prior mean of regression coefficients of the basis functions.
#' @param s_wgamma Prior standard deviation of regression coefficients of the
#'   basis functions.
#' @param a_sgamma Gamma prior (shape) for standard deviation for dispersion
#'   parameter `gamma`.
#' @param b_sgamma Gamma prior (rate) for standard deviation for dispersion
#'   parameter `gamma`.
#' @param rbf_c Scale parameter for empirically computing the variance of the
#'   RBFs.
#' @param init_using_eb Logical, initial values of parameters for STAN posterior
#'   inference. Preferably this should be set always to TRUE, to lower the
#'   chances of VB/MCMC initialisations being far away from posterior mass.
#' @param tol_rel_obj If VB algorithm, the convergence tolerance on the relative
#'   norm of the objective.
#' @param n_cores Total number of cores.
#' @param lambda The penalty term to fit the RBF coefficients for the
#'   mean-overdispersion trend when initialising hyper-parameter with EB.
#' @param seed The seed for random number generation.
#' @param ... Additional parameters passed to \code{Stan} fitting functions.
#'
#' @return An object of class \code{scmet_mcmc} or \code{scmet_vb} with the
#'   following elements: \itemize{ \item{ \code{posterior}: A list of matrices
#'   containing the samples from the posterior. Each matrix corresponds to a
#'   different parameter returned from scMET.} \item{ \code{Y}: The observed
#'   data Y. } \item{ \code{feature_names}: A vector of feature names.} \item{
#'   \code{theta_priors}: A list with all prior parameter values, for
#'   reproducibility purposes. } \item{\code{opts}: A list of all additional
#'   parameters when running scMET. For reproducibility purposes.}}
#'
#' @seealso \code{\link{scmet_differential}}, \code{\link{scmet_hvf_lvf}}
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @examples
#' # Fit scMET (in practice 'iter' should be much larger)
#' obj <- scmet(Y = scmet_dt$Y, X = scmet_dt$X, L = 4, iter = 300)
#'
#' @importFrom stats qlogis
#' @export
scmet <- function(Y, X = NULL, L = 4, use_mcmc = FALSE, use_eb = TRUE,
                  iter = 5000, algorithm = "meanfield", output_samples = 2000,
                  chains = 4, m_wmu = rep(0,NCOL(X)), s_wmu = 2, s_mu = 1.5,
                  m_wgamma = rep(0, L), s_wgamma = 2, a_sgamma = 2,
                  b_sgamma = 3, rbf_c = 1, init_using_eb = TRUE,
                  tol_rel_obj = 1e-4, n_cores = 2, lambda = 4,
                  seed = sample.int(.Machine$integer.max, 1), ...) {
  # So RMD check does not complain
  init_vals = Feature = total_reads = met_reads <- NULL

  #-------------------
  # Checking parameter intialization
  #-------------------
  # Check that Y has 4 columns and ensure it is a data.table
  assertthat::assert_that(NCOL(Y) == 4)
  if (!is.data.table(Y)) {
    message("Converting Y matrix to data.table.\n")
    Y <- as.data.table(Y)
    colnames(Y) <- c("Feature", "Cell", "total_reads", "met_reads")
  }
  if (chains < 1) { stop("You should have at least one chain.") }
  # If no X covariates, add a dummy column of ones
  if (is.null(X)) {
    X <- matrix(1, nrow = length(unique(Y$Feature)), ncol = 1)
    rownames(X) <- unique(Y$Feature)
  } else {
    if (NROW(X) != length(unique(Y$Feature))) {
      stop("Number of X covariates does not match number of features.")
    } else if (is.null(rownames(X))) {
      stop("X should have feture names as rownames(X).")
    } else if (!all(rownames(X) %in% unique(Y$Feature)) ) {
      stop("Rownames in X do not match feature names in Y.")
    }
  }
  # Check that dimensions match
  if ( NCOL(X) != length(m_wmu) ) {
    stop("Number of covariates X does not match length of coefficients.")
  }
  # Avoid issue with Stan considering 1-dim vector as vector and not real.
  if (NCOL(X) == 1) { m_wmu <- as.array(m_wmu) }

  #-------------------
  # Tests for correct parameter setting for mean-overdispersion trend
  #-------------------
  if (length(m_wgamma) != L) {
    stop("Number of RBFs should match length of `m_wgamma` parameter.")
  }
  # Avoid issue with Stan considering 1-dim vector as vector and not real
  if (L == 1) { m_wgamma <- as.array(m_wgamma) }

  ## Order by Feature and cell id, so we can the segment function inside stan.
  cat("Sorting features.\n")
  Y <- Y[order(Feature), , drop = FALSE]
  X <- as.matrix(X[order(rownames(X)), , drop = FALSE])

  cat("Total # of cells: ", length(unique(Y$Cell)), ".\n")
  cat("Total # of features: ", length(unique(Y$Feature)), ".\n")

  # Total number of observed cells in each feature
  N_cells <- Y[, .N, by = c("Feature")]$N
  # Do we have enough cells for each feature to perform downstream analysis
  if (!all(N_cells > 3) ) {
    stop("Each feature should have at least 3 cells, to perform inference.\n",
         "Perform a different filtering or remove those features!\n")
  }

  if (use_eb) {
    cat(date(), ": Using EB to set model priors.\n")

    # For efficiency, we should subsample features, since we just want a broad prior.
    # Also, keep only those features for which we could compute the MLE?

    # For each feature compute MLE or MM estimates of BB
    bb_mle_fit <- Y[, bb_mle(cbind(total_reads, met_reads))[c("gamma", "mu")],
                    by = c("Feature")]
    # Check outlier values and put a threshold
    idx <- bb_mle_fit$gamma < 1e-3
    if (sum(idx > 0)) {
      bb_mle_fit$gamma[idx] <- 1e-3 + stats::runif(n = sum(idx), min = 0,
                                                   max = 1e-4)
    }
    idx <- bb_mle_fit$gamma > 1 - 1e-3
    if (sum(idx > 0)) {
      bb_mle_fit$gamma[idx] <- 1 - 1e-3 - stats::runif(n = sum(idx), min = 0,
                                                max = 1e-4)
    }
    idx <- bb_mle_fit$mu < 1e-3
    if (sum(idx > 0)) {
      bb_mle_fit$mu[idx] <- 1e-3 + stats::runif(n = sum(idx), min = 0,
                                                max = 1e-4)
    }
    idx <- bb_mle_fit$mu > 1 - 1e-3
    if (sum(idx > 0)) {
      bb_mle_fit$mu[idx] <- 1 - 1e-3 - stats::runif(n = sum(idx), min = 0,
                                                    max = 1e-4)
    }

    # Compute prior for mean parameters from the data, empirical Bayes.
    df <- data.frame(X = X, y = stats::qlogis(bb_mle_fit$mu))
    m_wmu <- stats::coef(stats::lm(formula = y ~ X + 0, data = df))
    names(m_wmu) <- NULL
    if (NCOL(X) == 1) { m_wmu <- as.array(m_wmu) }

    # Compute prior for dispersion parameters from the data, empirical Bayes.
    H <- create_design_matrix(L = L, X = bb_mle_fit$mu, c = rbf_c)
    m_wgamma <- .lm_mle_penalized(y = stats::qlogis(bb_mle_fit$gamma),
                                  H = H, lambda = lambda)
    if (L == 1) { m_wgamma <- as.array(m_wgamma) }

    # If we also require to set initial values based on EB estimates
    if (init_using_eb) {
      cat("Usin EB estimates as initial values for posterior inference.\n")
      # Function form
      init_func <- function() {
        list(w_mu = m_wmu, w_gamma = m_wgamma,
             logit_mu = stats::qlogis(bb_mle_fit$mu),
             logit_gamma = stats::qlogis(bb_mle_fit$gamma))
      }
      # generate a list of lists to specify initial values
      init_vals <- lapply(seq_len(chains), function(id) init_func())
    }
    cat(date(), ": Finished EB.\n")
  }

  # To pass RMD check
  N = J = N_X = y = n = C <- NULL
  # Create data object for Stan
  dt <- within(list(), {
    N <- NROW(Y)
    J <- length(N_cells)
    N_X <- NCOL(X)
    L <- L
    y <- Y[, met_reads]
    n <- Y[, total_reads]
    C <- N_cells
    rbf_c <- rbf_c
    X <- X
    m_wmu <- m_wmu
    s_wmu <- s_wmu
    s_mu <- s_mu
    m_wgamma <- m_wgamma
    s_wgamma <- s_wgamma
    a_sgamma <- a_sgamma
    b_sgamma <- b_sgamma
  })

  # Perform inference
  cat(date(), ": Starting inference.\n")
  #-----------------
  # Handle possible errors when Stan initialization fails during inference.
  # Maybe should write a better code for this, at some point
  #-----------------
  # Make the posterior object as fake error
  posterior <- "Fake error to start"
  class(posterior) <- "try-error"
  # Total reruns until we find a valid initialisation
  total_attempts <- 20
  counter <- 1
  while (inherits(posterior, "try-error") && counter <= total_attempts) {
    posterior <- try(
      .stan_inference(dt = dt, use_mcmc = use_mcmc, use_eb = use_eb,
                      init_using_eb = init_using_eb, init_vals = init_vals,
                      chains = chains, iter = iter, algorithm = algorithm,
                      n_cores = n_cores, tol_rel_obj = tol_rel_obj,
                      output_samples = output_samples, seed = seed, ...)
    )
    if (inherits(posterior, "try-error") ) {
      message("Warning: Inference failed with seed: ", seed,
              "\n Reruning with different seed...\n")
      seed <- sample.int(.Machine$integer.max, 1)
    }
    counter <- counter + 1
  }
  # If the model failed all the attempts, stop and return a message
  if (inherits(posterior, "try-error")) {
    stop("ERROR: Inference failed with these settings.
         Try running again or perfrom some filtering on the data...")
  }
  cat(date(), ": Posterior inference finished.\n")

  ## Round simulated object to reduce memory storage
  tmp <- rstan::extract(posterior)
  posterior <- list(
    mu = round(tmp$mu, 3),
    gamma = round(tmp$gamma, 3),
    epsilon = round(tmp$epsilon, 3),
    w_mu = round(tmp$w_mu, 3),
    w_gamma = round(tmp$w_gamma, 3),
    s_gamma = round(tmp$s_gamma, 3),
    lp__ = tmp$lp__
  )
  colnames(posterior$mu) = colnames(posterior$gamma) =
    colnames(posterior$epsilon) <- unique(Y$Feature)
  colnames(posterior$w_mu) <- colnames(X)
  colnames(posterior$w_gamma) <- paste0("rbf", seq_len(L))
  # Set object class
  cl <- ifelse(use_mcmc, "scmet_mcmc", "scmet_vb")
  # Store list objects
  theta_priors <- list(m_wmu = m_wmu, s_wmu = s_wmu, s_mu = s_mu,
                       m_wgamma = m_wgamma, s_wgamma = s_wgamma,
                       a_sgamma = a_sgamma, b_sgamma = b_sgamma)
  opts <- list(L = L, use_mcmc = use_mcmc, use_eb = use_eb,
               algorithm = algorithm, iter = iter, chains = chains,
               rbf_c = rbf_c, tol_rel_obj = tol_rel_obj,
               output_samples = output_samples, seed = seed)
  obj <- structure(list(posterior = posterior, Y = Y, X = X,
                        feature_names = unique(Y$Feature),
                        theta_priors = theta_priors, opts = opts),
                   class = cl)
  return(obj)
}


.stan_inference <- function(dt, use_mcmc, use_eb, init_using_eb, init_vals,
                            chains, iter, algorithm, n_cores, tol_rel_obj,
                            output_samples, seed, ...) {
  if (use_mcmc) {
    if (iter > 10000) {
      message("Potentially large number of posterior samples for MCMC.\n",
              "Model might take a while to finish. \n")
    }
    # Run MCMC sampling using Stan
    if (!algorithm %in% c("NUTS", "HMC")) { algorithm <- "NUTS" }
    if (use_eb & init_using_eb) {
      posterior <- rstan::sampling(object = stanmodels$scmet, data = dt,
                                   chains = chains, init = init_vals,
                                   seed = seed,
                                   pars = c("mu", "gamma", "w_mu", "w_gamma",
                                            "s_gamma", "epsilon"),
                                   iter = iter, algorithm = algorithm,
                                   cores = n_cores, ...)
    } else {
      posterior <- rstan::sampling(object = stanmodels$scmet, data = dt,
                                   chains = chains, seed = seed,
                                   pars = c("mu", "gamma", "w_mu", "w_gamma",
                                            "s_gamma", "epsilon"),
                                   iter = iter, algorithm = algorithm,
                                   cores = n_cores, ...)
    }
  } else {
    if (!algorithm %in% c("meanfield", "fullrank")) {
      algorithm <- "meanfield"
    }
    if (iter < 2000) {
      message("Small number of VB iterations. scMET might not converge.\n")
    }
    if (use_eb & init_using_eb) {
      posterior <- rstan::vb(object = stanmodels$scmet, data = dt, seed = seed,
                             init = init_vals[[1]], algorithm = algorithm,
                             pars = c("mu", "gamma", "w_mu", "w_gamma",
                                      "s_gamma", "epsilon"),
                             tol_rel_obj = tol_rel_obj,
                             output_samples = output_samples, iter = iter, ...)
    } else {
      posterior <- rstan::vb(object = stanmodels$scmet, data = dt, seed = seed,
                             algorithm = algorithm, tol_rel_obj = tol_rel_obj,
                             pars = c("mu", "gamma", "w_mu", "w_gamma",
                                      "s_gamma", "epsilon"),
                             output_samples = output_samples, iter = iter, ...)
    }
  }
  return(posterior)
}
