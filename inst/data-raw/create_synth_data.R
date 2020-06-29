library(scMET)
set.seed(123)
###-----------------------------------
# Generate synthetic data
###-----------------------------------
scmet_dt <- scmet_simulate(N_feat = 150, N_cells = 50, N_cpgs = 15, L = 4, X = NULL,
                           w_mu = c(-0.5, -1.5), s_mu = 1,
                           w_gamma = c(-1.2, -.3, 1.1, -.9), s_gamma = 0.2, rbf_c = 1,
                           fit_linear_trend = FALSE, cells_range = c(0.4, 0.8),
                           cpgs_range = c(0.4, 0.8), seed = 123)
usethis::use_data(scmet_dt, overwrite = TRUE)


###-----------------------------------
# Generate synthetic differential mean data
###-----------------------------------
scmet_diff_dt <- scmet_simulate_diff(N_feat = 200, N_cells = 100, N_cpgs = 15, L = 4,
                                     diff_feat_prcg_mu = 0.05, diff_feat_prcg_gamma = 0.2,
                                     OR_change_mu = 3, OR_change_gamma = 3, X = NULL,
                                     w_mu = c(-0.5, -1.5), s_mu = 1,
                                     w_gamma = c(-1.2, -.3, 1.1, -.9), s_gamma = 0.2,
                                     rbf_c = 1, fit_linear_trend = FALSE,
                                     cells_range = c(0.4, 0.8), cpgs_range = c(0.4, 0.8),
                                     seed = 123)
usethis::use_data(scmet_diff_dt, overwrite = TRUE)


# ###-----------------------------------
# # Generate synthetic differential variability data
# ###-----------------------------------
# scmet_diffvar_dt <- scmet_simulate_diff(N_feat = 300, N_cells = 50, N_cpgs = 15, L = 4,
#                                         diff_feat_prcg_mu = 0, diff_feat_prcg_gamma = 0.2,
#                                         OR_change_mu = 5, OR_change_gamma = 5, X = NULL,
#                                         w_mu = c(-0.5, -1.5), s_mu = 1,
#                                         w_gamma = c(-1.2, -.3, 1.1, -.9), s_gamma = 0.2,
#                                         rbf_c = 1,  fit_linear_trend = FALSE,
#                                         cells_range = c(0.4, 0.8), cpgs_range = c(0.4, 0.8),
#                                         seed = 123)
# usethis::use_data(scmet_diffvar_dt, overwrite = TRUE)
