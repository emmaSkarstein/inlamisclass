test_that("modelling, plotting and summary works", {

  # Note: since getting accurate estimates would require running the importance
  # sampling procedure for ~100000 iterations, which takes a long time, the
  # tests in this file only run the IS for 2 iterations to ensure that it
  # actually runs. Therefore, there are not many actual tests in this document,
  # but if all the code runs then that is a great sign.

  library(inlamisclass)
  library(ggplot2)
  library(INLA)
  library(dplyr)

  inla.setOption(num.threads = 1)

  # Birthweight analysis -------------------------------------------------------
  # Misclassification in smoking status assumed
  MC_matrix <- matrix(c(0.95, 0.05, 0.2, 0.8), nrow = 2, byrow = T)

  # Test IS with no imp. covariates
  birthweight_model1 <- inla_is_misclass(formula_moi = bwt ~ smoke + lwt,
                                formula_imp = smoke ~ 1,
                                alpha = c(log(0.4/(1-0.4)), 0),
                                MC_matrix = MC_matrix,
                                data = birthweight, niter = 2, ncores = 2)

  # Test plot_compare_inlamisclass
  naive_mod <- INLA::inla(bwt ~ smoke + lwt, data = birthweight)
  plot_compare_inlamisclass(birthweight_model1, naive_mod = naive_mod)

  # Test IS with one covariate
  birthweight_model2 <- inla_is_misclass(formula_moi = bwt ~ smoke + lwt,
                                formula_imp = smoke ~ lwt,
                                alpha = c(-3, 0.02),
                                MC_matrix = MC_matrix,
                                data = birthweight, niter = 2, ncores = 2)

  # Test plot_compare_inlamisclass with two inla_is_misclass models
  plot_compare_inlamisclass(is_results = list(inlamisclass1 = birthweight_model1,
                                                inlamisclass2 = birthweight_model2),
                            naive_mod = naive_mod, niter = 2, num_inlamisclass_models = 2)



  # Case-control analysis ------------------------------------------------------
  validation <- filter(cervical_cancer, !is.na(x))
  incomplete_data <- filter(cervical_cancer, is.na(x))

  # Estimates misclassification matrix without accounting for y
  # w = 1 given x = 0
  pi10 <- sum((validation$w-validation$x) == 1)/sum(validation$x==0)
  # w = 0, given x = 1, y = 1
  pi01 <- sum(validation$w-validation$x == -1)/sum(validation$x==1)

  M <- matrix(c(1-pi10, pi10, pi01, 1-pi01), byrow = TRUE, nrow = 2)

  exp_mod <- inla(x~1, data = validation, family = "binomial", Ntrials = 1)

  case_control_model <- inla_is_misclass(formula_moi = y ~ w,
                                formula_imp = w ~ 1,
                                alpha = exp_mod$summary.fixed[1,1],
                                MC_matrix = M,
                                data = incomplete_data,
                                niter = 2,
                                family = "binomial", Ntrials = 1, ncores = 2)

  naive_cc <- inla(y ~ w, data = incomplete_data, family = "binomial", Ntrial = 1)

  plot_compare_inlamisclass(case_control_model, naive_mod = naive_cc, plot_intercept = TRUE)

  # Conditional model to deal with differential MC in case-control data
  validation <- filter(cervical_cancer, !is.na(x))
  incomplete_data <- filter(cervical_cancer, is.na(x))

  validation1 <- filter(validation, y == 1)
  validation0 <- filter(validation, y == 0)

  # w = 1 given x = 0, y = 1
  pi10_y1 <- sum((validation1$w-validation1$x) == 1)/sum(validation1$x==0)
  # w = 0, given x = 1, y = 1
  pi01_y1 <- sum(validation1$w-validation1$x == -1)/sum(validation1$x==1)

  # w = 1, x = 0, given y = 0
  pi10_y0 <- sum((validation0$w-validation0$x) == 1)/sum(validation0$x==0)
  # w = 0, x = 1, given y = 0
  pi01_y0 <- sum(validation0$w-validation0$x == -1)/sum(validation0$x==1)

  M1 <- matrix(c(1-pi10_y1, pi10_y1, pi01_y1, 1-pi01_y1), byrow = TRUE, nrow = 2)
  M0 <- matrix(c(1-pi10_y0, pi10_y0, pi01_y0, 1-pi01_y0), byrow = TRUE, nrow = 2)

  exp_glm <- glm(x~y, family = "binomial", data = validation)$coef
  alphas <- exp_glm["(Intercept)"] + c(0, exp_glm["y"])

  # Test sampling x conditional on y
  pi_cond <- new_pi_conditional(alpha = alphas, z = cervical_cancer$y,
                                MC_matrix = list(MC_0 = M0, MC_1 = M1),
                                w = cervical_cancer$w,
                                conditioning_var = cervical_cancer$y)

  case_control_model <- inla_is_misclass(formula_moi = y ~ w,
                                formula_imp = w ~ y,
                                alpha = alphas,
                                MC_matrix = list(MC_1 = M1, MC_0 = M0),
                                data = incomplete_data,
                                niter = 2,
                                conditional = "y",
                                family = "binomial", Ntrials = 1, ncores = 2)

  case_control_model3 <- inla_is_misclass(formula_moi = y ~ w,
                                 formula_imp = w ~ 1,
                                 alpha = exp_glm["(Intercept)"],
                                 MC_matrix = list(MC_1 = M1, MC_0 = M0),
                                 data = incomplete_data,
                                 niter = 2,
                                 conditional = "y",
                                 family = "binomial", Ntrials = 1, ncores = 2)


  # Simulations ----------------------------------------------------------------
  MC_matrix <- matrix(c(0.9, 0.1, 0.2, 0.8), nrow = 2, byrow = T)

  set.seed(1)
  data2 <- generate_misclassified(n = 200, p = 2, MC_matrix = MC_matrix,
                                  betas = c(1, 1, 1),
                                  alphas = c(-0.5, 0.25))
  # Test sampling x
  pi <- new_pi(alpha = c(-0.5, 0.25), z = data2$z, MC_matrix = MC_matrix, w = data2$w)

  model2 <- inla_is_misclass(formula_moi = y ~ w + z,
                    formula_imp = w ~ z,
                    alpha = c(-0.5, 0.25),
                    MC_matrix = MC_matrix,
                    data = data2,
                    niter = 2, ncores = 2)

  naive2 <- inla(y ~ w + z, data = data2)
  correct2 <- inla(y ~ x + z, data = data2)

  plot_compare_inlamisclass(is_results = model2, naive_mod = naive2,
                            correct_mod = correct2,
                            plot_intercept = TRUE, niter = 2)


  # Missing observations -------------------------------------------------------

  # Assuming misclassification and missingness
  MC_matrix <- matrix(c(0.9, 0.1, 0.2, 0.8), nrow = 2, byrow = T)

  set.seed(1)
  data_miss <- generate_misclassified(n = 200, p = 2, MC_matrix = MC_matrix,
                                  betas = c(1, 1, 1),
                                  alphas = c(-0.5, 0.25))
  data_miss$w[1:10] <- NA

  model_miss <- inla_is_misclass(formula_moi = y ~ w + z,
                             formula_imp = w ~ z,
                             alpha = c(-0.5, 0.25),
                             MC_matrix = MC_matrix,
                             data = data_miss,
                             niter = 2, ncores = 2)

  # Assuming only missingness, no MC
  model_miss2 <- inla_is_misclass(formula_moi = y ~ w + z,
                                 formula_imp = w ~ z,
                                 alpha = c(-0.5, 0.25),
                                 data = data_miss,
                                 missing_only = TRUE,
                                 niter = 2, ncores = 2)

  # Check if no MC matrix and missing_only = FALSE
  expect_error(inla_is_misclass(formula_moi = y ~ w + z,
                                formula_imp = w ~ z,
                                alpha = c(-0.5, 0.25),
                                data = data_miss,
                                niter = 2, ncores = 2))

  # Check if MC_matrix and missing_only = TRUE
  expect_error(inla_is_misclass(formula_moi = y ~ w + z,
                                formula_imp = w ~ z,
                                alpha = c(-0.5, 0.25),
                                MC_matrix = MC_matrix,
                                data = data_miss,
                                missing_only = TRUE,
                                niter = 2, ncores = 2))

})
