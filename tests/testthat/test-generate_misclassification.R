test_that("generate_misclassified works", {
  M <- matrix(c(0.9, 0.1, 0.2, 0.8), nrow = 2, byrow = T)

  set.seed(1)
  data <- generate_misclassified(n = 5000, p = 2, MC_matrix = M,
                                  betas = c(1, 1, 1),
                                  alphas = c(-0.5, 0.25))


  # pi10
  sum((data$w-data$x) == 1)/sum(data$x==0)
  M[1,2]
  expect_equal(sum((data$w-data$x) == 1)/sum(data$x==0), M[1,2], tolerance = 0.1)

  # pi01
  sum((data$w-data$x) == -1)/sum(data$x==1)
  M[2,1]
  expect_equal(sum((data$w-data$x) == -1)/sum(data$x==1), M[2,1], tolerance = 0.1)

  # Total prob. of misclass
  sum(abs(data$w - data$x)) / nrow(data)
  M[1,2]*sum(data$x==0)/nrow(data) + M[2,1]*sum(data$x==1)/nrow(data)

  # Check correct model
  correct_coef <- summary(lm(y~x + z, data = data))$coef
  correct_coef

  # Attenuated version
  naive_coef <- summary(lm(y~w + z, data = data))$coef
  naive_coef
})
