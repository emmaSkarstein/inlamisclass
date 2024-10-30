# Create data set described in Carroll, Gail and Lubin (1993) (https://www.jstor.org/stable/pdf/2290713.pdf)
library(dplyr)

case_control_data_raw1 <- matrix(c(
  rep(c(1,0,0), 13),
  rep(c(1,0,1), 3),
  rep(c(1,1,0), 5),
  rep(c(1,1,1), 18),
  rep(c(0,0,0), 33),
  rep(c(0,0,1), 11),
  rep(c(0,1,0), 16),
  rep(c(0,1,1), 16),
  rep(c(1,NA,0), 318),
  rep(c(1,NA,1), 375),
  rep(c(0,NA,0), 701),
  rep(c(0,NA,1), 535)),
  ncol = 3, byrow = TRUE)


case_control_data_raw2 <- data.frame(case_control_data_raw1)

colnames(case_control_data_raw2) <- c("y", "x", "w")

cervical_cancer <- case_control_data_raw2 #%>%
  #mutate(complete = ifelse(is.na(x), 0, 1))

usethis::use_data(cervical_cancer, overwrite = TRUE)

