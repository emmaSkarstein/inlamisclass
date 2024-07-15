## code to prepare `birthweight` dataset goes here

library(aplore3)

dd <- lowbwt

# Recode smoking as 1 (yes) and 0 (no)
dd$smoke <- ifelse(dd$smoke == "Yes", 1, 0)

# Rescale birth weight to kg
dd$bwt <- dd$bwt/1000

birthweight <- data.frame(bwt = scale(dd$bwt, scale=F),
                 lwt = scale(dd$lwt, scale=F),
                 smoke = dd$smoke)



usethis::use_data(birthweight, overwrite = TRUE)
