library(readxl)

icb_cohort <- read_excel(
  "IFNG_ICB/1-s2.0-S1535610823000417-mmc2.xlsx",
  sheet = "ICB cohort"
)

studies <- unique(icb_cohort[, 1])