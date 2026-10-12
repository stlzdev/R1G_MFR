# imports
library(lme4)
library(car)

# Sum-to-zero contrasts: required for meaningful Type III tests of main effects
# (with R's default treatment coding, Type III 'main effects' are simple effects at the reference cell)
options(contrasts = c('contr.sum', 'contr.poly'))


add_partial_eta_squared <- function(res_anova) {
  residual_ss <- res_anova["Residuals", "Sum Sq"]
  effect_rows <- !(rownames(res_anova) %in% c("(Intercept)", "Residuals"))

  res_anova$eta_p_sq <- NA_real_
  res_anova[effect_rows, "eta_p_sq"] <- res_anova[effect_rows, "Sum Sq"] /
    (res_anova[effect_rows, "Sum Sq"] + residual_ss)

  return(res_anova)
}

# 1-factor ANOVA of list length within one presentation rate
# (presentation rate varies only for 20-item lists, so list length is tested within each rate)
run_ll_anova <- function(mwr_data) {
  # treat categorical variables as factors
  mwr_data$subject <- factor(mwr_data$subject)
  mwr_data$l_length <- factor(mwr_data$l_length)

  model <- lm(mwr ~ l_length, data=mwr_data)
  res_anova <- Anova(model, type='III')

  # calculate mean squared errors
  res_anova$mse <- res_anova$`Sum Sq` / res_anova$Df
  res_anova <- add_partial_eta_squared(res_anova)

  return (res_anova)
}

# mean words recalled, before partitioning by recall initiation group
mwr_data_bsa <- read.csv('analyses/dataframes/mwr_data_bsa.csv')

# 2 s per item: 10, 15, 20 items
mwr_ll_2s_anova <- run_ll_anova(subset(mwr_data_bsa, pres_rate == 2000))
write.csv(mwr_ll_2s_anova, 'statistics/dataframes/mwr_ll_2s_anova.csv')

# 1 s per item: 20, 30, 40 items
mwr_ll_1s_anova <- run_ll_anova(subset(mwr_data_bsa, pres_rate == 1000))
write.csv(mwr_ll_1s_anova, 'statistics/dataframes/mwr_ll_1s_anova.csv')
