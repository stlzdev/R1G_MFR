# imports
library(lme4)
library(car)

add_partial_eta_squared <- function(res_anova) {
  residual_ss <- res_anova["Residuals", "Sum Sq"]
  effect_rows <- !(rownames(res_anova) %in% c("(Intercept)", "Residuals"))
  
  res_anova$eta_p_sq <- NA_real_
  res_anova[effect_rows, "eta_p_sq"] <- res_anova[effect_rows, "Sum Sq"] /
    (res_anova[effect_rows, "Sum Sq"] + residual_ss)
  
  return(res_anova)
}

# recall initiation bias
r1_data <- read.csv('analyses/dataframes/prim_rec_pfr.csv')

# treat categorical variables as factors
r1_data$subject <- factor(r1_data$subject)
r1_data$l_length <- factor(r1_data$l_length)
r1_data$pres_rate <- factor(r1_data$pres_rate)

# 2-factor ANOVA
model <- lm(rec_prim_bias ~ l_length + pres_rate, data=r1_data)
r1_anova <- Anova(model, type='III')

# calculate mean squared errors
r1_anova$mse <- r1_anova$`Sum Sq` / r1_anova$Df
r1_anova <- add_partial_eta_squared(r1_anova)

# save out results
write.csv(r1_anova, 'statistics/dataframes/r1_anova.csv')
