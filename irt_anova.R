# imports
library(lme4)
library(car)
library(emmeans)

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

# 4-factor ANOVA
# inter-response times for final 4 transitions
run_irt_anova <- function(irt_data) {
  # treat categorical variables as factors
  irt_data$subject <- factor(irt_data$subject)
  irt_data$strategy <- factor(irt_data$strategy)
  irt_data$rot <- factor(irt_data$rot)
  irt_data$l_length <- factor(irt_data$l_length)
  irt_data$pres_rate <- factor(irt_data$pres_rate)
  
  model <- lm(irt ~ strategy + rot + l_length + pres_rate + strategy:rot + strategy:l_length + strategy:pres_rate + rot:l_length + rot:pres_rate, data=irt_data)
  res_anova <- Anova(model, type='III')
  
  # calculate mean squared errors
  res_anova$mse <- res_anova$`Sum Sq` / res_anova$Df
  res_anova <- add_partial_eta_squared(res_anova)
  
  return (list(res_anova=res_anova, model=model))
}

# post-hoc pairwise tests
post_hoc_pairwise <- function(res_anova, model) {
  me_pval <- res_anova["strategy", "Pr(>F)"]
  if (TRUE) {  # always run so output files never go stale; interpret alongside the omnibus p (me_pval)
    emm <- emmeans(model, ~ strategy)
    res_tukey <- pairs(emm, adjust='tukey')
    res_tukey_s <- as.data.frame(summary(res_tukey))
    res_tukey_s$cohens_d <- res_tukey_s$estimate / sigma(model)
    return (res_tukey_s)
  } else {
    print('No significant main effect of strategy')
    return (NULL)
  }
}

irt_data <- read.csv('analyses/dataframes/irt_final_4_data_bsa.csv')
irt_res <- run_irt_anova(irt_data)
irt_anova <- irt_res$res_anova
irt_model <- irt_res$model
write.csv(irt_anova, 'statistics/dataframes/irt_anova.csv')


irt_tukey_s <- post_hoc_pairwise(irt_anova, irt_model)
if (!is.null(irt_tukey_s)) {
  write.csv(irt_tukey_s, 'statistics/dataframes/irt_tukey_s.csv', row.names=FALSE)
}
