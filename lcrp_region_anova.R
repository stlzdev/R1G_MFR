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

# three-factor ANOVA (same model as strategy_anova.R)
run_anova <- function(data_bsa, dv) {
  # treat categorical variables as factors
  data_bsa$subject <- factor(data_bsa$subject)
  data_bsa$strategy <- factor(data_bsa$strategy)
  data_bsa$l_length <- factor(data_bsa$l_length)
  data_bsa$pres_rate <- factor(data_bsa$pres_rate)

  formula_str <- paste(dv, '~ strategy*l_length + strategy*pres_rate')
  formula <- as.formula(formula_str)
  model <- lm(formula, data=data_bsa)
  res_anova <- Anova(model, type='III')

  # calculate mean squared errors
  res_anova$mse <- res_anova$`Sum Sq` / res_anova$Df
  res_anova <- add_partial_eta_squared(res_anova)

  return (list(res_anova=res_anova, model=model))
}

# post-hoc pairwise tests
post_hoc_pairwise <- function(model) {
  emm <- emmeans(model, ~ strategy)
  res_tukey_s <- as.data.frame(summary(pairs(emm, adjust='tukey')))
  res_tukey_s$cohens_d <- res_tukey_s$estimate / sigma(model)
  return (res_tukey_s)
}


# forward asymmetry (lag +1 - lag -1) for transitions from middle-list items
# (controls for serial position of transitioning item, Reviewer 3)
lcrp_region_data_bsa <- read.csv('analyses/dataframes/lcrp_region_data_bsa.csv')
fwd_asym_mid_res <- run_anova(subset(lcrp_region_data_bsa, region == 'middle'), 'asym')
fwd_asym_mid_anova <- fwd_asym_mid_res$res_anova
write.csv(fwd_asym_mid_anova, 'statistics/dataframes/fwd_asym_mid_anova.csv')

fwd_asym_mid_tukey_s <- post_hoc_pairwise(fwd_asym_mid_res$model)
write.csv(fwd_asym_mid_tukey_s, 'statistics/dataframes/fwd_asym_mid_tukey_s.csv', row.names=FALSE)
