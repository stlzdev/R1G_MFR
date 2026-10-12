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

# mean words recalled in practice session
# analyzed participants v. participants who did not complete at least three sessions
mwr_practice_data <- read.csv('analyses/dataframes/mwr_practice_data.csv')
mwr_practice_data <- subset(mwr_practice_data, status != 'other')

# treat categorical variables as factors
mwr_practice_data$subject <- factor(mwr_practice_data$subject)
mwr_practice_data$status <- factor(mwr_practice_data$status)
mwr_practice_data$condition <- factor(mwr_practice_data$condition)

# 2-factor ANOVA (practice session condition)
model <- lm(mwr ~ status*condition, data=mwr_practice_data)
attrition_anova <- Anova(model, type='III')

# calculate mean squared errors
attrition_anova$mse <- attrition_anova$`Sum Sq` / attrition_anova$Df
attrition_anova <- add_partial_eta_squared(attrition_anova)

# save out results
write.csv(attrition_anova, 'statistics/dataframes/attrition_anova.csv')

# estimated marginal means and pairwise comparison
attrition_emm <- as.data.frame(summary(emmeans(model, ~ status)))
write.csv(attrition_emm, 'statistics/dataframes/attrition_emm.csv', row.names=FALSE)

attrition_pairs <- as.data.frame(summary(pairs(emmeans(model, ~ status))))
attrition_pairs$cohens_d <- attrition_pairs$estimate / sigma(model)
write.csv(attrition_pairs, 'statistics/dataframes/attrition_pairs.csv', row.names=FALSE)
