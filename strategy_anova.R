# imports
library(lme4)
library(car)
library(emmeans)
library(glue)

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

# three-factor ANOVA
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


# mean words recalled
mwr_data_bsa <- read.csv('analyses/dataframes/mwr_data_bsa.csv')
mwr_res <- run_anova(mwr_data_bsa, 'mwr')
mwr_anova <- mwr_res$res_anova
mwr_model <- mwr_res$model
write.csv(mwr_anova, 'statistics/dataframes/mwr_anova.csv')

mwr_tukey_s <- post_hoc_pairwise(mwr_anova, mwr_model)
if (!is.null(mwr_tukey_s)) {
  write.csv(mwr_tukey_s, 'statistics/dataframes/mwr_tukey_s.csv', row.names=FALSE)
}


# R1 intrusions
r1_intr_data_bsa <- read.csv('analyses/dataframes/r1_intr_data_bsa.csv')
r1_intr_res <- run_anova(r1_intr_data_bsa, 'prop_wrong')
r1_intr_anova <- r1_intr_res$res_anova
r1_intr_model <- r1_intr_res$model
write.csv(r1_intr_anova, 'statistics/dataframes/r1_intr_anova.csv')

r1_intr_tukey_s <- post_hoc_pairwise(r1_intr_anova, r1_intr_model)
if (!is.null(r1_intr_tukey_s)) {
  write.csv(r1_intr_tukey_s, 'statistics/dataframes/r1_intr_tukey_s.csv', row.names=FALSE)
}


# initial response times
rti_data_only_cr_bsa <- read.csv('analyses/dataframes/rti_data_only_cr_bsa.csv')
rti_res <- run_anova(rti_data_only_cr_bsa, 'rt')
rti_anova <- rti_res$res_anova
rti_model <- rti_res$model
write.csv(rti_anova, 'statistics/dataframes/rti_anova.csv')

rti_tukey_s <- post_hoc_pairwise(rti_anova, rti_model)
if (!is.null(rti_tukey_s)) {
  write.csv(rti_tukey_s, 'statistics/dataframes/rti_tukey_s.csv', row.names=FALSE)
}


# PLI and ELI rates
intr_data_only_cr_bsa <- read.csv('analyses/dataframes/intr_data_only_cr_bsa.csv')

pli_res <- run_anova(intr_data_only_cr_bsa, 'pli_rate')
pli_anova <- pli_res$res_anova
pli_model <- pli_res$model
write.csv(pli_anova, 'statistics/dataframes/pli_anova.csv')

pli_tukey_s <- post_hoc_pairwise(pli_anova, pli_model)
if (!is.null(pli_tukey_s)) {
  write.csv(pli_tukey_s, 'statistics/dataframes/pli_tukey_s.csv', row.names=FALSE)
}

eli_res <- run_anova(intr_data_only_cr_bsa, 'eli_rate')
eli_anova <- eli_res$res_anova
eli_model <- eli_res$model
write.csv(eli_anova, 'statistics/dataframes/eli_anova.csv')

eli_tukey_s <- post_hoc_pairwise(eli_anova, eli_model)
if (!is.null(eli_tukey_s)) {
  write.csv(eli_tukey_s, 'statistics/dataframes/eli_tukey_s.csv', row.names=FALSE)
}


# temporal clustering score (w/ H2-H1)
tcl_data_bsa <- read.csv('analyses/dataframes/tcl_data_bsa.csv')
tcl_res <- run_anova(tcl_data_bsa, 'tcl')
tcl_anova <- tcl_res$res_anova
tcl_model <- tcl_res$model
write.csv(tcl_anova, 'statistics/dataframes/tcl_anova.csv')

tcl_tukey_s <- post_hoc_pairwise(tcl_anova, tcl_model)
if (!is.null(tcl_tukey_s)) {
  write.csv(tcl_tukey_s, 'statistics/dataframes/tcl_tukey_s.csv', row.names=FALSE)
}

tcl_h_data_bsa <- read.csv('analyses/dataframes/tcl_h_data_bsa.csv')
tcl_h_res <- run_anova(tcl_h_data_bsa, 'tcl_delta')
tcl_h_anova <- tcl_h_res$res_anova
tcl_h_model <- tcl_h_res$model
write.csv(tcl_h_anova, 'statistics/dataframes/tcl_h_anova.csv')

tcl_h_tukey_s <- post_hoc_pairwise(tcl_h_anova, tcl_h_model)
if (!is.null(tcl_h_tukey_s)) {
  write.csv(tcl_h_tukey_s, 'statistics/dataframes/tcl_h_tukey_s.csv', row.names=FALSE)
}


# forward asymmetry (+1 lag - -1 lag)
fwd_asym_data_bsa <- read.csv('analyses/dataframes/fwd_asym_data_bsa.csv')
fwd_asym_res <- run_anova(fwd_asym_data_bsa, 'asym')
fwd_asym_anova <- fwd_asym_res$res_anova
fwd_asym_model <- fwd_asym_res$model
write.csv(fwd_asym_anova, 'statistics/dataframes/fwd_asym_anova.csv')

fwd_asym_tukey_s <- post_hoc_pairwise(fwd_asym_anova, fwd_asym_model)
if (!is.null(fwd_asym_tukey_s)) {
  write.csv(fwd_asym_tukey_s, 'statistics/dataframes/fwd_asym_tukey_s.csv', row.names=FALSE)
}


# semantic clustering score
scl_data_bsa <- read.csv('analyses/dataframes/scl_data_bsa.csv')
scl_res <- run_anova(scl_data_bsa, 'scl')
scl_anova <- scl_res$res_anova
scl_model <- scl_res$model
write.csv(scl_anova, 'statistics/dataframes/scl_anova.csv')

scl_tukey_s <- post_hoc_pairwise(scl_anova, scl_model)
if (!is.null(scl_tukey_s)) {
  write.csv(scl_tukey_s, 'statistics/dataframes/scl_tukey_s.csv', row.names=FALSE)
}
