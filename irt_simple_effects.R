# imports
library(lme4)
library(car)
library(emmeans)

# subject-level (repeated-measures) follow-up for the final-4 IRT analysis
# subject is a random intercept, so the 4 transitions per subject are not treated as independent
# between subjects: strategy, l_length, pres_rate; within subjects: rot
irt_data <- read.csv('analyses/dataframes/irt_final_4_data_bsa.csv')
irt_data$subject <- factor(irt_data$subject)
irt_data$strategy <- factor(irt_data$strategy)
irt_data$rot <- factor(irt_data$rot)
irt_data$l_length <- factor(irt_data$l_length)
irt_data$pres_rate <- factor(irt_data$pres_rate)

# sanity check: each subject has one strategy, one condition, and one row per transition
stopifnot(all(tapply(irt_data$strategy, irt_data$subject, function(x) length(unique(x))) == 1))
stopifnot(all(table(irt_data$subject) == nlevels(irt_data$rot)))

# same fixed effects as irt_anova.R, plus a random intercept for subject
options(contrasts = c('contr.sum', 'contr.poly'))
model <- lmer(irt ~ strategy + rot + l_length + pres_rate + strategy:rot + strategy:l_length +
                strategy:pres_rate + rot:l_length + rot:pres_rate + (1 | subject), data=irt_data)

# type III F tests with Kenward-Roger denominator df
mixed_anova <- as.data.frame(Anova(model, type='III', test.statistic='F'))
write.csv(mixed_anova, 'statistics/dataframes/irt_mixed_anova.csv')

# simple effects: recall initiation group at each relative output transition
# Tukey-adjusted within each transition
simple_rot <- as.data.frame(summary(pairs(emmeans(model, ~ strategy | rot), adjust='tukey')))
write.csv(simple_rot, 'statistics/dataframes/irt_simple_effects_rot.csv', row.names=FALSE)

# simple effects: group at each list length, and at each presentation rate
simple_ll <- as.data.frame(summary(pairs(emmeans(model, ~ strategy | l_length), adjust='tukey')))
write.csv(simple_ll, 'statistics/dataframes/irt_simple_effects_l_length.csv', row.names=FALSE)

simple_pr <- as.data.frame(summary(pairs(emmeans(model, ~ strategy | pres_rate), adjust='tukey')))
write.csv(simple_pr, 'statistics/dataframes/irt_simple_effects_pres_rate.csv', row.names=FALSE)

# pairwise group comparisons, averaged over transitions, list lengths, and presentation rates
# Cohen's d standardized by total SD (residual + subject variance), comparable to the lm-based d
vc <- as.data.frame(VarCorr(model))
total_sd <- sqrt(sum(vc$vcov))
main_pairs <- as.data.frame(summary(pairs(emmeans(model, ~ strategy), adjust='tukey')))
main_pairs$cohens_d <- main_pairs$estimate / total_sd
write.csv(main_pairs, 'statistics/dataframes/irt_mixed_tukey.csv', row.names=FALSE)
