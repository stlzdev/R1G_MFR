# imports
library(car)
library(emmeans)

# follow-up for the forward asymmetry interactions (group x presentation rate, group x list length)
d <- read.csv('analyses/dataframes/fwd_asym_data_bsa.csv')
d$subject <- factor(d$subject)
d$strategy <- factor(d$strategy, levels=c('prim', 'ns', 'rec'))
d$l_length <- factor(d$l_length)
d$pres_rate <- factor(d$pres_rate)

# same model as strategy_anova.R
model <- lm(asym ~ strategy*l_length + strategy*pres_rate, data=d)

# (1) group differences at each presentation rate (Tukey within rate)
grp_by_rate <- as.data.frame(summary(pairs(emmeans(model, ~ strategy | pres_rate), adjust='tukey')))
grp_by_rate$cohens_d <- grp_by_rate$estimate / sigma(model)
write.csv(grp_by_rate, 'statistics/dataframes/fwd_asym_group_by_rate.csv', row.names=FALSE)

# (2) effect of presentation rate within each group (Holm across the 3 groups)
rate_by_grp <- as.data.frame(summary(pairs(emmeans(model, ~ pres_rate | strategy)), adjust='holm'))
rate_by_grp$cohens_d <- rate_by_grp$estimate / sigma(model)
write.csv(rate_by_grp, 'statistics/dataframes/fwd_asym_rate_by_group.csv', row.names=FALSE)

# (3) design-clean check: list length 20 is the only length run at both rates (20-1 vs 20-2)
d20 <- d[d$l_length == '20', ]
rows <- list()
for (s in levels(d$strategy)) {
  x <- d20[d20$strategy == s, ]
  tt <- t.test(asym ~ pres_rate, data=x)
  rows[[s]] <- data.frame(strategy=s, n_1000=sum(x$pres_rate == '1000'), n_2000=sum(x$pres_rate == '2000'),
                          mean_1000=mean(x$asym[x$pres_rate == '1000']), mean_2000=mean(x$asym[x$pres_rate == '2000']),
                          t=unname(tt$statistic), df=unname(tt$parameter), p=tt$p.value)
}
rate_ll20 <- do.call(rbind, rows)
rate_ll20$p_holm <- p.adjust(rate_ll20$p, method='holm')
write.csv(rate_ll20, 'statistics/dataframes/fwd_asym_rate_ll20.csv', row.names=FALSE)
