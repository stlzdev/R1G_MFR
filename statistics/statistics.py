### Statistics code
### FDR correction with Benjamini-Yekutieli method

# imports
import pandas as pd
import numpy as np
import scipy.stats


# ---------- Utility ----------

# sort dataframe by condition
# list length and presentation rate or total time
def sort_by_condition(df, style='ll_pr'):
    if style=='total_time':
        conds = ['10-2', '20-1', '15-2', '30-1', '20-2', '40-1']
    elif style=='ll_pr':
        conds = ['10-2', '15-2', '20-2', '20-1', '30-1', '40-1']
    else:
        raise ValueError(f"Invalid category sort style: {style}.")
        
    df['condition'] = pd.Categorical(df['condition'], categories = conds)
    df = df.sort_values(by='condition', ignore_index=True)
    
    return df


# sort by condition and strategy
def sort_by_condition_strategy(df):
    strats = ['prim', 'ns', 'rec']
    conds_ll_pr = ['10-2', '15-2', '20-2', '20-1', '30-1', '40-1']
    df['condition'] = pd.Categorical(df['condition'], categories = conds_ll_pr)
    df['strategy'] = pd.Categorical(df['strategy'], categories = strats)
    df = df.sort_values(by=['condition', 'strategy'], ignore_index=True)
    
    return df


def cohens_d_1samp(sample, popmean=0):
    sample = pd.Series(sample).dropna()
    sample_std = sample.std(ddof=1)

    if len(sample) < 2 or sample_std == 0:
        return np.nan

    return (sample.mean() - popmean) / sample_std


# ---------- Hypothesis Testing ----------

# primacy and recency effect
# 1-sample t-test against mean of 0
def prim_rec_slopes_statistics(spc_prim_rec_lr_all):
    stats = []
    for (c, ll, pr), cond_data in spc_prim_rec_lr_all.groupby(['condition', 'l_length', 'pres_rate']):
        prim_res = scipy.stats.ttest_1samp(cond_data.prim_slope, popmean=0, alternative='two-sided')
        rec_res = scipy.stats.ttest_1samp(cond_data.rec_slope, popmean=0, alternative='two-sided')
        prim_d = cohens_d_1samp(cond_data.prim_slope, popmean=0)
        rec_d = cohens_d_1samp(cond_data.rec_slope, popmean=0)
        stats.append((c, ll, pr, prim_res.statistic, prim_res.pvalue, prim_res.df, prim_d,
                      rec_res.statistic, rec_res.pvalue, rec_res.df, rec_d))
        
    stats = pd.DataFrame(stats, columns=['condition', 'l_length', 'pres_rate', 
                                         'prim_t_stat', 'prim_p_val', 'prim_dof', 'prim_cohens_d',
                                         'rec_t_stat', 'rec_p_val', 'rec_dof', 'rec_cohens_d'])
    
    # FDR correction
    all_pvals = list(stats.prim_p_val) + list(stats.rec_p_val)
    fdr_pvals = scipy.stats.false_discovery_control(all_pvals, method='by')
    fdr_prim, fdr_rec = np.split(fdr_pvals, 2)
    
    stats['prim_p_val_fdr'] = fdr_prim
    stats['rec_p_val_fdr'] = fdr_rec
    
    # sort by condition
    stats = sort_by_condition(stats, style='ll_pr')
    
    return stats


# primacy and recency initiation bias
# 1-sample t-test against mean of 0
def prim_rec_pfr_statistics(prim_rec_pfr):
    stats = []
    for (c, ll, pr), cond_data in prim_rec_pfr.groupby(['condition', 'l_length', 'pres_rate']):
        res = scipy.stats.ttest_1samp(cond_data.rec_prim_bias, popmean=0, alternative='two-sided')
        stats.append((c, ll, pr, res.statistic, res.pvalue, res.df))
        
    stats = pd.DataFrame(stats, columns=['condition', 'l_length', 'pres_rate', 't_stat', 'p_val', 'dof'])
    
    # FDR correction
    fdr_pvals = scipy.stats.false_discovery_control(stats.p_val, method='by')
    stats['p_val_fdr'] = fdr_pvals
    
    # sort by condition
    stats = sort_by_condition(stats, style='ll_pr')
    
    return stats


# correlation of PFR and SPC (at each serial position)
def pfr_spc_correlation_statistics(pfr_spc_corrs):
    stats = pfr_spc_corrs.copy()
    
    # FDR correction (all comparisons)
    fdr_pvals = scipy.stats.false_discovery_control(pfr_spc_corrs.p_value, method='by')
    stats['p_val_fdr'] = fdr_pvals
    
    # sort by condition
    stats = sort_by_condition(stats, style='total_time')
    
    return stats


# within session recall initiation variance
# paired t-test
def r1_variance_statistics(r1_var_data_bsa):
    stats = []
    for (c, ll, pr), cond_data in r1_var_data_bsa.groupby(['condition', 'l_length', 'pres_rate']):
        res = scipy.stats.ttest_rel(cond_data.sp_sem, cond_data.permutation_sp_sem, nan_policy='omit', alternative='two-sided')
        diff = cond_data.sp_sem - cond_data.permutation_sp_sem
        cohens_dz = cohens_d_1samp(diff, popmean=0)
        stats.append((c, ll, pr, res.statistic, res.pvalue, res.df, cohens_dz))
        
    stats = pd.DataFrame(stats, columns=['condition', 'l_length', 'pres_rate', 't_stat', 'p_val', 'dof', 'cohens_dz'])
    
    # FDR correction
    fdr_pvals = scipy.stats.false_discovery_control(stats.p_val, method='by')
    stats['p_val_fdr'] = fdr_pvals
    
    # sort by condition
    stats = sort_by_condition(stats, style='ll_pr')
    
    return stats


# change in recall initiation serial position across sessions
# 1-sample t-test against mean of 0
def r1_sp_statistics(r1_sp_dec_data_bsa):
    stats = []
    for (c, ll, pr), cond_data in r1_sp_dec_data_bsa.groupby(['condition', 'l_length', 'pres_rate']):
        res = scipy.stats.ttest_1samp(cond_data.r1_sp_slope, popmean=0, alternative='two-sided')
        cohens_d = cohens_d_1samp(cond_data.r1_sp_slope, popmean=0)
        stats.append((c, ll, pr, res.statistic, res.pvalue, res.df, cohens_d))
        
    stats = pd.DataFrame(stats, columns=['condition', 'l_length', 'pres_rate', 't_stat', 'p_val', 'dof', 'cohens_d'])
    
    # FDR corrections
    fdr_pvals = scipy.stats.false_discovery_control(stats.p_val, method='by')
    stats['p_val_fdr'] = fdr_pvals
    
    # sort by condition
    stats = sort_by_condition(stats, style='ll_pr')
    
    return stats

# effect of presentation rate on mean words recalled (before partitioning by recall initiation group)
# independent samples t-test, 20-2 v. 20-1 (presentation rate varies only for 20-item lists)
def mwr_pr_statistics(mwr_data_bsa):
    slow = mwr_data_bsa[mwr_data_bsa.condition == '20-2'].mwr
    fast = mwr_data_bsa[mwr_data_bsa.condition == '20-1'].mwr
    res = scipy.stats.ttest_ind(slow, fast, alternative='two-sided')

    # Cohen's d with pooled standard deviation
    pooled_sd = np.sqrt(((len(slow) - 1) * slow.var(ddof=1) + (len(fast) - 1) * fast.var(ddof=1)) / (len(slow) + len(fast) - 2))
    cohens_d = (np.mean(slow) - np.mean(fast)) / pooled_sd

    stats = pd.DataFrame([(np.mean(slow), np.mean(fast), res.statistic, res.pvalue, res.df, cohens_d)],
                         columns=['mwr_20_2', 'mwr_20_1', 't_stat', 'p_val', 'dof', 'cohens_d'])

    return stats

# change in recall initiation from early to late lists within session
# paired t-test (all conditions together)
def r1_list_statistics(r1_list_data_bsa):
    data = r1_list_data_bsa.dropna(subset=['prim_early', 'prim_late', 'rec_early', 'rec_late'])
    stats = []
    for r1 in ['prim', 'rec']:
        res = scipy.stats.ttest_rel(data[f'{r1}_late'], data[f'{r1}_early'], alternative='two-sided')
        cohens_dz = cohens_d_1samp(data[f'{r1}_late'] - data[f'{r1}_early'], popmean=0)
        stats.append((r1, np.mean(data[f'{r1}_early']), np.mean(data[f'{r1}_late']), res.statistic, res.pvalue, res.df, cohens_dz))

    stats = pd.DataFrame(stats, columns=['initiation', 'early', 'late', 't_stat', 'p_val', 'dof', 'cohens_dz'])

    return stats

# semantic clustering score
# 1-sample t-test against mean of 0.5
def scl_statistics(scl_data_bsa):
    stats = []
    for (strat, c, ll, pr), data in scl_data_bsa.groupby(['strategy', 'condition', 'l_length', 'pres_rate']):
        res = scipy.stats.ttest_1samp(data.scl, popmean=0.5, nan_policy='omit', alternative='two-sided')
        cohens_d = cohens_d_1samp(data.scl, popmean=0.5)
        stats.append((strat, c, ll, pr, res.statistic, res.pvalue, res.df, cohens_d))
        
    stats = pd.DataFrame(stats, columns=['strategy', 'condition', 'l_length', 'pres_rate', 't_stat', 'p_val', 'dof', 'cohens_d'])
    
    # FDR correction
    fdr_pvals = scipy.stats.false_discovery_control(stats.p_val, method='by')
    stats['p_val_fdr'] = fdr_pvals
    
    # sort by condition and strategy
    stats = sort_by_condition_strategy(stats)
    
    return stats
