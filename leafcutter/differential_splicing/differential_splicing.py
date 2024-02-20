import pandas as pd
import numpy as np
import torch
from collections import defaultdict, OrderedDict

from leafcutter.differential_splicing.dm_glm import dirichlet_multinomial_anova, SimpleGuide, CleverGuide

try: 
    from scipy.stats import false_discovery_control
except: 
    print("Warning: your scipy version is older than 1.11.0 so scipy.stats.false_discovery_control isn't available.")
    def false_discovery_control(ps, method = "bh"):
        order = np.argsort(ps) # put p values in ascending order
        ps = ps[order]
        m = len(ps)
        i = np.arange(1, m+1)
        ps *= m / i # adjust
        if method == 'by': ps *= np.sum(1. / i)
        ps = ps[np.argsort(order)] # put back in original order
        return np.clip(ps, 0, 1)
    
def robust_fdr(ps, method = "bh"): 
    mask = ~np.isnan(ps)
    p_adjust = ps.copy()
    p_adjust[mask] = false_discovery_control(ps[mask], method = method)
    return p_adjust

def differential_splicing(counts, x, confounders = None, max_cluster_size=10, min_samples_per_intron=5, min_samples_per_group=4, min_coverage=0, min_unique_vals = 10, device = "cpu", **kwargs):
    '''Perform pairwise differential splicing analysis.

    counts: An [introns] x [samples] dataframe of counts. The rownames must be of the form chr:start:end:cluid. If the counts file comes from the leafcutter clustering code this should be the case already.
    x: A [samples] numeric pandas Series, should typically be 0s and 1s, although in principle scaling shouldn't matter.
    confounders: A [samples] x [confounders] pandas dataframe to be controlled for in the GLM. Factors should already have been converted to a 1-of-(K-1) encoding
    ####, e.g. using model.matrix (see scripts/leafcutter_ds.R for how to do this). Can be None, implying no covariates are controlled for.
    max_cluster_size: Don't test clusters with more introns than this (Default = 10)
    min_samples_per_intron: Ignore introns used (i.e. at least one supporting read) in fewer than n samples (Default = 5)
    min_samples_per_group: Require this many samples in each group to have at least min_coverage reads (Default = 4)
    min_coverage: Require min_samples_per_group samples in each group to have at least this many reads (Default = 20)
    kwargs: keyword arguments passed to dirichlet_multinomial_anova
    '''
    
    normalize = lambda g: g/g.sum()

    results = OrderedDict()
    junc_results = OrderedDict()
    
    statuses = OrderedDict()

    junc_meta = counts.index.to_series().str.split(':',expand=True).rename(columns = {0:"chr", 1:"start", 2:"end", 3:"cluster"})
    cluster_ids = junc_meta.cluster.unique()
    
    torch_types = { "device" : device, "dtype" : torch.float } # would we ever want float64? 
    
    for clu in cluster_ids:
        idx = clu == junc_meta.cluster
        cluster_size = idx.sum()

        if cluster_size>max_cluster_size: 
            statuses[clu] = "Too many introns in cluster"
            continue
        if cluster_size <= 1: 
            statuses[clu] = "<=1 junction in cluster"
            continue
        cluster_counts=np.array(counts.loc[ idx,: ]).transpose()
        sample_totals=cluster_counts.sum(1)
        samples_to_use=sample_totals>0
        if samples_to_use.sum()<=1: 
            statuses[clu] = "<=1 sample with coverage>0"
            continue
        sample_totals=sample_totals[samples_to_use]
        if (sample_totals>=min_coverage).sum()<=1:
            statuses[clu] = "<=1 sample with coverage>min_coverage"
            continue
        # this might be cleaner using anndata
        x_subset=x[samples_to_use] # assumes one covariate? ah no, covariates handled later 
        cluster_counts=cluster_counts[samples_to_use,]
        introns_to_use=(cluster_counts>0).sum(0)>=min_samples_per_intron # only look at introns used by at least 5 samples
        if introns_to_use.sum()<2:
            statuses[clu] = "<2 introns used in >=min_samples_per_intron samples"
            continue
        cluster_counts=cluster_counts[:,introns_to_use]
        
        # this is the only part that depends on x
        unique_vals, ta = np.unique(x_subset[sample_totals>=min_coverage], return_counts = True)
        if x_subset.dtype.kind in 'OUS': # categorical x
            if (ta >= min_samples_per_group).sum()<2: # at least two groups have this
                statuses[clu] = "Not enough valid samples"
                continue
            x_subset = pd.get_dummies(x_subset, drop_first = True)
            to_drop = (x_subset == 0).all()
            x_subset = x_subset.loc[:, ~to_drop ] # remove empty groups
        else: # continuous x
            if len(unique_vals) < min_unique_vals:
                statuses[clu] = "Not enough valid samples"
                continue
            x_subset = pd.DataFrame( {"x" : x_subset} )
        
        x_only = torch.tensor(x_subset.to_numpy(), **torch_types)
        intercept = torch.ones(x_only.shape[0], 1, **torch_types) 
        if confounders is None:
            x_full = torch.cat((intercept, x_only), axis = 1)
            x_null = intercept
        else:
            these_confounders = torch.tensor(confounders.iloc[samples_to_use].to_numpy(), **torch_types)
            #filter out confounders with no standard deviation
            these_confounders = these_confounders[:,torch.std(these_confounders, dim = 0) != 0.]
            x_full = torch.cat((intercept, these_confounders, x_only), axis = 1)
            x_null = torch.cat((intercept, these_confounders), axis = 1)
        
        y = torch.tensor(cluster_counts, **torch_types)

        #np.savetxt("cached_xy/x_full_%s.tsv" % clu, x_full.numpy(), delimiter="\t")
        #np.savetxt("cached_xy/x_null_%s.tsv" % clu, x_null.numpy(), delimiter="\t")
        #np.savetxt("cached_xy/y_%s.tsv" % clu, y.numpy(), delimiter="\t")

        #if (clu != 'clu_711_NA'): continue
        try: 
            loglr, df, lrtp, null_fit, full_fit, refit_null_flag = dirichlet_multinomial_anova(
                x_full, 
                x_null, 
                y, 
                **kwargs)
        except ValueError as value_error: # catch rare numerical errors
            statuses[clu] = "ValueError: " + str(value_error).replace('\n', ' ')
            continue
        
        statuses[clu] = "Success"

        results[clu] = {
            'loglr': loglr,
            'null_ll': -null_fit.loss,
            'null_exit_status': null_fit.exit_status,
            'full_ll': -full_fit.loss,
            'full_exit_status': full_fit.exit_status,
            'df': df,
            'p': lrtp
        }
        
        x_dim = x_subset.shape[1]
        P_null = null_fit.beta.shape[0]
        
        # extract effect sizes
        logef = pd.DataFrame( full_fit.beta[-x_dim:,:].cpu().numpy().T, columns = x_subset.columns ).add_prefix('ef_')
        
        # calculate model-based PSI for each category. For continuous x this will correspond to being +1s.d. from the mean. 
        perturbed = torch.stack( [ normalize( (full_fit.beta[0,:] + full_fit.beta[P_null+i,:]).softmax(0) * full_fit.conc) for i in range(x_dim) ]).T.cpu().numpy()
        perturbed = pd.DataFrame( perturbed, columns = x_subset.columns ).add_prefix('psi_')
        
        junc_results[clu] = pd.concat(
            [pd.DataFrame({
                'cluster' : [clu] * y.shape[1], 
                'intron' : idx[idx][introns_to_use].index, 
                'psi_0' : normalize(full_fit.beta[0,:].softmax(0) * full_fit.conc).cpu().numpy()}), 
            logef, 
            perturbed], axis = 1)      
    
    status_df = pd.DataFrame(statuses.values()) 
    status_df.index = statuses.keys()
    
    cluster_table = pd.DataFrame(results.values()) 
    cluster_table.index = results.keys()
    cluster_table['p.adjust'] = robust_fdr(cluster_table['p'], method = 'bh')
    
    junc_table = pd.concat(junc_results.values(), axis=0) # note this should handle missing categories fine
    #junc_table["deltapsi"] = junc_table['perturbed'] - junc_table['baseline']
    
    return cluster_table, junc_table, status_df