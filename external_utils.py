from sccoda.util import data_visualization as viz
from sccoda.util import cell_composition_data as dat
import anndata as ad
import pandas as pd

def plot_proportions(adata, label_category, groupby, **kwargs):
    """Create an sccoda plot of proportions of labels in a groupby.
    
    Params
    ------
    adata : ad.AnnData
    label_category : str
        Categorical label from which to plot proportions.
    groupby : str
    kwargs
        Passed into `sccoda.util.data_visualization.stacked_barplot`.

    Returns
    -------
    df_adata : ad.AnnData
    """
    df_adata = ad.AnnData(pd.crosstab(adata.obs[groupby], adata.obs[label_category]))
    df_adata.obs = df_adata.obs.reset_index()

    g = viz.stacked_barplot(df_adata, feature_name=groupby, **kwargs)
    g.set_xticklabels(g.get_xticklabels(), rotation=90, ha='right')

    return df_adata

import numpy as np
import matplotlib.pyplot as plt

def define_marker_genes(adata, clf, n_markers=5000):
    """Add marker genes to adata.var according to sklearn model coefficients.
    
    Params
    ------
    adata : AnnData
    clf : `sklearn.linear_model`
        `clf.coef_` must be defined.
    n_markers : int (default: 5000)
    """
    # distinguish between linear and logistic regression
    if len(clf.coef_.shape) == 1:  # linear
        coef = clf.coef_
    else:
        coef = clf.coef_[0]
    # quantiles required to reach n_markers genes
    q = n_markers/len(coef)/2
    
    # plotting
    plt.hist(coef, bins=50)
    upper = np.quantile(coef, 1-q)
    lower = np.quantile(coef, q)
    plt.axvline(upper, c='cyan')
    plt.axvline(lower, c='cyan')

    plt.yscale('log')
    plt.xlabel('model coefficients')
    plt.grid(b=None)
    
    # add to .var
    adata.var['coef'] = coef
    adata.var['marker_genes'] = [v < lower or v > upper for v in coef]
