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

import scanpy as sc
def perform_label_transfer(ref_emb, query_emb, cell_type_column, k=15):
    # calculate an object representing the joing neighbor graph of ref + query
    ing = sc.tl.Ingest(ref_emb)
    ing.fit(query_emb)
    ing.neighbors(k=k)
    # calculate distances to top k neighbors for each cell and store indices
    # of neighbor cells
    top_k_distances, top_k_indices = (
        ing._distances,
        ing._indices,
    )
    # transform distances with Gaussian kernel (?)
    stds = np.std(top_k_distances, axis=1)
    stds = (2.0 / stds) ** 2  # don't know why the first 2.0
    stds = stds.reshape(-1, 1)
    top_k_distances_tilda = np.exp(-np.true_divide(top_k_distances, stds))
    # normalize so that transformed distances sum to 1
    weights = top_k_distances_tilda / np.sum(
        top_k_distances_tilda, axis=1, keepdims=True
    )
    # initialize empty series to store predicted labels and matching
    # uncertaintites for every query cell
    uncertainties = pd.Series(index=query_emb.obs_names, dtype="float64")
    pred_labels = pd.Series(index=query_emb.obs_names, dtype="object")
    # now loop through query cells
    y_train_labels = ref_emb.obs[cell_type_column].values
    for i in range(len(weights)):
        # store cell types present among neighbors in reference
        unique_labels = np.unique(y_train_labels[top_k_indices[i]])
        # store best label and matching probability so far
        best_label, best_prob = None, 0.0
        # now loop through all cell types present among the cell's neighbors:
        for candidate_label in unique_labels:
            candidate_prob = weights[
                i, y_train_labels[top_k_indices[i]] == candidate_label
            ].sum()
            if best_prob < candidate_prob:
                best_prob = candidate_prob
                best_label = candidate_label
        else:
            pred_label = best_label
        # store best label and matching uncertainty
        uncertainties.iloc[i] = max(1 - best_prob, 0)
        pred_labels.iloc[i] = pred_label
    # print info
    print(
        "Storing transferred labels in your query adata under .obs column:",
        f"transf_{cell_type_column}",
    )
    print(
        "Storing label transfer uncertainties in your query adata under .obs column:",
        f"transf_{cell_type_column}_unc",
    )
    # store results
    query_emb.obs[f"transf_{cell_type_column}"] = pred_labels
    query_emb.obs[f"transf_{cell_type_column}_unc"] = uncertainties
