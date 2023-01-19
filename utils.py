import pandas as pd
from gprofiler import GProfiler
def enrich(query, background, return_full=False, organism='hsapiens', sources=['GO:BP']):
    gp = GProfiler(return_dataframe=True, user_agent='g:GOSt')

    df = gp.profile(
        organism=organism, sources=sources, user_threshold=0.05,
        significance_threshold_method='fdr', 
        background=list(background),
        query=list(query), 
        no_evidences=False)

    if return_full:
        return df
    else:
        return df[['name', 'p_value', 'intersections']]

# Enrichr setup
import json
import requests

def enrichr(genes, library='KEGG_2019_Human'):
    ENRICHR_URL = 'http://maayanlab.cloud/Enrichr/addList'
    genes_str = '\n'.join(list(genes))
    # name of analysis or list
    description = 'enrichment'
    
    ##run enrichment
    payload = {'list': (None, genes_str), 'description': (None, description)}
    response = requests.post(ENRICHR_URL, files=payload)
    if not response.ok:
        raise Exception('Error analyzing gene list')
    job_id = json.loads(response.text)

    ## get enrichment results
    ENRICHR_URL = 'http://maayanlab.cloud/Enrichr/enrich'
    query_string = '?userListId=%s&backgroundType=%s'
    response = requests.get(
        ENRICHR_URL + query_string % (str(job_id['userListId']), library)
     )
    if not response.ok:
        raise Exception('Error fetching enrichment results')
    print('Enrichr API : Get enrichment results: Job Id:', job_id)
    
    # construct readable dataframe
    out = pd.DataFrame(list(json.loads(response.text).values())[0])
    out.columns = [
        'Rank',
        'Term name',
        'P-value',
        'Z-score',
        'Combined score',
        'Overlapping genes',
        'Adjusted p-value',
        'Old p-value',
        'Old adjusted p-value'
    ]
    return out


import scipy
import numpy as np
import matplotlib.pyplot as plt

def plot_correlation(df, col1, col2, margin_thresh = 30, expression_thresh = 25, plot_legend=False):
    """Scatterplot of two columns of a dataframe containing scores or pvalues.
    Points with significant correlation are colored in blue.

    Params
    ------
    df : pd.DataFrame
        Dataframe with numerical values to be plotted.
    col1 : str
        x-axis values.
    col2 : str
        y-axis values.
    margin_thresh : int (default: 30)
        Threshold by which to color significant correlations blue.
        The default value expects overinflated -log10(pvalue).
    expression_thresh : int (default: 25)
        Threshold by which to remove insignificant correlations
        depending on their absolute difference.
    plot_legend : bool (default: False)
        Whether or not to plot the pearsonr and correlation pvalue.

    Returns
    -------
    A subset of the original dataframe containing only the highly correlated rows.
    """
    from adjustText import adjust_text
    p_r2 = scipy.stats.pearsonr(df[col1], df[col2])
    plt.scatter(df[col1], df[col2], c='grey',
        label='R2 = {:.2f}\n'.format(p_r2[0]))

    # calculate high correlation points
    df['margin'] = np.abs(df[col1] - df[col2])
    # df['corr_score'] = np.abs(np.cross(  # for when the fit line is not y=x
    #     np.array([1000, 1000]),
    #     df[[col1, col2]].values
    # ))
    df['expr'] = np.abs(np.mean(df[[col1, col2]], axis=1))

    # plot high correlation points
    high_corr = df[(df.margin < margin_thresh) & (df.expr > expression_thresh)]
    plt.scatter(high_corr[col1], high_corr[col2], color='blue')
#    texts = [plt.text(c2, c3, c1, size='x-small') for i, (c1, c2, c3) in high_corr[[col1, col2]].reset_index().iterrows() if i < 20]
#    adjust_text(texts)

    # touch-ups
    plt.axline(xy1=(0, 0), slope=1)
    if plot_legend:
        plt.legend(handlelength=0, markerscale=0)

    return high_corr
