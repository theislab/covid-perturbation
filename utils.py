import pandas as pd
from gprofiler import GProfiler
def enrich(query, background, return_full=False):
    gp = GProfiler(return_dataframe=True, user_agent='g:GOSt')

    df = gp.profile(
        organism='hsapiens', sources=['GO:BP'], user_threshold=0.05,
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
