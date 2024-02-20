import pandas as pd

def map_clusters_to_genes(intron_meta, exons_table):
    """ Work out which gene each cluster belongs to. 
    
    Note the chromosome names used in the two inputs must match. 
    
    Parameter: 
        intron_meta Data.frame describing the introns, usually from get_intron_meta
        exons_table Data.frame of exons, see e.g. /data/gencode19_exons.txt.gz
    Returns: Data.frame with cluster ids and genes separated by commas
    """
    
    exons_table['temp'] = exons_table['start'].astype(int)
    intron_meta['temp'] = intron_meta['end'].astype(int)
    three_prime_matches = pd.merge(exons_table, intron_meta, on = ['chr', 'temp'], how = 'inner')

    exons_table['temp'] = exons_table['end'].astype(int)
    intron_meta['temp'] = intron_meta['start'].astype(int)
    five_prime_matches = pd.merge(exons_table, intron_meta, on = ['chr', 'temp'], how = 'inner')

    gene_df = pd.concat([three_prime_matches, five_prime_matches]).drop_duplicates().loc[:,('chr', 'cluster', 'gene_name')]
    
    if gene_df.shape[0] == 0:
        return None

    clu_df = gene_df.groupby('cluster')['gene_name'].apply(lambda x: ','.join(list(set(x)))).reset_index()
    clu_df.rename(columns={'gene_name': 'genes'}, inplace=True)
    return clu_df

def get_intron_meta(introns):
    """ Make a data.frame of meta data about the introns
    
    introns Names of the introns
    Returns: Data.frame with chr, start, end, cluster_name
    """
    intron_meta = pd.Series(introns).str.split(":", expand = True)
    intron_meta.columns = ['chr', 'start', 'end', 'cluster']
    return intron_meta

def add_chr(chrs):
    """Add `chr` as a prefix if not currently there. """ 
    
    if not ('chr' in str(chrs[0])): 
        chrs = 'chr' + chrs
    return chrs


