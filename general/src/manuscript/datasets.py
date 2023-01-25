import os

import numpy as np
import pandas as pd

from manuscript import inout

import sys
sys.path.append(inout.gustav_path())
from gustav import ncbi


def _get_material_path(user, p, date=False):
    """
    Takes extension path p, and makes it as a subfile
    within material folder
    """

    p = os.path.join(
        inout.get_internal_path(f'materials/{user}'),
        p)

    if date:
        [fo, fn] = os.path.split(p)
        [fb, ext] = os.path.splitext(fn)
        dt = datetime.datetime.today().strftime('%y%m%d_%H%M')
        p = os.path.join(fo, fb + '_' + dt + ext)

    return p



def reference_genes(taxon, flavor):
    """
    Standardized sets of genes to consider:
    
     Input:
        taxon   ncbi taxonomy ID (e.g.: 9606 for humans)
        flavor  e.g.: 'orp' where
                      o will require official gene symbol
                      p will limit to protein-coding
                      u will require unambiguous 1:1 mapping with ensembl
    
    Output:
        ref_genes   list
    
    """
    
    
    columns = ['gene_ncbi', 'dbxrefs', 'symbol_from_nomenclature_authority', 'type_of_gene', 'taxon_ncbi']
      
    if 'u' in flavor:
        gi = ncbi.gene_info(taxon, mode='unambiguous_ensembl', usecols=columns)
    else:
        gi = ncbi.gene_info(taxon, usecols=columns)

    if 'o' in flavor:
        f = gi['symbol_from_nomenclature_authority'] == '-'
        gi = gi.loc[~f, :].copy()

    if 'p' in flavor:
        f = gi['type_of_gene'] == 'protein-coding'
        gi = gi.loc[f, :].copy()
        
    return sorted(list(set(gi['gene_ncbi'])))


def reference_gene2lit(taxon, flavor):
    """
    Mapping between papers and genes
    
    Input:
        taxon   ncbi taxonomy ID (e.g.: 9606 for humans)
        flavor  One of following options:
                    gene2pubmed
                    pubtator_title_or_abstract
                    pubtator_title_or_abstract_in_any_gene2pubmed_paper
    
    Output:
        gene2lit   df with ['pubmed_id', 'gene_ncbi']
    """
    
    def _from_pubtator(taxon):
        pubtator = ncbi.pubtator_medline(       
                'pooled_counts_gene_unambiguous')

        gene2pubtator = pubtator[
            (pubtator['bioconcept_title_or_abstract']>0)][
            ['pubmed_id', 'gene_ncbi']
        ].drop_duplicates()

        gi = ncbi.gene_info('all', ['taxon_ncbi', 'gene_ncbi'])

        gene2pubtator = gene2pubtator[gene2pubtator['gene_ncbi'].isin(
            gi[gi['taxon_ncbi']==taxon]['gene_ncbi']
        )]
        return gene2pubtator
        
    
    if flavor=='gene2pubmed':
        g2lit = ncbi.gene2pubmed(taxon)[
            ['pubmed_id', 'gene_ncbi']].drop_duplicates().copy()
    elif flavor=='pubtator_title_or_abstract':
        g2lit = _from_pubtator(taxon)
    elif flavor=='pubtator_title_or_abstract_in_any_gene2pubmed_paper':
        g2lit = _from_pubtator(taxon)
                
        g2lit = g2lit[g2lit['pubmed_id'].isin(
            ncbi.gene2pubmed('all')['pubmed_id']
        )].copy()
        
        
    return g2lit.reset_index(drop=True)
    
      

def reference_publications(taxon='any'):
    """
    Will import reference publications
    
    Input:
        taxon   str, default: 'any', other: 
                NCBI taxonomy code (e.g.: 9606 for human)
        
    Output:
        pubmed_ids set
    
    """

    taxa_to_dataset = {
        559292: 'scerevisiae',
        6239: 'celegans',
        7227: 'dmelanogaster',
        7955: 'drerio',
        9606: 'human',
        10090: 'mouse',
        10116: 'rat'
    }
    
    
    if taxon=='any':
        file = 'references/literature.parquet'
    elif taxon in taxa_to_dataset.keys():
        name = taxa_to_dataset[taxon]
        file = f'references/literature_{name}.parquet'
    else:
        raise AssertionError('taxon not supported.')
    
    p = _get_material_path('general', file)
    ref_publications = set(pd.read_parquet(p)['pubmed_id'])
    return ref_publications
