"""
Created on Mon Feb  19 14:22:31 2019

@author: Murat Cem Köse
"""

import numpy as np
import pandas as pd
import scipy
import tqdm

def annotateCells(sc_data, refDataset, refAnnot, de):
    """ 
    The wrapper function to annotate cell types.
    
    Parameters
    ----------
    sc_data : DataFrame
        Single cell dataset gene expression matrix.
        
    refDataset : DataFrame
        Reference dataset gene expression matrix.
        
    refAnnot : DataFrame
        Annotations for samples in reference dataset.
        
    de : Multi Indexed DataFrame
        A table containing pairwise differentially expressed genes between subgroups.

    Returns
    -------
    round_res : DataFrame
        A matrix, including annotations of the subtypes.
        
    """
    unique_types = np.unique(refAnnot.cellType)
    n = 1
    combinations = [(i,j) for i in unique_types  for j in unique_types if i != j]
    pbar = tqdm.tqdm(total=100)

    while(n < len(unique_types)):
        if n == 1:
            type1 = unique_types[n-1]
            type2 = unique_types[n]
            round_res = paiwiseComparison(sc_data, refDataset, refAnnot, de, combinations, type1, type2)
            n+=1
        else:
            round_res = comparisonRound(sc_data, refDataset, refAnnot, unique_types, de, combinations, round_res, n)
            n+=1
        pbar.update(100/len(unique_types))
    pbar.update(100/len(unique_types))
    pbar.close()
    return round_res

def comparisonRound(sc_data, refDataset, refAnnot, unique_types, de, combinations, round_res, n):
    """ 
    A function to make comparison of correlations between the current annotation and next possible annotation.
    
    Parameters
    ----------
    sc_data : DataFrame
        Single cell dataset gene expression matrix.
        
    refDataset : DataFrame
        Reference dataset gene expression matrix.
        
    refAnnot : DataFrame
        Annotations for samples in reference dataset.
        
    unique_types : List
        A list containing unique cell types (subtypes) belonging to group of interest (maintype).
        
    de : Multi Indexed DataFrame
        A table containing pairwise differentially expressed genes between subgroups.
        
    combinations : List
        A list containing possible cell type combinations.
        
    round res : DataFrame
        The results form the previous round.
        
    n : Int
        The number of the current round.
        
    Returns
    -------
    temp : DataFrame
        A temporary matrix, including annotations of the subtypes for the current round.
        
    """
    groups = set(round_res.groupby("cellType").groups)
    temp = pd.DataFrame()
    for i in groups:
        type1 = i
        type2 = unique_types[n]
        comparison_res = paiwiseComparison(sc_data.loc[:,round_res[round_res.cellType == i].index], refDataset, refAnnot, de, combinations, type1, type2)
        temp = pd.concat([temp,comparison_res],axis = 0)
    return temp

def paiwiseComparison(sc_data, refDataset, refAnnot, de, combinations, type1, type2):
    """ 
    A function to make pairwise comparison of correlations between cell types.
    
    Parameters
    ----------
    sc_data : DataFrame
        Single cell dataset gene expression matrix.
        
    refDataset : DataFrame
        Reference dataset gene expression matrix.
        
    refAnnot : DataFrame
        Annotations for samples in reference dataset.
        
    de : Multi Indexed DataFrame
        A table containing pairwise differentially expressed genes between subgroups.

    combinations : List
        A list containing possible cell type combinations.
        
    type1 : Str
        First type to compare.
        
    type2 : Str
        Second type to compare.
        
    Returns
    -------
    comparison_res : DataFrame
        A temporary matrix, including annotations of the cell types for this comparison.
        
    """
    key = [i for i in de.columns if type1 in i and type2 in i][0]
    de_genes = de.loc[:,key].values

    cells = []
    cells.extend(refAnnot[refAnnot.cellType == type1].index.values)
    cells.extend(refAnnot[refAnnot.cellType == type2].index.values)
    
    cor = scipy.stats.spearmanr(sc_data.loc[de_genes], refDataset.loc[de_genes, cells])
    cor = pd.DataFrame(cor[0]).iloc[:,0:len(sc_data.columns)][-len(cells):].replace(float("NaN"),0)
    cor.columns = sc_data.columns
    cor["cellType"] = refAnnot.cellType.loc[cells].values
    comparison_res = pd.DataFrame(cor.groupby("cellType").quantile(q = 0.8).idxmax(), columns = ["cellType"])
    return comparison_res


