#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from statsmodels.stats import multitest
import collections
from itertools import combinations

def getDEgenes(refDataset, refAnnot, refName, scName, n=20):
    """ 
    Creates a dictionary for differentially highly expressed genes for all pairwise cell types in the a reference data set.
    
    Parameters
    ----------
    refDataset : DataFrame
        Reference dataset gene expression matrix.
        
    refAnnot : DataFrame
        Annotations for samples in reference dataset.
        
    refName : String
        The name of the reference dataset.
    
    scName : String
        The name of the single cell dataset.
      
    n : Int
        Number of genes selecected for each comparison.
        
    Returns
    -------
    deGenes : DataFrame
        Table containing differentially highly expressed genes for each combination of cell types. 
    """
    deAll = {}
    try:
        de = pd.read_csv(refName+"_"+scName+"_DE.csv",index_col=0,header=[0, 1])
    except:
        de = rankDE(refDataset,refAnnot,n)
        de.to_csv(refName+"_"+scName+"_DE.csv")
    return de

def rankDE(refDataset,refAnnot,n):
    """ 
    Creates a data frame for differentially highly expressed genes for all pairwise cell types in the a reference data set.
    
    Parameters
    ----------
    refDataset : DataFrame
        Reference dataset gene expression matrix.
        
    refAnnot : DataFrame
        Annotations for samples in reference dataset.
        
    n : Int
        Number of genes selecected for each comparison.
        
    Returns
    -------
    deGenes : DataFrame
        Table containing differentially highly expressed genes for each combination of cell types. 
    """
    types = set(refAnnot.cellType)
    comb = combinations(types, 2) 
    deGenes={}
    for i in comb:
        rank1 = refDataset.loc[:,refAnnot[refAnnot.cellType == i[0]].index].rank(method = "min").median(axis = 1)
        rank2 = refDataset.loc[:,refAnnot[refAnnot.cellType == i[1]].index].rank(method = "min").median(axis = 1)
        rank_diff = rank1 - rank2
        rank_diff = rank_diff.sort_values(ascending = False)
        genes = []
        genes.extend(rank_diff.index.values[:n])
        genes.extend(rank_diff.index.values[-n:])
        deGenes.update({i:genes})
    return pd.DataFrame(deGenes)

