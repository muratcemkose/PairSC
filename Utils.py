
"""
Created on Mon Dec  3 14:45:33 2018

@author: Murat Cem Köse
"""

import os
import scanpy.api as sc
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib
import collections

def readCountMartix(path,min_genes):
    """ 
    Reads, precesses and returns single cell data. 
    
    Parameters
    ----------
    path : Str
        Diractory path, the location of single cell data.
        
    min_genes : Int
        The minimum number of genes for a cell to have in order to participate the analysis.
        
    Returns
    -------
    scData : AnnData
        Single cell data. 
        
    """
    
    result = sc.read_csv(path).transpose() #, cache=True
    result.var_names_make_unique()
    n_counts = np.sum(result.X, axis=1)
    result.obs['n_counts'] = n_counts
    sc.pp.filter_cells(result, min_genes=min_genes)
    
    return result

def read10xData(path,min_genes):
    """ 
    Reads, precesses and returns single cell data. 
    
    Parameters
    ----------
    path : Str
        Diractory path, the location of single cell data.
        
    min_genes : Int
        The minimum number of genes for a cell to have in order to participate the analysis.
        
    Returns
    -------
    scData : AnnData
        Single cell data. 
        
    """
    
    result = sc.read(path + 'matrix.mtx').transpose() #, cache=True
    result.var_names = np.genfromtxt(path + 'genes.tsv', dtype=str)[:, 1]
    result.obs_names = np.genfromtxt(path + 'barcodes.tsv', dtype=str)
    result.var_names_make_unique()
    result.obs['n_counts'] = np.sum(result.X, axis=1).A1
    sc.pp.filter_cells(result, min_genes=min_genes)
    
    return result
def readSCData(path,min_genes=500):
    """ 
    Reads, precesses and returns single cell data. 
    
    Parameters
    ----------
    path : Str
        Diractory path, the location of single cell data.
        
    min_genes : Int
        The minimum number of genes for a cell to have in order to participate the analysis.
        
    Returns
    -------
    scData : AnnData
        Single cell data. 
        
    """
    
    if "matrix.mtx" in os.listdir(path):
        result = read10xData(path, min_genes)
        return result
    else:
        try:
            data_file = [i for i in os.listdir(path) if "expression" in i][0]
            result = readCountMartix(path+data_file, min_genes)
            return result
        except:
            print("Please end your single cell data count matrix file with word 'expression'.")

def convertAnnDataToDf(scData):
    """ 
    Converts AnnData object obtained from scanpy into a pandas dataframe.
    
    Parameters
    ----------
    scData : AnnData
        Single cell data. 
        
    Returns
    -------
    result : DataFrame
        Single cell data as pandas dataframe. 
        
    """
    try:
        result = pd.DataFrame(scData.X.toarray()) # If data is 10x data
    except:
        result = pd.DataFrame(scData.X[:]) # If data is Digital Gene expression matrix
        
    result.index = scData.obs_names.values
    result.columns = scData.var_names.values
    return result.T


def _findDrivers(scAnnot, deDict):
    """ 
    Finds drivers genes for each cell type at each level.
    
    Parameters
    ----------
    sc_annot : DataFrame
        Single cell annotations. 
    de_dict : Dict
        Dictionary containing DE genes generated while annotating cell types.
    Returns
    -------
    driver_genes : Dict
        A dictionary including driver genes for cell types based on level. 
        
    """
    cand_lvl = pd.DataFrame()
    unique_types = [i for i in scAnnot.groupby("cellType").groups]
    for type1 in unique_types:
        keys = [i for i in deDict.keys() if type1 in i]
        de_all = []
        [de_all.extend(deDict.get(i)) for i in keys]
        type1_drivers = pd.DataFrame(collections.Counter(de_all),index=[type1]).T.sort_values(by=type1, ascending = False) 
        cand_lvl = pd.concat([cand_lvl, type1_drivers],axis =1).replace(float("NaN"),0)
    driver_genes = {}
    for type1 in unique_types:
        candidates = cand_lvl[cand_lvl[type1]> np.sqrt(len(cand_lvl.columns))].index.values
        temp = cand_lvl.drop(columns=type1).loc[candidates]
        driver_genes.update({type1:temp[temp < np.sqrt(len(cand_lvl.columns))].dropna(how = "any").index.values})
    return driver_genes


    
    
    
