
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
import collections
import matplotlib

# This library is not complete yet, the parameters will be modified.

def _performSamplePreparation(data, minCounts, minGenes, minCells, maxPercentMito = .1):
    matplotlib.rcParams['figure.figsize'] = [3, 3]    
    data.var_names_make_unique()
    data.raw = data.copy()
    mito_genes = [name for name in data.var_names if name.startswith("MT-")]
    data.obs['percent_mito'] = np.ravel(np.sum(
    data[:, mito_genes].X, axis=1)) / np.ravel(np.sum(data.X, axis=1))
    data.obs['n_counts'] = np.ravel(np.sum(data.X, axis=1))
    sc.pp.filter_cells(data, min_genes=0) # Hack to generate the n_genes column
    sc.pp.filter_genes(data, min_cells=0)
    priorFilteringCellCount = data.n_obs
    priorFilteringGeneCount = len(data.var[data.var["n_cells"] > 0])

    sc.pp.filter_cells(data, min_genes=minGenes)
    sc.pp.filter_genes(data, min_cells=minCells) 
    
    data = data[data.obs['n_counts'] > minCounts, :]
    data = data[data.obs['percent_mito'] < maxPercentMito, :]
#   data_raw = sc.pp.log1p(data, copy=True)
    postFilteringCellCount = data.n_obs
    postFilteringGeneCount = len(data.var[data.var["n_cells"] > 0])
    sc.pp.normalize_per_cell(data, counts_per_cell_after=1e4)
    filter_result = sc.pp.filter_genes_dispersion(
        data.X, min_mean=0.0125, max_mean=3, min_disp=0.5)
    sc.pl.filter_genes_dispersion(filter_result)
    data = data[:, filter_result.gene_subset]
    hgvFilteringGeneCount = len(filter_result.gene_subset)
    hgvFilteringCellCount = len(data.var[data.var["n_cells"] > 0])
    sc.pp.log1p(data)
    sc.tl.pca(data)

    sc.pp.filter_genes(data, min_counts=1) # Remove unwanted zero count genes
    sc.pp.regress_out(data, ['n_counts'])
    sc.pp.scale(data, max_value=10)
    sc.tl.pca(data)
    sc.pl.pca(data)
    sc.pl.violin(data, ['n_genes', 'n_counts', 'percent_mito'], jitter=0.4, multi_panel=True, use_raw = False)       
    return data

def _cluster(data):
    sc.tl.pca(data)
    sc.pp.neighbors(data)
    sc.logging.print_memory_usage()
    sc.tl.umap(data)
    sc.logging.print_memory_usage()

def _preProcess(data):
    _performSamplePreparation(data, 500, 0, 3)
    _cluster(data)
    return data
    
def _UMAP(data, processed = True):
    if processed == False:
        _preProcess(data)
    matplotlib.rcParams['figure.figsize'] = [8,8]
    sc.pl.umap(data, color = 'cellType')
    return data
def _getDispersion(scData, geneList):
    try:
        matplotlib.rcParams['figure.figsize'] = [10, 10]  
        matplotlib.rcParams['axes.titlesize'] = 28
        matplotlib.rcParams['axes.labelsize'] = 20
        sc.pl.umap(scData, color = geneList,size=40,color_map="Blues")
    except:
        _preProcess(scData)
        matplotlib.rcParams['figure.figsize'] = [10, 10]  
        matplotlib.rcParams['axes.titlesize'] = 28
        matplotlib.rcParams['axes.labelsize'] = 20
        sc.pl.umap(scData, color = geneList,size=40,color_map="Blues")
    return scData
