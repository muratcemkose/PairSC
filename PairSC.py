"""
Created on Mon Feb  19 14:22:31 2019

@author: Murat Cem Köse
"""

import numpy as np
import pandas as pd
import Utils
import Annotate
import DE
import ScProcess

class createObject:
    def __init__(self, refDataset, refAnnot, refName):
        """
        Contructor function for HierarchicalSC class.
    
        Parameters
        ----------            
        refDataset : DataFrame
            The reference dataset gene expression matrix. Columns representing sample names, rows representing gene symbols.
            
        refAnnot : DataFrame
            Annotations for samples in reference dataset. This data frame should contain only one column that represents cell types. The rows represents cell names in reference dataset. 
        
        refName : String
            Name of the reference dataset.
        """

        self.refDataset = refDataset.astype(float)
        self.refAnnot = refAnnot
        self.refAnnot.columns = ["cellType"]
        self.refName = refName
        
    def addScData(self, scData, scName):
        """
        Contructor function for HierarchicalSC class.
    
        Parameters
        ----------
        scData : DataFrame
            Single cell data matrix. Columns representing sample names, rows representing gene symbols.
        
        scName : String
            Name of the single cell dataset.
        """

        self.scData = scData
        self.scName = scName
        self.scData.var_names = self.scData.var_names.str.upper()
        self.refDataset.index = self.refDataset.index.str.upper()
        self._genesIntersect = np.intersect1d(self.scData.var_names, self.refDataset.index)
        scData = self.scData[:,self._genesIntersect].copy()
        refDataset = self.refDataset.loc[self._genesIntersect].copy()
        self.deAll = DE.getDEgenes(refDataset=refDataset, refAnnot=self.refAnnot, refName = self.refName, scName = self.scName)
    def getCellTypes(self):
        """
        Annotates single cell types at each level and adds the result to the object.
        
        """
        scData = Utils.convertAnnDataToDf(self.scData)
        scData = scData.loc[self._genesIntersect]
        refDataset = self.refDataset.loc[self._genesIntersect]
   
        self.scAnnot = Annotate.annotateCells(scData, refDataset, self.refAnnot, self.deAll)
        self.scData.obs["cellType"] = self.scAnnot["cellType"]

    def getDriverGenes(self):
        """
        Finds driver genes and adds the result to the object.
        
        """
        try:
            self.scAnnot
        except:
             print("Please run getCellTypes first to get cell annotations. This step is needed for driver gene finding.")
        self.driverGenes = Utils._findDrivers(self.scAnnot, self.deAll)
            
    def getUMAP(self, processed = True):
        """
        Gerating UMAP for annotated data.
        
        """

        self.scData = ScProcess._UMAP(self.scData, processed)

    def getDriverDispersion(self, cellType = None, geneList = None):
        """Plots gene expression dispersion. Either a cell type for driver gene expression or a gene list can be given. 

        Parameters
        ----------
        scData : DataFrame
            Single cell data. 
            
        cellType : String
            The cell type for driver genes.
            
        geneList : List
            The list of genes.
        """
        if cellType != None:
            geneList = self.driverGenes.get(cellType)
        self.scData = ScProcess._getDispersion(self.scData, geneList)
        