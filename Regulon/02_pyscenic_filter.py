#2022.04.04
#pyscenic step2

#https://docs.python.org/3/library/getopt.html
##matrix shhould have this exac name, and gzipped
# import dependencies
import os
import numpy as np
import pandas as pd
import scanpy as sc
import loompy as lp
import pyscenic
import sys, getopt


#set workin dir
#wdir='/diskmnt/Projects/Users/ritalu/20220331_HTAN_BRCA/01_SCENIC/object'
wdir = '/diskmnt/Projects/Users/ritalu/20220331_HTAN_BRCA/02_pySCENIC/20220418/object/'
os.chdir( wdir )

# Set maximum number of jobs for Scanpy.
sc.settings.njobs = 40

# read unfiltered data from a loom file
adata = sc.read_loom('snRNA_obj.loom')

# compute the number of genes per cell (computes ‘n_genes' column)
sc.pp.filter_cells(adata, min_genes=0)

#
#nCells=adata.X.shape[0]
#cutoff=.001*nCells
cutoff=20
sc.pp.filter_genes(adata, min_cells=cutoff )
nCountsPerGene = np.sum(adata.X, axis=0)
nCellsPerGene = np.sum(adata.X>0, axis=0)
print("Number of counts (in the dataset units) per gene:", nCountsPerGene.min(), " - " ,nCountsPerGene.max())
print("Number of cells in which each gene is detected:", nCellsPerGene.min(), " - " ,nCellsPerGene.max())

# create basic row and column attributes for the loom file:
row_attrs = {
   "Gene": np.array(adata.var_names) ,
}
col_attrs = {
   "CellID": np.array(adata.obs_names) ,
   "nGene": np.array( np.sum(adata.X.transpose()>0 , axis=0)).flatten() ,
   "nUMI": np.array( np.sum(adata.X.transpose() , axis=0)).flatten() ,
}
lp.create( "snRNA_obj.v2.loom", adata.X.transpose(), row_attrs, col_attrs)

