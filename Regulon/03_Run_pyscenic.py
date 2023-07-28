#!/usr/bin/env python
# coding: utf-8

# new loom with cellinfo
# 

# Ref: https://github.com/aertslab/SCENICprotocol/blob/master/notebooks/PBMC10k_SCENIC-protocol-CLI.ipynb

# In[2]:


import os
import numpy as np
import pandas as pd
import scanpy as sc
import loompy as lp
from sklearn.manifold import TSNE


# In[9]:


wdir='/diskmnt/Projects/Users/ritalu/20220331_HTAN_BRCA/02_pySCENIC/20220418/object/'
os.chdir( wdir )


# In[20]:


# path to unfiltered loom file (this will be created in the optional steps below)
f_loom_path_unfilt = "snRNA_obj.loom" 

# path to loom file with basic filtering applied (this will be created in the "initial filtering" step below). Optional.
f_loom_path_scenic = "snRNA_obj_v2.loom"

# path to anndata object, which will be updated to store Scanpy results as they are generated below
f_anndata_path = "anndata.h5ad"

# path to pyscenic output
f_pyscenic_output = "snRNA_pyscenic_output.loom"

# loom output, generated from a combination of Scanpy and pySCENIC results:
f_final_loom = 'pyscenic_output_integrated-output.loom'



# ## Initial/basic filtering

# In[34]:


adata = sc.read_loom(f_loom_path_unfilt)


# In[35]:


sc.pp.filter_cells(adata, min_genes=0)


# In[36]:


cutoff=20
sc.pp.filter_genes(adata, min_cells=cutoff )
nCountsPerGene = np.sum(adata.X, axis=0)
nCellsPerGene = np.sum(adata.X>0, axis=0)
print("Number of counts (in the dataset units) per gene:", nCountsPerGene.min(), " - " ,nCountsPerGene.max())
print("Number of cells in which each gene is detected:", nCellsPerGene.min(), " - " ,nCellsPerGene.max())


# In[43]:


adata.obs['cell_type_specific']


# ### Finalize the selected filters
# Update the anndata file, to be used for further processing, clustering, visualization, etc..
# 
# 

# In[40]:


row_attrs = {
   "Gene": np.array(adata.var_names) ,
}


# In[44]:


col_attrs = {
   "CellID": np.array(adata.obs_names) ,
    
   "nGene": np.array( np.sum(adata.X.transpose()>0 , axis=0)).flatten() ,
   "nUMI": np.array( np.sum(adata.X.transpose() , axis=0)).flatten() ,
}


# In[46]:


np.array(adata.obs['nGene'])


# In[47]:


lp.create( f_loom_path_scenic, adata.X.transpose(), row_attrs, col_attrs)


# In[ ]:





# In[ ]:


# # create basic row and column attributes for the loom file:
# row_attrs = {
#    "Gene": np.array(adata.var_names) ,
# }
# col_attrs = {
#    "CellID": np.array(adata.obs_names) ,
#    "nGene": np.array( np.sum(adata.X.transpose()>0 , axis=0)).flatten() ,
#    "nUMI": np.array( np.sum(adata.X.transpose() , axis=0)).flatten() ,
# }
# lp.create( f_loom_path_scenic, adata.X.transpose(), row_attrs, col_attrs)


# ## Further pre-processing of expression data

# In[15]:


# # read filtered data from a loom file
adata = sc.read_loom( f_loom_path_scenic )
adata


# ### 1 PCA

# In[62]:


# adata = sc.read_h5ad( f_anndata_path )
# principal component analysis
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata, log=True)
adata.write( f_anndata_path )


# ### 2 Visualization of highly variable genes

# In[63]:


# neighborhood graph of cells (determine optimal number of PCs here)
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
# compute UMAP
sc.tl.umap(adata)
# tSNE
tsne = TSNE( n_jobs=20 )
adata.obsm['X_tsne'] = tsne.fit_transform( adata.X )
adata.write( f_anndata_path )


# In[65]:


adata


# ### 3 Clustering
# 

# In[64]:


# cluster the neighbourhood graph
sc.tl.louvain(adata,resolution=0.4)

sc.pl.umap(adata, color=['louvain'] )


# In[66]:


# find marker genes
sc.tl.rank_genes_groups(adata, 'louvain', method='t-test')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)

# sc.tl.rank_genes_groups(adata, 'louvain', method='logreg')
# sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(10)
adata.write( f_anndata_path )


# In[30]:


wdir


# # SCENIC steps
# ## STEP 1: Gene regulatory network inference, and generation of co-expression modules
# ### Phase Ia: GRN inference using the GRNBoost2 algorithm
# 

# In[32]:


# pyscenic grn \
# /diskmnt/Projects/Users/ritalu/20220331_HTAN_BRCA/02_pySCENIC/20220418/object/snRNA_obj.v2.loom \
# /diskmnt/Projects/Users/ritalu/20220331_HTAN_BRCA/01_SCENIC/Resource/JAPSAR2020_unique_motifs.txt \
# -o /diskmnt/Projects/Users/ritalu/20220331_HTAN_BRCA/02_pySCENIC/20220418/object/adj.csv \
# --num_workers 45 --sparse 


# In[ ]:


pyscenic ctx /diskmnt/Projects/Users/ritalu/20220331_HTAN_BRCA/02_pySCENIC/20220418/object/adj.csv     /diskmnt/Projects/Users/ritalu/20220331_HTAN_BRCA/01_SCENIC/Resource/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather /diskmnt/Projects/Users/ritalu/20220331_HTAN_BRCA/01_SCENIC/Resource/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather     --annotations_fname /diskmnt/Projects/Users/ritalu/20220331_HTAN_BRCA/01_SCENIC/Resource/motifs-v9-nr.hgnc-m0.001-o0.0.tbl     --expression_mtx_fname /diskmnt/Projects/Users/ritalu/20220331_HTAN_BRCA/02_pySCENIC/20220418/object/snRNA_obj.v2.loom     --output /diskmnt/Projects/Users/ritalu/20220331_HTAN_BRCA/02_pySCENIC/20220418/object/reg.csv     --mask_dropouts     --num_workers 50


# In[ ]:


pyscenic aucell    /diskmnt/Projects/Users/ritalu/20220331_HTAN_BRCA/02_pySCENIC/20220418/object/snRNA_obj.v2.loom    /diskmnt/Projects/Users/ritalu/20220331_HTAN_BRCA/02_pySCENIC/20220418/object/reg.csv    --output /diskmnt/Projects/Users/ritalu/20220331_HTAN_BRCA/02_pySCENIC/20220418/object/snRNA_pyscenic_output.loom    --num_workers 50


# In[ ]:





# In[ ]:





# ## Visualization of SCENIC's AUC matrix

# In[67]:


import json
import zlib
import base64

# collect SCENIC AUCell output
lf = lp.connect( f_pyscenic_output, mode='r+', validate=False )
auc_mtx = pd.DataFrame( lf.ca.RegulonsAUC, index=lf.ca.CellID)
lf.close()


# In[68]:


import umap

# UMAP
runUmap = umap.UMAP(n_neighbors=10, min_dist=0.4, metric='correlation').fit_transform
dr_umap = runUmap( auc_mtx )
pd.DataFrame(dr_umap, columns=['X', 'Y'], index=auc_mtx.index).to_csv( "scenic_umap.txt", sep='\t')
# tSNE
tsne = TSNE( n_jobs=20 )
dr_tsne = tsne.fit_transform( auc_mtx )
pd.DataFrame(dr_tsne, columns=['X', 'Y'], index=auc_mtx.index).to_csv( "scenic_tsne.txt", sep='\t')


# ### Integrate the output
# Here, we combine the results from SCENIC and the Scanpy analysis into a SCope-compatible loom file
# 
# 

# In[108]:


# scenic output
lf = lp.connect( f_pyscenic_output, mode='r+', validate=False )
meta = json.loads(zlib.decompress(base64.b64decode( lf.attrs.MetaData )))
#exprMat = pd.DataFrame( lf[:,:], index=lf.ra.Gene, columns=lf.ca.CellID)
auc_mtx = pd.DataFrame( lf.ca.RegulonsAUC, index=lf.ca.CellID)
regulons = lf.ra.Regulons
dr_umap = pd.read_csv( 'scenic_umap.txt', sep='\t', header=0, index_col=0 )
dr_tsne = pd.read_csv( 'scenic_tsne.txt', sep='\t', header=0, index_col=0 )
###


# Fix regulon objects to display properly in SCope:

# In[109]:


auc_mtx.columns = auc_mtx.columns.str.replace('\(','_(')
regulons.dtype.names = tuple( [ x.replace("(","_(") for x in regulons.dtype.names ] )
# regulon thresholds
rt = meta['regulonThresholds']
for i,x in enumerate(rt):
    tmp = x.get('regulon').replace("(","_(")
    x.update( {'regulon': tmp} )


# Concatenate embeddings (tSNE, UMAP, etc.)

# In[110]:


tsneDF = pd.DataFrame(adata.obsm['X_tsne'], columns=['_X', '_Y'])

Embeddings_X = pd.DataFrame( index=lf.ca.CellID )
Embeddings_X = pd.concat( [
        pd.DataFrame(adata.obsm['X_umap'],index=adata.obs.index)[0] ,
        pd.DataFrame(adata.obsm['X_pca'],index=adata.obs.index)[0] ,
        dr_tsne['X'] ,
        dr_umap['X']
    ], sort=False, axis=1, join='outer' )
Embeddings_X.columns = ['1','2','3','4']

Embeddings_Y = pd.DataFrame( index=lf.ca.CellID )
Embeddings_Y = pd.concat( [
        pd.DataFrame(adata.obsm['X_umap'],index=adata.obs.index)[1] ,
        pd.DataFrame(adata.obsm['X_pca'],index=adata.obs.index)[1] ,
        dr_tsne['Y'] ,
        dr_umap['Y']
    ], sort=False, axis=1, join='outer' )
Embeddings_Y.columns = ['1','2','3','4']


# Metadata:
# 
# 

# In[111]:


### metadata
metaJson = {}

metaJson['embeddings'] = [
    {
        "id": -1,
        "name": f"Scanpy t-SNE (highly variable genes)"
    },
    {
        "id": 1,
        "name": f"Scanpy UMAP  (highly variable genes)"
    },
    {
        "id": 2,
        "name": "Scanpy PC1/PC2"
    },
    {
        "id": 3,
        "name": "SCENIC AUC t-SNE"
    },
    {
        "id": 4,
        "name": "SCENIC AUC UMAP"
    },
]

metaJson["clusterings"] = [{
            "id": 0,
            "group": "Scanpy",
            "name": "Scanpy louvain default resolution",
            "clusters": [],
        }]

metaJson["metrics"] = [
        {
            "name": "nUMI"
        }, {
            "name": "nGene"
        }]

metaJson["annotations"] = [
    {
        "name": "Louvain_clusters_Scanpy",
        "values": list(set( adata.obs['louvain'].astype(np.str) ))
    },
    #{
    #    "name": "Genotype",
    #    "values": list(set(adata.obs['Genotype'].values))
    #},
    #{
    #    "name": "Timepoint",
    #    "values": list(set(adata.obs['Timepoint'].values))
    #},
    #{
    #    "name": "Sample",
    #    "values": list(set(adata.obs['Sample'].values))
    #}
]

# SCENIC regulon thresholds:
metaJson["regulonThresholds"] = rt

for i in range(max(set([int(x) for x in adata.obs['louvain']])) + 1):
    clustDict = {}
    clustDict['id'] = i
    clustDict['description'] = f'Unannotated Cluster {i + 1}'
    metaJson['clusterings'][0]['clusters'].append(clustDict)
    
clusterings = pd.DataFrame()
clusterings["0"] = adata.obs['louvain'].values.astype(np.int64)


# Assemble loom file row and column attributes

# In[112]:


def dfToNamedMatrix(df):
    arr_ip = [tuple(i) for i in df.values]
    dtyp = np.dtype(list(zip(df.dtypes.index, df.dtypes)))
    arr = np.array(arr_ip, dtype=dtyp)
    return arr


# In[113]:


adata


# In[114]:


# add annotation form cellID table
cellID = pd.read_csv('/diskmnt/Projects/Users/ritalu/20220331_HTAN_BRCA/02_pySCENIC/object/cellInfo_snRNA_0418.tsv', sep='\t')
cellID


# In[115]:


adata.obs['celltype'] = adata.obs.index.map(cellID.set_index('Barcode')['cell_type_specific'])


# In[116]:


adata.obs['Clinical_subtype'] = adata.obs.index.map(cellID.set_index('Barcode')['Clinical_subtype'])


# In[117]:


adata.obs['platform'] = adata.obs.index.map(cellID.set_index('Barcode')['Chemistry'])


# In[118]:


adata.obs['Piece_ID'] = adata.obs.index.map(cellID.set_index('Barcode')['Piece_ID'])


# In[119]:


col_attrs = {
    "CellID": np.array(adata.obs.index),
    "nUMI": np.array(adata.obs['nUMI'].values),
    "nGene": np.array(adata.obs['nGene'].values),
    "Louvain_clusters_Scanpy": np.array( adata.obs['louvain'].values ),
    "Celltype": np.array(adata.obs['celltype'].values),
    "Clinical_subtype": np.array(adata.obs['Clinical_subtype'].values),
    "PieceID": np.array(adata.obs['Piece_ID'].values),
    "Platform": np.array(adata.obs['platform'].values),
    "Embedding": dfToNamedMatrix(tsneDF),
    "Embeddings_X": dfToNamedMatrix(Embeddings_X),
    "Embeddings_Y": dfToNamedMatrix(Embeddings_Y),
    "RegulonsAUC": dfToNamedMatrix(auc_mtx),
    "Clusterings": dfToNamedMatrix(clusterings),
    "ClusterID": np.array(adata.obs['louvain'].values)
}


# In[120]:


lf


# In[121]:


row_attrs = {
    "Gene": lf.ra.Gene,
    "Regulons": regulons,
}


# In[122]:




attrs = {
    "title": "sampleTitle",
    "MetaData": json.dumps(metaJson),
    "Genome": 'hg38',
    "SCopeTreeL1": "",
    "SCopeTreeL2": "",
    "SCopeTreeL3": ""
}

# compress the metadata field:
attrs['MetaData'] = base64.b64encode(zlib.compress(json.dumps(metaJson).encode('ascii'))).decode('ascii')


# Create a new loom file, copying the expression matrix from the open loom connection:

# In[124]:


lp.create(
    filename = f_final_loom ,
    layers=lf[:,:],
    row_attrs=row_attrs, 
    col_attrs=col_attrs, 
    file_attrs=attrs
)
lf.close() # close original pyscenic loom file


# In[125]:


##saveadata
adata.write('anndata_integrated.h5ad' )


# In[126]:


adata


# In[ ]:




