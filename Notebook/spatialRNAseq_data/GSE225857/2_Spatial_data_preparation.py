#!/usr/bin/env python
# coding: utf-8

# In[1]:


import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import scvi
import seaborn as sns
import pandas as pd


# Path setup

# In[7]:


slice_name = "C2"
spa_path = "/parallel_scratch/mp01950/Hybrid_cell_spatial/"
acc_no = "GSE225857/"
slice_loc = "spatial/"+slice_name+"/"
full_slice_path = spa_path+acc_no+slice_loc
sub_acc_no = "GSM7058757_"


# Spatial data loading & modification

# In[8]:


#Loading spatial data
st_adata=sc.read_mtx(full_slice_path+sub_acc_no+slice_name+".matrix.mtx.gz")
st_adata_bc=pd.read_csv(full_slice_path+sub_acc_no+slice_name+".barcodes.tsv.gz", header=None)
st_adata_features=pd.read_csv(full_slice_path+sub_acc_no+slice_name+".features.tsv.gz",header=None, sep='\t')
st_position=pd.read_csv(full_slice_path+sub_acc_no+slice_name+'_tissue_positions_list.csv', header=None)
st_position=st_position[st_position.iloc[:, 1] == 1]
st_position = st_position.sort_values(by=st_position.columns[0])
st_adata = st_adata.T
st_adata.var_names = st_adata_features[1].tolist()
st_adata.var['gene_ids']= st_adata_features[0].tolist()
st_adata.var['feature_types']= st_adata_features[2].tolist()
st_adata.obs.index = st_adata_bc.iloc[:,0]
st_adata.var_names_make_unique()
st_adata.obs["array_row"]=st_position[2].tolist()
st_adata.obs["array_col"]=st_position[3].tolist()
st_position.index=st_position[0].tolist()
st_adata.obsm["spatial"]=st_position.iloc[:, [4,5]].values
st_adata.obs.index.name = None


# In[9]:


st_adata.write_h5ad(full_slice_path+slice_name+".h5ad.gz",compression="gzip")


# In[ ]:





# In[ ]:





# In[ ]:




