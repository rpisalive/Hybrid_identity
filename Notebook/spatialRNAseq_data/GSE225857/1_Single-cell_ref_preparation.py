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

# In[2]:


slice_name = "C1"
spa_path = "/parallel_scratch/mp01950/spatial_study/"
acc_no = "GSE225857/"
slice_loc = "spatial/"+slice_name+"/"
full_slice_path = spa_path+acc_no+slice_loc


# #sc data loading and concatenation

# In[3]:


sc_immune_ref = sc.read_text(spa_path+acc_no+"scRNA/GSM7058754_immune_counts.txt.gz", delimiter='\t', first_column_names=True)
sc_nonim_ref = sc.read_text(spa_path+acc_no+"scRNA/GSM7058755_non_immune_counts.txt.gz", delimiter='\t', first_column_names=True)
#Transpose cell_id & features
sc_immune_ref = sc_immune_ref.T
sc_nonim_ref = sc_nonim_ref.T
sc_immune_ref.var_names_make_unique()
sc_nonim_ref.var_names_make_unique()
#Concatenate two anndata objects
sc_ref = sc_immune_ref.concatenate(sc_nonim_ref,batch_key = "Sample", join = "outer")
#Rename to remove the last two characters generated during concatenation
rename = []
for i in range(len(sc_ref)):
    rename.append(sc_ref.obs.index[i][0:-2])
sc_ref.obs.index = rename
# Renaming the "Sample" column to Dataset_no
sc_ref.obs = sc_ref.obs.rename(columns={'Sample': 'Dataset_no'})
# Replace '.' with '-' in the index names
sc_ref.obs.index = sc_ref.obs.index.str.replace('.', '-', regex=False)


# Metadata transfer

# In[4]:


immune_meta = pd.read_csv(spa_path+acc_no+"scRNA/GSM7058754_immune_meta.txt.gz", sep='\t', compression='gzip', index_col = 0)
nonim_meta = pd.read_csv(spa_path+acc_no+"scRNA/GSM7058755_non_immune_meta.txt.gz", sep='\t', compression='gzip', index_col = 0)
#Modify metadata
immune_meta = immune_meta.drop(['nCount_antibody', 'nFeature_antibody'], axis=1)
nonim_meta = nonim_meta.drop(['integrated_snn_res.0.1'], axis=1)
nonim_meta["samples"] = "NA"
#Merging two metadata
combined_meta = pd.concat([immune_meta, nonim_meta], ignore_index=False)
# Ensure both DataFrames have the same index
sc_ref.obs = sc_ref.obs.join(combined_meta)


# General cell type annotation

# In[5]:


# Function to annotate general cell types based on "cluster"
def annotate_general_cell_type(cluster_name):
    first_char = cluster_name[0]
    first_two_chars = cluster_name[:2]
    fifth_char = cluster_name[4] if len(cluster_name) >= 5 else None

    if first_two_chars == "Tu":
        return "Tumor"
    elif first_char == "T":
        return "T"
    elif first_char == "N":
        return "NK"
    elif first_char == "B":
        if fifth_char == "p":
            return "Plasma"
        elif fifth_char in ["B", "G"]:
            return "B"
    elif first_char == "E":
        return "Endothelial_cell"
    elif first_char == "F":
        return "Fibroblast"
    else:
        return "Myeloid"

sc_ref.obs['general_cell_type'] = sc_ref.obs['cluster'].apply(annotate_general_cell_type)


# Anndata object saving

# In[8]:


sc_ref.write_h5ad(spa_path+acc_no+"scRNA/sc_ref.h5ad.gz",compression="gzip")


# In[ ]:




