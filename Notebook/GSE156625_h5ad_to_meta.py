#!/usr/bin/env python
# coding: utf-8

# In[9]:


import scanpy as sc
import pandas as pd


# In[10]:


adata = sc.read_h5ad("/parallel_scratch/mp01950/raw_data/liver/saved_RDS/metadata/GSE156625_HCCscanpyobj.h5ad")


# In[11]:


columns_to_subset = ['patient_id', 'patient_tumorsection','NormalvsTumor','patientno','PNC','PIC','ViralvsNonViral']


# In[16]:


subset_df = adata.obs[columns_to_subset].copy()


# In[17]:


subset_df


# In[18]:


subset_df.to_csv('subset_data.csv', index=True)

