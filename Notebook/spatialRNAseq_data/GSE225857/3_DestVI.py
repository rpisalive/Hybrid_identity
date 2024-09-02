#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#Modules loading


# In[1]:


import tempfile
import destvi_utils
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from adjustText import adjust_text
from scipy.stats import ks_2samp, ttest_ind
from statsmodels.stats.multitest import multipletests
from matplotlib.patches import Patch
import numpy as np
import scanpy as sc
import umap
import scvi
import seaborn as sns
import torch
import pandas as pd
import jax
import hotspot
from scvi.model import CondSCVI, DestVI
from torch.distributions import Gamma
import gseapy
import base64
from io import BytesIO
get_ipython().run_line_magic('matplotlib', 'inline')


# In[ ]:


#Check GPU availability & scvi-tools version


# In[2]:


print(f'Jax backend: {jax.default_backend()}')
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"Using device: {device}")
scvi.settings.seed = 0
print("Last run with scvi-tools version:", scvi.__version__)


# In[ ]:


#Path & plot specification setup


# In[3]:


slice_name = "C1"
spa_path = "/parallel_scratch/mp01950/spatial_study/"
acc_no = "GSE225857/"
slice_loc = "spatial/"+slice_name+"/"
full_slice_path = spa_path+acc_no+slice_loc
save_dir = full_slice_path+"/DestVI_output/"
sc.set_figure_params(figsize=(6, 6), frameon=False)
sns.set_theme()
torch.set_float32_matmul_precision("high")
get_ipython().run_line_magic('config', 'InlineBackend.print_figure_kwargs={"facecolor": "w"}')
get_ipython().run_line_magic('config', 'InlineBackend.figure_format="retina"')


# In[ ]:


#Single cell reference loading 


# In[4]:


sc_ref = sc.read_h5ad(spa_path+acc_no+"scRNA/sc_ref.h5ad.gz")


# In[ ]:


#Genes filtering on scRNA


# In[5]:


# let us filter some genes
G = 3000 #This determies how much genes to include to infer the cell type proportion, could lead to huge differences if less/more genes are to be included
sc.pp.filter_genes(sc_ref, min_counts=10)

sc_ref.layers["counts"] = sc_ref.X.copy()

sc.pp.highly_variable_genes(
    sc_ref, n_top_genes=G, subset=True, layer="counts", flavor="seurat_v3"
)

sc.pp.normalize_total(sc_ref, target_sum=10e4)
sc.pp.log1p(sc_ref)
sc_ref.raw = sc_ref
sc.pp.scale(sc_ref, max_value=10)
sc.tl.pca(sc_ref, svd_solver='arpack')
sc.pl.pca_variance_ratio(sc_ref, log=True, n_pcs = 50)
sc.pp.neighbors(sc_ref, n_pcs = 45)
sc.tl.umap(sc_ref)
sc.pl.umap(sc_ref, color = "general_cell_type")


# In[6]:


if "CD3D" in sc_ref.var_names:
    print("CD3D is present in var_names.")
else:
    print("CD3D is not present in var_names.")


# In[ ]:


#Spatial data loading & layer setting


# In[7]:


st_adata = sc.read_h5ad(full_slice_path+slice_name+".h5ad.gz")
st_adata.layers["counts"] = st_adata.X.copy()


# In[ ]:


#Intersecting genes in both scRNA and spatial RNA


# In[8]:


# filter genes to be the same on the spatial data
intersect = np.intersect1d(sc_ref.var_names, st_adata.var_names)
st_adata = st_adata[:, intersect].copy()
sc_ref = sc_ref[:, intersect].copy()
G = len(intersect)
sc.pl.embedding(st_adata, basis="spatial", s=80)


# In[ ]:


#Fitting scRNA-seq model


# In[9]:


CondSCVI.setup_anndata(sc_ref, layer="counts", labels_key="general_cell_type")


# In[10]:


train = False
if train:
    # train the conditional VAE
    sc_model = CondSCVI(sc_ref, weight_obs=True)
    sc_model.train(max_epochs=250)
    sc_model.history["elbo_train"].plot()
    sc_model.save(save_dir+"sc_model", overwrite=True)
else:
    sc_model = CondSCVI.load(save_dir+"sc_model/", sc_ref)


# In[11]:


sc_ref.obsm["X_CondSCVI"] = sc_model.get_latent_representation()


# In[ ]:


#Deconvolution


# In[12]:


sc.pp.normalize_total(st_adata, target_sum=10e4)
sc.pp.log1p(st_adata)
st_adata.raw = st_adata


# In[13]:


# get dataset ready
DestVI.setup_anndata(st_adata, layer="counts")


# In[14]:


# add here number of cell type
train_st = False
if train_st:
    st_model = DestVI.from_rna_model(st_adata, sc_model, vamp_prior_p=100)
    st_model.train(max_epochs=2500, plan_kwargs={"lr":0.005})
    st_model.history["elbo_train"].plot()
    st_model.save(save_dir+"st_model", overwrite=True)
else:
    st_model = DestVI.load(save_dir+"st_model", st_adata)


# In[15]:


st_adata.obsm["proportions"] = st_model.get_proportions()


# In[16]:


ct_list = ["B", "T", "Endothelial_cell", "Myeloid", "NK", "Plasma", "Tumor", "Fibroblast"]
for ct in ct_list:
    data = st_adata.obsm["proportions"][ct].values
    st_adata.obs[ct] = np.clip(data, 0, np.quantile(data, 0.99))


# In[17]:


sc.pl.embedding(st_adata, basis="spatial", color=ct_list, cmap="magma", s=100)


# In[18]:


ct_thresholds = destvi_utils.automatic_proportion_threshold(
    st_adata, ct_list=ct_list, kind_threshold="secondary"
)


# In[ ]:


#Intra cell type information


# In[19]:


# more globally, the values of the gamma are all summarized in this dictionary of data frames
for ct, g in st_model.get_gamma().items():
    st_adata.obsm[f"{ct}_gamma"] = g


# In[20]:


destvi_utils.explore_gamma_space(st_model, sc_model, ct_list=ct_list, ct_thresholds=ct_thresholds)


# In[ ]:


#Expression of gene signatures of other cells in T cells


# In[21]:


#Tumor cell
plt.figure(figsize=(8, 8))

ct_name = "T"
gene_name = ["EPCAM","SOX9"]

# we must filter spots with low abundance (consult the paper for an automatic procedure)
indices = np.where(st_adata.obsm["proportions"][ct_name].values > ct_thresholds["T"])[0]

# impute genes and combine them
specific_expression = np.sum(st_model.get_scale_for_ct(ct_name, indices=indices)[gene_name], 1)
specific_expression = np.log(1 + 1e4 * specific_expression)

# plot (i) background (ii) g
plt.scatter(st_adata.obsm["spatial"][:, 0], st_adata.obsm["spatial"][:, 1], alpha=0.05)
plt.scatter(
    st_adata.obsm["spatial"][indices][:, 0],
    st_adata.obsm["spatial"][indices][:, 1],
    c=specific_expression,
    s=10,
    cmap="Reds",
)
plt.colorbar()
plt.title(f"Imputation of {gene_name} in {ct_name}")
plt.show()


# In[22]:


#Fibroblasts
plt.figure(figsize=(8, 8))

ct_name = "T"
gene_name = ["COL1A1", "COL1A2"]

# we must filter spots with low abundance (consult the paper for an automatic procedure)
indices = np.where(st_adata.obsm["proportions"][ct_name].values > ct_thresholds["T"])[0]

# impute genes and combine them
specific_expression = np.sum(st_model.get_scale_for_ct(ct_name, indices=indices)[gene_name], 1)
specific_expression = np.log(1 + 1e4 * specific_expression)

# plot (i) background (ii) g
plt.scatter(st_adata.obsm["spatial"][:, 0], st_adata.obsm["spatial"][:, 1], alpha=0.05)
plt.scatter(
    st_adata.obsm["spatial"][indices][:, 0],
    st_adata.obsm["spatial"][indices][:, 1],
    c=specific_expression,
    s=10,
    cmap="Reds",
)
plt.colorbar()
plt.title(f"Imputation of {gene_name} in {ct_name}")
plt.show()


# In[23]:


#Endothelial cells
plt.figure(figsize=(8, 8))

ct_name = "T"
gene_name = ["CD34","PECAM1"]

# we must filter spots with low abundance (consult the paper for an automatic procedure)
indices = np.where(st_adata.obsm["proportions"][ct_name].values > ct_thresholds["T"])[0]

# impute genes and combine them
specific_expression = np.sum(st_model.get_scale_for_ct(ct_name, indices=indices)[gene_name], 1)
specific_expression = np.log(1 + 1e4 * specific_expression)

# plot (i) background (ii) g
plt.scatter(st_adata.obsm["spatial"][:, 0], st_adata.obsm["spatial"][:, 1], alpha=0.05)
plt.scatter(
    st_adata.obsm["spatial"][indices][:, 0],
    st_adata.obsm["spatial"][indices][:, 1],
    c=specific_expression,
    s=10,
    cmap="Reds",
)
plt.colorbar()
plt.title(f"Imputation of {gene_name} in {ct_name}")
plt.show()


# In[24]:


#Monocytes/macrophages
plt.figure(figsize=(8, 8))

ct_name = "T"
gene_name = ["CD14"]

# we must filter spots with low abundance (consult the paper for an automatic procedure)
indices = np.where(st_adata.obsm["proportions"][ct_name].values > ct_thresholds["T"])[0]

# impute genes and combine them
specific_expression = np.sum(st_model.get_scale_for_ct(ct_name, indices=indices)[gene_name], 1)
specific_expression = np.log(1 + 1e4 * specific_expression)

# plot (i) background (ii) g
plt.scatter(st_adata.obsm["spatial"][:, 0], st_adata.obsm["spatial"][:, 1], alpha=0.05)
plt.scatter(
    st_adata.obsm["spatial"][indices][:, 0],
    st_adata.obsm["spatial"][indices][:, 1],
    c=specific_expression,
    s=10,
    cmap="Reds",
)
plt.colorbar()
plt.title(f"Imputation of {gene_name} in {ct_name}")
plt.show()


# In[25]:


#B/Plasma
plt.figure(figsize=(8, 8))

ct_name = "T"
gene_name = ["CD19"]

# we must filter spots with low abundance (consult the paper for an automatic procedure)
indices = np.where(st_adata.obsm["proportions"][ct_name].values > ct_thresholds["T"])[0]

# impute genes and combine them
specific_expression = np.sum(st_model.get_scale_for_ct(ct_name, indices=indices)[gene_name], 1)
specific_expression = np.log(1 + 1e4 * specific_expression)

# plot (i) background (ii) g
plt.scatter(st_adata.obsm["spatial"][:, 0], st_adata.obsm["spatial"][:, 1], alpha=0.05)
plt.scatter(
    st_adata.obsm["spatial"][indices][:, 0],
    st_adata.obsm["spatial"][indices][:, 1],
    c=specific_expression,
    s=10,
    cmap="Reds",
)
plt.colorbar()
plt.title(f"Imputation of {gene_name} in {ct_name}")
plt.show()


# In[26]:


ct = "T"
imputation = st_model.get_scale_for_ct(ct)
color = np.log(1 + 1e4 * imputation["CD19"].values)
threshold = 0.3

mask = color > threshold
#mask2 = np.logical_and(st_adata.obs["_indices"].isin(np.where(st_adata.obsm["proportions"][ct].values > 0.02)[0]),color < threshold,).values

_ = destvi_utils.de_genes(st_model, mask=mask, threshold=ct_thresholds[ct], ct=ct, key="B_signature")

display(st_adata.uns["B_signature"]["de_results"].head(10))

destvi_utils.plot_de_genes(
    st_adata,
    interesting_genes=["CD19"],
    key="B_signature",
)


# In[ ]:


#Locate the indices of the spots of interest


# In[27]:


# Find the indices where the values are larger than 0.3
indices = np.where(color > 0.3)[0] #Modify the criteria to highlight spots of interest
indices


# In[ ]:


#Subset the cell type proportion for spots of interest


# In[28]:


indices_list = []
for i in range(len(indices)):
    globals()["index_"+str(i)] = imputation.index[indices[i]]
    indices_list.append(globals()["index_"+str(i)])
spots_of_interest = st_adata.obs[st_adata.obs_names.isin(indices_list)]
spots_of_interest = spots_of_interest.drop(columns=["array_row", "array_col", "_indices"])
spots_of_interest["ct_sum"] = spots_of_interest.sum(axis=1)
spots_of_interest["others"] = 1-spots_of_interest["ct_sum"]


# In[29]:


spots_of_interest


# In[ ]:


#Create data frame for gene of interest


# In[55]:


gene_of_interest = "CD19" #One gene supported at the moment
for ct in ct_list:
    globals()["imputation_"+ct] = st_model.get_scale_for_ct(ct)
    globals()["imputation_"+ct] = globals()["imputation_"+ct][globals()["imputation_"+ct].index.isin(indices_list)]
    globals()["imputation_"+ct] = globals()["imputation_"+ct][[gene_of_interest]]
    globals()["imputation_"+ct] = globals()["imputation_"+ct].rename(columns={gene_of_interest: ct})
globals()[gene_of_interest] = pd.concat([imputation_T, imputation_B, imputation_Endothelial_cell, imputation_NK, imputation_Myeloid, imputation_Plasma, imputation_Tumor, imputation_Fibroblast], axis=1)
globals()[gene_of_interest]["expression_sum"] = globals()[gene_of_interest].sum(axis=1)
globals()[gene_of_interest]["others"] = 0


# In[ ]:


#Plotting


# In[57]:


# Define colors and cell types
cell_types = ['B', 'T', 'Endothelial_cell', 'Myeloid', 'NK', 'Plasma', 'Tumor', 'Fibroblast', 'others']
colors = ['#4daf4a', '#377eb8', '#ff7f00', '#984ea3', '#e41a1c', '#ffff33', '#a65628', '#f781bf', '#07090d']

# Normalize CD19 expression values to proportions
CD19_proportions = CD19[cell_types].div(CD19['expression_sum'], axis=0)

# Set up the plot
fig, ax = plt.subplots(figsize=(6, 6))

# Bar width and positions
bar_width = 0.2
indices = np.arange(len(spots_of_interest.index))

# Plot stacked bars for cell type proportions with solid border
bars1 = []  # To store bar handles for legend
for i, cell_type in enumerate(cell_types):
    bar = ax.bar(
        indices - bar_width / 2,
        spots_of_interest[cell_type],
        bar_width,
        bottom=spots_of_interest[cell_types[:i]].sum(axis=1) if i > 0 else 0,
        color=colors[i],
        edgecolor='black',  # Solid border
        linewidth=1.5
    )
    bars1.append(bar[0])

# Plot stacked bars for CD19 expression proportions with dashed border
bars2 = []  # To store bar handles for legend
for i, cell_type in enumerate(cell_types):
    bar = ax.bar(
        indices + bar_width / 2,
        CD19_proportions[cell_type],
        bar_width,
        bottom=CD19_proportions[cell_types[:i]].sum(axis=1) if i > 0 else 0,
        color=colors[i],
        edgecolor='black',  # Dashed border
        linestyle='--',
        linewidth=1.5
    )
    bars2.append(bar[0])

# Labeling
ax.set_xlabel('Spots of Interest')
ax.set_ylabel('Proportion')
ax.set_title('Cell Type Proportion and CD19 Expression Proportion')
ax.set_xticks(indices)
ax.set_xticklabels(spots_of_interest.index, rotation=45)

# Create a legend for the cell types
cell_type_legend = ax.legend(
    handles=bars1,
    labels=cell_types,
    loc='upper left',
    bbox_to_anchor=(1.05, 1),  # Position the legend to the right of the plot
    title="Cell Types"
)

# Create custom legend entries for bar styles
proportion_patch = Patch(facecolor='none', edgecolor='black', linestyle='-', label='Cell Type Proportion')
expression_patch = Patch(facecolor='none', edgecolor='black', linestyle='--', label='CD19 Expression')

# Add custom legend for the bar types, positioning it directly below the cell type legend
bar_type_legend = ax.legend(
    handles=[proportion_patch, expression_patch],
    loc='upper left',
    bbox_to_anchor=(1.05, 0.3),  # Position this below the Cell Types legend
    title="Bar Types"
)

# Add the cell type legend back to avoid it being overwritten by the second legend
ax.add_artist(cell_type_legend)

# Adjust layout to accommodate both legends
plt.tight_layout()
plt.subplots_adjust(right=0.8)  # Adjust right margin to make room for legends

plt.show()


# CD45 - immune vs. non-immune
# tummor cells - EPCAM & SOX9
# Fibroblast - COL1A1 & COL1A2
# Endothelial cells PECAM1 & CD34
# T - CD3 & pAbO
# NK - CD56 & pAbO
# B/plasma - CD19 & pAbO
# Monocytes/macrophages - CD14 & pAbO
# DCs - HLA & DRA
# Mast cells - TPSAB1

# In[ ]:





# In[ ]:





# In[ ]:




