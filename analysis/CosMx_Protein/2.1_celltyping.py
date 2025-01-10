#Import
import numpy as np
import pandas as pd
from scipy.io import mmread

import matplotlib.pyplot as plt
plt.rcParams["figure.figsize"] = (6, 4)

from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay

import anndata as ad
import scanpy as sc
import maxfuse as mf

import seaborn as sns

# Create Anndata for RNA SINGLE CELL
count_matrix = pd.read_csv('/home/mmoro/SPATIAL/Mackensy_analysis/CosMx_RNA/1_Annotation/Objects/counts_SC_cut.csv', index_col=0)  
metadata = pd.read_csv('/home/mmoro/SPATIAL/Mackensy_analysis/CosMx_RNA/1_Annotation/Objects/meta_SC_cut.csv', index_col=0)  
rna_adata = ad.AnnData(X=count_matrix.T)
rna_adata.obs = metadata

#Create Anndata for protein
count_protein = pd.read_csv('/home/mmoro/SPATIAL/Mackensy_analysis/CosMx_Prot/0_Curation/Objects/counts.csv', index_col=0)  
norm_protein = pd.read_csv('/home/mmoro/SPATIAL/Mackensy_analysis/CosMx_Prot/02_celltyping/counts.csv', index_col=0)  
protein_adata = ad.AnnData(X= norm_protein.T )

correspondence = pd.read_csv('/home/mmoro/SPATIAL/Mackensy_analysis/CosMx_Prot/protein_gene_conversion.csv')

rna_protein_correspondence = []

for i in range(correspondence.shape[0]):
    curr_protein_name, curr_rna_names = correspondence.iloc[i]
    if curr_protein_name not in protein_adata.var_names:
        continue
    if curr_rna_names.find('Ignore') != -1: # some correspondence ignored eg. protein isoform to one gene
        continue
    curr_rna_names = curr_rna_names.split('/') # eg. one protein to multiple genes
    for r in curr_rna_names:
        if r in rna_adata.var_names:
            rna_protein_correspondence.append([r, curr_protein_name])
            
rna_protein_correspondence = np.array(rna_protein_correspondence)

rna_shared = rna_adata[:, rna_protein_correspondence[:, 0]].copy()
protein_shared = protein_adata[:, rna_protein_correspondence[:, 1]].copy()

# Make sure no column is static
mask = (
    (rna_shared.X.std(axis=0) > 0.5) 
    & (protein_shared.X.std(axis=0) > 0.1)
)
rna_shared = rna_shared[:, mask].copy()
protein_shared = protein_shared[:, mask].copy()

# process rna_shared
sc.pp.normalize_total(rna_shared)
sc.pp.log1p(rna_shared)
sc.pp.scale(rna_shared,max_value=10)

# plot UMAP of rna cells based only on rna markers with protein correspondence
sc.pp.neighbors(rna_shared, n_neighbors=15)
sc.tl.umap(rna_shared)
sc.pl.umap(rna_shared, color= 'anot_maxfuse')

#plot UMAP of protein cells based on protein markers with rna correspondence
sc.pp.neighbors(protein_shared, n_neighbors=15)
sc.tl.umap(protein_shared)
sc.pl.umap(protein_shared)

#Observe RNA adata
sc.pp.normalize_total(rna_adata)
sc.pp.log1p(rna_adata)
sc.pp.scale(rna_adata)
sc.pp.neighbors(rna_adata, n_neighbors=15)
sc.tl.umap(rna_adata)
sc.pl.umap(rna_adata, color='anot_maxfuse')

#Observe protein adata
sc.pp.neighbors(protein_adata, n_neighbors=15)
sc.tl.umap(protein_adata)

# make sure no feature is static
rna_active = rna_adata.X
protein_active = protein_adata.X
rna_active = rna_active[:, rna_active.std(axis=0) > 1e-5] 
protein_active = protein_active[:, protein_active.std(axis=0) > 1e-5] 

labels1_array = rna_adata.obs['subset'].to_numpy()
labels2_array = np.resize(labels1_array, 77057)

rna_shared = rna_shared.X.copy()
protein_shared = protein_shared.X.copy()

# inspect shape of the four matrices
print(rna_active.shape)
print(protein_active.shape)
print(rna_shared.shape)
print(protein_shared.shape)

# call constructor for Fusor object
# which is the main object for running MaxFuse pipeline
fusor = mf.model.Fusor(
    shared_arr1=rna_shared,
    shared_arr2=protein_shared,
    active_arr1=rna_active,
    active_arr2=protein_active,
    labels1=None,
    labels2=None
)

fusor.split_into_batches(
    max_outward_size=8000,
    matching_ratio=4,
    metacell_size=2,
    verbose=True
)

# plot top singular values of avtive_arr1 on a random batch
fusor.plot_singular_values(
    target='active_arr1',
    n_components=None # can also explicitly specify the number of components
)

# plot top singular values of avtive_arr2 on a random batch
fusor.plot_singular_values(
    target='active_arr2',
    n_components=None
)

fusor.construct_graphs(
    n_neighbors1=15,
    n_neighbors2=15,
    svd_components1=50,
    svd_components2=20,
    resolution1=2,
    resolution2=2,
    # if two resolutions differ less than resolution_tol
    # then we do not distinguish between then
    resolution_tol=0.1,
    verbose=True
)

# plot top singular values of shared_arr1 on a random batch
fusor.plot_singular_values(
    target='shared_arr1',
    n_components=None
)

# plot top singular values of shared_arr2 on a random batch
fusor.plot_singular_values(
    target='shared_arr2',
    n_components=None
)

fusor.find_initial_pivots(
    wt1=0.3, wt2=0.3,
    svd_components1=25, svd_components2=20
)

# plot top canonical correlations in a random batch
fusor.plot_canonical_correlations(
    svd_components1=30,
    svd_components2=30,
    cca_components=30
)

fusor.refine_pivots(
    wt1=0.3, wt2=0.3,
    svd_components1=30, svd_components2=30,
    cca_components=20,
    n_iters=1,
    randomized_svd=False, 
    svd_runs=1,
    verbose=True
)

fusor.filter_bad_matches(target='pivot', filter_prop=0.3)
pivot_matching = fusor.get_matching(target='pivot')

fusor.propagate(
    svd_components1=30, 
    svd_components2=30, 
    wt1=0.7,
    wt2=0.7,
)

fusor.filter_bad_matches(
    target='propagated',
    filter_prop=0
)

full_matching = fusor.get_matching(order=(2, 1), target='full_data')

connections = pd.DataFrame(list(zip(full_matching[0], full_matching[1], full_matching[2])), 
             columns = ['mod1_indx', 'mod2_indx', 'score'])
# columns: cell idx in mod1, cell idx in mod2, and matching scores

connections.to_csv("connections.csv")

