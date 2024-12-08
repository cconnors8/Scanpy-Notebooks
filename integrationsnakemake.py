import pandas as pd
import scanpy as sc
import scanpy.external as sce
import matplotlib.pyplot as plt
import seaborn as sns

import os
import sys
import argparse

parser = argparse.ArgumentParser(description="Process single-cell RNA-seq data")
parser.add_argument("--base_dir", required=True, help="Base directory containing sample subdirectories")
parser.add_argument("--output_name", required=True, help="Output file name prefix")
parser.add_argument("--resolution", type=float, required=True, help="Resolution for clustering")

args = parser.parse_args()

base_dir = args.base_dir
name = args.output_name
res = args.resolution

# Read sample directories dynamically
sample_dirs = [os.path.join(base_dir, d) for d in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, d))]
sample_list = [os.path.basename(d) for d in sample_dirs]




#Message logging and figure settings
sc.settings.verbosity = 3            
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')

adatas = []
for sample_dir, sample_name in zip(sample_dirs, sample_list):
    adata2 = sc.read_10x_mtx(sample_dir, var_names='gene_symbols', cache=True)
    adata2.var_names_make_unique()
    adata2.obs['sample'] = sample_name
    adatas.append(adata2)

# Merge all samples into one AnnData object
adata = adatas[0].concatenate(adatas[1:], batch_key="sample")

results_file = "merged_" + name + '.h5ad'  # the file that will store the unintegrated analysis results



#Filters out cells with gene counts less than 200, over 150,000 and genes in less than 3 cells
#Can be adjusted based on sample, but these are well working baselines in our testing
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.filter_cells(adata, max_counts = 150000)

adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

#Filters out cells with less than 600 genes with 1 count in a cell and more than 20% mitochondrial counts
initial_cell_count = adata.n_obs

adata = adata[adata.obs.n_genes_by_counts > 600, :]

filtered_cell_count = adata.n_obs

adata = adata[adata.obs.pct_counts_mt < 20, :]

filtered_cell_count2 = adata.n_obs

cells_filtered_out = initial_cell_count - filtered_cell_count
cells_filtered_out2 = filtered_cell_count - filtered_cell_count2

print(f"Initial number of cells: {initial_cell_count}")
print(f"Number of cells after counts filtering: {filtered_cell_count}")
print(f"Number of cells filtered out in counts filtering: {cells_filtered_out}")
print(f"Number of cells after mt filtering: {filtered_cell_count2}")
print(f"Number of cells filtered out in mt filtering: {cells_filtered_out2}")


#runs scrublet to identify doublets in the data
sc.external.pp.scrublet(adata, batch_key = 'sample')

#filters out doublets based on scrublets reccomendation
adata = (adata[adata.obs['predicted_doublet'] == False])
#preserves the counts layer for scvi. not needed for any other tools.
adata.layers["counts"] = adata.X.copy()

#plots the genes by counts, total counts, and mitochondrial counts violin plots
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)

#normalizes and logarithmizes the data
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

#preserves the normalized and logarithmized raw gene expression
adata.raw = adata


#finds and plots highly variable genes
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5, batch_key = 'sample')
#benchmark adata pre-integration utilized for scib as a reference, 
benchmark = adata.copy()
sc.pl.highly_variable_genes(adata)

#filters the adata to include only the highly variable genes
adata = adata[:, adata.var.highly_variable]

#regresses out the total counts and mitochondrial counts before doing principal component analysis
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])

#Scales each gene to unit variance and clips values exceeding standard deviation 10.
sc.pp.scale(adata, max_value=10)

#performs principal coponent analysis and plots the variance ratio
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata, log=True)

#preserves a copy of the adata object
copydata = adata.copy()


#calls the standard dimensionality reduction and clustering steps on the adata object 
# neighbors, UMAP, and leiden (clustering)
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=30)
sc.tl.umap(adata)
#clustering resolution should be adjusted to ensure that it matches the elbow plot produced above
#cut off should be around the "elbow" or corner of the plot where intra-PC variance becomes minimal
res = .25
sc.tl.leiden(adata, resolution= res, key_added = "leiden")

#plots the umap based on both sample and cluster
sc.pl.umap(adata, color = ['leiden'], hspace = 1)
sc.pl.umap(adata, color = ['sample'], hspace = 1)

#plots the umap with just one sample colored in at a time -> useful for sample by sample location visualization
for sample in sample_list:
    sc.pl.umap(adata, color = 'sample', groups = [sample])
#writes the adata to a results file in the working directory
adata.write(results_file)

#loads in the saved adata copy
bbknnadata = copydata.copy()

#runs bbknn on the adata
sc.external.pp.bbknn(bbknnadata, batch_key='sample', n_pcs=30) 

# run umaps and clusters and on the integrated space
sc.tl.umap(bbknnadata)
#adjust resolution to match the cell above
sc.tl.leiden(bbknnadata, resolution= res, key_added = "leiden")

#plots the umap colored by both cluster and sample
sc.pl.umap(bbknnadata, color="leiden", title="BBKNN Corrected umap", show=False)
sc.pl.umap(bbknnadata, color="sample", title="BBKNN Corrected umap", show=False)
plt.show()

#plots the umap with just one sample colored in at a time -> useful for sample by sample location visualization
for sample in sample_list:
    sc.pl.umap(bbknnadata, color = 'sample', groups = [sample])

plt.show()


#loads in a copy of the adata
scadata = copydata.copy()
#runs scanorama on the data
sce.pp.scanorama_integrate(scadata, 'batch')
#calls neighbors, umap, and leiden on the data
sc.pp.neighbors(scadata, n_neighbors= 10, n_pcs =30, use_rep = "X_scanorama")
sc.tl.umap(scadata)
sc.tl.leiden(scadata, resolution= res, key_added = "leiden")

#plots the umap colored by both cluster and sample
sc.pl.umap(scadata, color="leiden", title="Scanorama umap", show=False)
sc.pl.umap(scadata, color="sample", title="Scanorama umap", show=False)

plt.show()
#plots the umap with just one sample colored in at a time -> useful for sample by sample location visualization
for sample in sample_list:
    sc.pl.umap(scadata, color = 'sample', groups = [sample])

plt.show()

#loads in a copy of the adata
harmadata = copydata.copy()

#integrates the data with harmony
sce.pp.harmony_integrate(harmadata, 'batch')

#replaces the X_pca category with the harmony integrated X_pca_harmony
harmadata.obsm['X_pca'] = harmadata.obsm['X_pca_harmony']

#calls neighbors, umap, leiden
sc.pp.neighbors(harmadata, n_neighbors=10, n_pcs=30)
sc.tl.umap(harmadata)
sc.tl.leiden(harmadata, resolution= res, key_added = "leiden")

#plots the umap colored by both the clusters and the sample
sc.pl.umap(harmadata, color="leiden", title="Harmony umap", show=False)
sc.pl.umap(harmadata, color="sample", title="Harmony umap", show=False)
plt.show()
#plots the umap with just one sample colored in at a time -> useful for sample by sample location visualization
for sample in sample_list:
    sc.pl.umap(harmadata, color = 'sample', groups = [sample])

plt.show()


import scvi

#loads in a copy of the adata object
scvidata = copydata.copy()

#sets up the SCVI model
scvi.model.SCVI.setup_anndata(scvidata, layer="counts", batch_key="sample")

#runs the model on the adata
model = scvi.model.SCVI(scvidata, n_layers=2, n_latent=30, gene_likelihood="nb")
model.train()

#adds the correction into the adata object
SCVI_LATENT_KEY = "X_scVI"
scvidata.obsm[SCVI_LATENT_KEY] = model.get_latent_representation()

#runs neighbors using the correction as opposed to the standard PCA output
sc.pp.neighbors(scvidata, use_rep=SCVI_LATENT_KEY)

#calls leiden clustering
sc.tl.leiden(scvidata, resolution = res, key_added = "leiden")


#runs the mde (gpu accelerated UMAP) on the adata object

SCVI_MDE_KEY = "X_scVI_MDE"
scvidata.obsm[SCVI_MDE_KEY] = scvi.model.utils.mde(scvidata.obsm[SCVI_LATENT_KEY])

sc.set_figure_params(figsize = (8,6))
#plots the embedding by sample and cluster
sc.pl.embedding(
    scvidata,
    basis=SCVI_MDE_KEY,
    color=["batch", "leiden"],
    frameon=False,
    ncols=1,
)

sc.tl.umap(scvidata)
sc.tl.leiden(scvidata, resolution = res, key_added = "leiden")
sc.pl.umap(scvidata, color="leiden", title="SCVI umap", show=False)
sc.pl.umap(scvidata, color="sample", title="SCVI umap", show=False)

plt.show()
#plots the umap with just one sample colored in at a time -> useful for sample by sample location visualization
for sample in sample_list:
    sc.pl.umap(scvidata, color = 'sample', groups = [sample])

plt.show()

import scib
metrics_scvi = scib.metrics.metrics(
    benchmark, scvidata, batch_key = 'sample', embed="X_scVI", label_key = "sample", hvg_score_ = True, pcr_ = True, graph_conn_ = True 
)
metrics_scanorama = scib.metrics.metrics(
    benchmark, scadata, batch_key = 'sample', embed="X_scanorama", label_key = "sample", hvg_score_ = True, pcr_ = True, graph_conn_ = True
)
metrics_bbknn = scib.metrics.metrics(benchmark, bbknnadata, batch_key = 'sample', label_key = "sample", hvg_score_ = True, pcr_ = True, graph_conn_ = True)
metrics_harmony = scib.metrics.metrics(benchmark, harmadata, batch_key = 'sample', label_key = "sample", hvg_score_ = True, pcr_ = True, graph_conn_ = True)
metrics = pd.concat(
    [metrics_scvi, metrics_scanorama, metrics_bbknn, metrics_harmony],
    axis="columns",
)
# Set methods as column names
metrics = metrics.set_axis(
    ["SCVI", "scanorama", "BBKNN", "harmony"], axis="columns"
)
# Select only the fast metrics
metrics = metrics.loc[
    [
        "PCR_batch",
        "graph_conn",
        "hvg_overlap",
    ],
    :,
]
# Transpose so that metrics are columns and methods are rows

metrics = metrics.T

metrics
metrics.style.background_gradient(cmap="Blues")
metrics_scaled = (metrics - metrics.min()) / (metrics.max() - metrics.min())
metrics_scaled.style.background_gradient(cmap="Blues")
metrics_scaled["Batch"] = metrics_scaled[
    ["hvg_overlap", "PCR_batch", "graph_conn"]
].mean(axis=1)

metrics
# fig, ax = plt.subplots()
# ax.set_xlim(0, 1)
# ax.set_ylim(0, 1)
# metrics_scaled.plot.scatter(
#     x="Batch",
#     y="Bio",
#     c=range(len(metrics_scaled)),
#     ax=ax,
# )

# for k, v in metrics_scaled[["Batch", "Bio"]].iterrows():
#     ax.annotate(
#         k,
#         v,
#         xytext=(6, -3),
#         textcoords="offset points",
#         family="sans-serif",
#         fontsize=12,
#     )
# metrics_scaled["Overall"] = 0.4 * metrics_scaled["Batch"] + 0.6 * metrics_scaled["Bio"]
# metrics_scaled.style.background_gradient(cmap="Blues")
# metrics_scaled.plot.bar(y="Overall")

bbknnadata.write_h5ad('BBKNN' + name + '.h5ad')
harmadata.write_h5ad('harmony' + name + '.h5ad')
scadata.write_h5ad('scanorama' + name + '.h5ad')
scvidata.write_h5ad('scvi' + name + '.h5ad')