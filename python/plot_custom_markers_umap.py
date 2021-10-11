"""
Plotting umaps of a small list of custom genes
CRG 2020-06-19
"""
import scanpy as sc
import pandas as pd
import argparse
sc.settings.autoshow = False

import matplotlib
matplotlib.use('agg')

# parse argumetns
parser = argparse.ArgumentParser()
parser.add_argument("--adata_object",
                    default="/Users/crg/Documents/Projects/combat/data/scanpy_pipe_test/pbmc3k.h5ad",
                    help="file name, format: .h5ad")
parser.add_argument("--adata_layer",
                    default="X",
                    help="layer containing data to be plotted, default=X")
parser.add_argument("--umap_coords", default=None, help="txt of umap cords, if not already found in the adata object")
parser.add_argument("--variables_file",
                    default="/Users/crg/Documents/Projects/combat/data/Marker_Lists-myeloid.csv",
                    help="marker_list csv, one column of gene_ids")
parser.add_argument("--figure_suffix", default="custom_markers.png",
                    help="figures filename suffix to be appended to figures/*_")
parser.add_argument("--figure_dir", default="figures/",
                    help="figures path")

args = parser.parse_args()
sc.settings.figdir = args.figure_dir

# script
adata = sc.read_h5ad(args.adata_object)
#adata = sc.read_h5ad("../../data/first_run/anndata-n5000_nneigh15_drharmony_metcosine_neighbors.h5ad")

if args.adata_layer != 'X':
    adata.X = adata.layers[args.adata_layer]

umap = pd.read_csv(args.umap_coords, sep='\t', index_col=0)
#umap = pd.read_csv("/well/combat/shared/citeseq/gex/lowdepth/1.0/baseline_analysis/harmony.seurat.dir/50_1_1_wilcox/umap.dir/umap.txt", sep='\t', header=None, index_col=None)
if "X_umap" in adata.obsm.keys(): 
    print("overwriting the calculated coordinates")
    adata.obsm["X_umap"] = umap.to_numpy()
else:
    adata.obsm["X_umap"] = umap.to_numpy()

fetches = pd.read_csv(args.variables_file, header=None)[0].tolist()
#fetches = pd.read_csv("../../data/Minimal_marker_list_for_full_data.csv", header=None)[0].tolist()

plot_genes = [gg for gg in fetches if gg in adata.var_names]
plot_genes = list(set(plot_genes))
print("plotting")
print(plot_genes)


sc.pl.umap(adata, color=plot_genes, save=args.figure_suffix)
