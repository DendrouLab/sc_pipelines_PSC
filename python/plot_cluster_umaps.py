"""
Scanpy inbuilt plots for custom markers
CRG 2020-06-19
"""

import scanpy as sc
sc.settings.autoshow = False
import pandas as pd
import argparse

import matplotlib
matplotlib.use('agg')

# parse argumetns
parser = argparse.ArgumentParser()
parser.add_argument("--adata_object",
                    default="/Users/crg/Documents/Projects/combat/data/scanpy_pipe_test/pbmc3k.h5ad",
                    help="file name, format: .h5ad")
parser.add_argument("--clusters_file",
                    default="/Users/crg/Documents/Projects/combat/data/scanpy_pipe_test/pbmc3k_nn30_algleiden_clusters_list.txt.gz",
                    help="compiled clusters dataframe")
parser.add_argument("--umap_file",
                    default="/Users/crg/Documents/Projects/combat/data/scanpy_pipe_test/pbmc3k_nn30_drharmony_metcosine_md0.1_umap.txt",
                    help = "txt file containing umap coordinates")
parser.add_argument("--figure_dir",
                    default="figures",
                    help="figures directory")
parser.add_argument("--figure_suffix",
                    default="_pbmc3k_nn30_algleiden.pdf",
                    help="figures filename suffix to be appended to figures/umap")
args = parser.parse_args()

sc.settings.figdir = args.figure_dir

adata = sc.read_h5ad(args.adata_object)
clusters = pd.read_csv(args.clusters_file, sep='\t', index_col=0)
umap = pd.read_csv(args.umap_file, sep='\t', index_col=0)


# convert clusters to type category
clusters = clusters.apply(lambda x: x.astype('category'))
# merge clusters with adata
adata.obs = adata.obs.merge(clusters, left_index=True, right_index=True)
res_cols = clusters.columns

# add umap coords
adata.obsm['X_umap'] = umap.to_numpy()
# plot umap

# auto saving which is annoying
sc.pl.umap(adata, color=res_cols, show=False, legend_loc='on data', save=args.figure_suffix)

