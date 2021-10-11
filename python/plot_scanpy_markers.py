"""
Fetching data for plotting
CRG 2020-06-19
"""
import scanpy as sc
import pandas as pd
import argparse
import os
import re
sc.settings.autoshow = False

import matplotlib
matplotlib.use('agg')

# parse argumetns
parser = argparse.ArgumentParser()
parser.add_argument("--adata_object",
                    default="ann-data-log1p.h5ad, full object, not filtered by HVG",
                    help="file name, format: .h5ad")
parser.add_argument("--adata_layer",
                    default="X",
                    help="layer containing data to be plotted, default=X")
parser.add_argument("--marker_file",
                    default="",
                    help="markers from sc.rank_genes_group() as txt(.gz)")
parser.add_argument("--cluster_file",
                    default="",
                    help="cluster info as txt (for the same resolution as marker_files, one column")
parser.add_argument("--figure_prefix", default="figures/markers_",
                    help="figures path")
parser.add_argument("--n", default=5,
                    help="number of genes per cluster")


args = parser.parse_args()
sc.settings.figdir = args.figure_prefix

# script
adata = sc.read_h5ad(args.adata_object)
if args.adata_layer != 'X':
    adata.X = adata.layers[args.adata_layer]

# run pca if not done already.( for dendorgram calc)
if("X_pca" in adata.obsm.keys()):
    sc.pp.pca(adata)


clusters = pd.read_csv(args.cluster_file, sep='\t', index_col=0)
# convert clusters to type category
clusters = clusters.apply(lambda x: x.astype('category'))
# merge clusters with adata
adata.obs = adata.obs.merge(clusters, left_index=True, right_index=True)

# read file
mf = pd.read_csv(args.marker_file, sep='\t')

df = mf.groupby('cluster').apply(lambda x: x.nlargest(args.n, ['scores'])).reset_index(drop=True)
marker_list={str(k): list(v) for k,v in df.groupby("cluster")["gene"]}


sc.pl.stacked_violin(adata,
              marker_list,
              groupby='clusters',
              save='top_markers.png',
              figsize=(24, 5))

sc.pl.matrixplot(adata,
                 marker_list,
                 groupby='clusters',
                 save='top_markers.png',
                 figsize=(24, 5))

sc.pl.dotplot(adata,
              marker_list,
              groupby='clusters',
              save='top_markers.png',
              figsize=(24, 5))
