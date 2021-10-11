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
                    default="/Users/crg/Documents/Projects/combat/data/scanpy_pipe_test/pbmc3k.h5ad",
                    help="file name, format: .h5ad")
parser.add_argument("--adata_layer",
                    default="X",
                    help="layer containing data to be plotted, default=X")
parser.add_argument("--marker_files",
                    default="/Users/crg/Documents/Projects/combat/data/Marker_Lists-myeloid.csv",
                    help="comma separateed list of marker_list csvs, one column of gene_ids")
parser.add_argument("--clusters_file",
                    default="",
                    help="cluster info, one resolution per column")
parser.add_argument("--n_genes",
                    default=10,
                    help="how many top genes to plot")
parser.add_argument("--figure_dir", default="figures/",
                    help="figures path")

args = parser.parse_args()
sc.settings.figdir = args.figure_dir

# script
adata = sc.read_h5ad(args.adata_object)
if args.adata_layer != 'X':
    adata.X = adata.layers[args.adata_layer]

# run pca if not done already.( for dendorgram calc)
if("X_pca" in adata.obsm.keys()):
    sc.pp.pca(adata)

clusters = pd.read_csv(args.clusters_file, sep='\t', index_col=0)

# convert clusters to type category
clusters = clusters.apply(lambda x: x.astype('category'))
# merge clusters with adata
adata.obs = adata.obs.merge(clusters, left_index=True, right_index=True)
res_cols = clusters.columns

marker_files = args.marker_files.split(',')
for mf in marker_files:
    uniq_id = os.path.basename(mf)
    uniq_id = os.path.splitext(uniq_id)[0]
    fetches = pd.read_csv(mf, header=None)[0].tolist()

    plot_genes = [gg for gg in fetches if gg in adata.var_names]
    plot_genes = list(set(plot_genes))

    not_found = [gg for gg in fetches if gg not in adata.var_names]
    not_found_file = uniq_id + "not_found.txt"
    with open(not_found_file, 'w') as f:
        for item in not_found:
            f.write("%s\n" % item)

    sc.tl.dendrogram(adata, 'clusters')
    sc.pl.dotplot(adata,
                  var_names=plot_genes,
                  groupby='clusters',
                  dendrogram=True,
                  save=uniq_id + '.png',
                  figsize=(24, 5))

    sc.pl.matrixplot(adata,
                     var_names=plot_genes,
                     groupby='clusters',
                     dendrogram=True,
                     save=uniq_id + '.png',
                     figsize=(24, 5))
