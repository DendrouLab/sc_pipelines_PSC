"""
# reads in an anndata object with neighbours already computed, runs scanpy umap
# save umap coordinates (for plotting elsewhere
# CRG 2020-06-12
"""
import scanpy as sc
import numpy as np
import pandas as pd
import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument("--infile", default="/Users/crg/Documents/Projects/combat/data/scanpy_pipe_test/pbmc3k_nn30_drharmony_metcosine_neighbors.h5ad", help="file name, format: .h5ad")
parser.add_argument("--outfile", default="/Users/crg/Documents/Projects/combat/data/scanpy_pipe_test/clusters.h5ad", help="file name, format: .h5ad")
parser.add_argument("--min_dist", default=0.1, help="no. neighbours parameters for sc.pp.neighbors()")

args = parser.parse_args()

# set seed
# seed = int(200612)

# read data
adata = sc.read_h5ad(args.infile)

# check sc.pp.neihgbours has been run
if 'neighbors' not in adata.uns.keys():
    sys.exit("Error: sc.pp.neighbours has not been run on this object")

# what parameters?
sc.tl.umap(adata, min_dist=float(args.min_dist))

# extract umap coordinates for plotting (in R??)
umap_coords = pd.DataFrame(adata.obsm['X_umap'])

# add in the rownames 
umap_coords.index = adata.obs_names

# save coordinates to file
# (note this saves values values up to 6 significant figures, because why save 20 for a plot
umap_coords.to_csv(args.outfile, sep = '\t')
