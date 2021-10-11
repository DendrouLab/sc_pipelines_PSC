##
# Scanorama is designed to be used in scRNA-seq pipelines downstream of noise-reduction methods,
# including those for imputation and highly-variable gene filtering.
##

import pandas as pd
import scanpy as sc
# import scanorama
import scanpy.external as sce
import argparse
import sys
import os
import logging

# parse arguments
parser = argparse.ArgumentParser()

parser.add_argument('--input_anndata',
                    default='adata_scaled.h5ad',
                    help='')
parser.add_argument('--output_csv', default='batch_correction/umap_bc_scanorama.csv',
                    help='')
parser.add_argument('--integration_col', default='batch',
                    help='')
parser.add_argument('--n_pcs', default=40,
                    help="n_comps")
parser.add_argument('--n_neighbors', default=40,
                    help="n_neighbours")
args = parser.parse_args()

logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)
L = logging.getLogger("scanorama_integration")

L.info("reading data and starting integration pipeline with script: ")
L.info(os.path.basename(__file__))
L.info("Running with options: %s", args)

# Scanorama is designed to be used in scRNA-seq pipelines downstream of noise-reduction methods,
# including those for imputation and highly-variable gene filtering.
L.info("reading data and starting integration pipeline with script: ")
L.info(os.path.basename(__file__))

adata = sc.read(args.input_anndata)
bcs = adata.obs_names.tolist()

# scanorama can't integrate on 2+ variables, so create a fake column with combined information
columns = [x.replace(" " , "") for x in args.integration_col.split(",")]
if len(columns) > 1:
    L.info("using 2 columns to integrate on more variables")
    # comb_columns = "_".join(columns)
    adata.obs["comb_columns"] = adata.obs[columns].apply(lambda x: '|'.join(x), axis=1)

    # make sure that batch is a categorical
    adata.obs["comb_columns"] = adata.obs["comb_columns"].astype("category")
    adata.obs['batch'] = adata.obs["comb_columns"]
else:
    if args.integration_col != "batch":
        adata.obs['batch'] = adata.obs[args.integration_col]

batches = adata.obs['batch'].unique()

# need contiguous batches
scanorama_input = adata[adata.obs.sort_values(by="batch").index.tolist(), :]

# filter by HVGs to make it equivalent to the old scripts,
# which inputted the scaled object after filtering by hvgs.
scanorama_input = scanorama_input[:, scanorama_input.var.highly_variable]
# run scanoramam using the scanpy integrated approach
sce.pp.scanorama_integrate(scanorama_input, key='batch')
# not integrated

# old method (simplified)
# alldata = {}
# for bb in batches:
#     alldata[bb] = scanorama_input[scanorama_input.obs['batch'] == bb]
#
# adatas = list(alldata.values())
# X_scanorama = scanorama.integrate_scanpy(adatas)
# scanorama_input.obsm['X_scanorama'] = np.concatenate(X_scanorama)

# put into the original order
scanorama_input = scanorama_input[bcs, :]
# check it worked
if all(scanorama_input.obs_names == bcs):
    L.info("barcode order is correct")
else:
    L.debug("barcode order in  batch corrected object is incorrect")
    sys.exit("barcode order in  batch corrected object is incorrect")


L.info(adata)
# add to the AnnData object
adata.obsm["X_scanorama"] = scanorama_input.obsm["X_scanorama"]

L.info(adata)
# add to the AnnData object
# adat.obsm["X_scanorama"] = all_anndata
L.info("integration run now calculate umap")
# umap
sc.pp.neighbors(adata, n_pcs=int(args.n_pcs), n_neighbors=int(args.n_neighbors), use_rep="X_scanorama")
sc.tl.umap(adata)
L.info("done umap, saving stuff")
# write out
umap = pd.DataFrame(adata.obsm['X_umap'], adata.obs.index)
umap.to_csv(args.output_csv)

# save the scanorama dim reduction in case scanorama is our favourite
adata.write("tmp/scanorama_scaled_adata.h5ad")
L.debug("done")
