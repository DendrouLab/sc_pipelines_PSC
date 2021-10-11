import pandas as pd
import scanpy as sc
import argparse
import pickle
import bbknn
import os 
import sys
import logging

## there is a problem with numvba pushong a warning about TBB which is causing the pipeline to crash
# we might need to need to suppress this warning

# parse arguments
parser = argparse.ArgumentParser()

parser.add_argument('--input_anndata',
                    default='adata_scaled.h5ad',
                    help='')
parser.add_argument('--output_csv', default='batch_correction/umap_bc_bbknn.csv',
                    help='')
parser.add_argument('--integration_col', default='batch',
                    help='')
parser.add_argument('--neighbors_within_batch', default='3',
                    help='How many top neighbours to report for each batch; total number of neighbours will be this number times the number of batches.')


args = parser.parse_args()

L = logging.getLogger("bbknn_integration")

L.info("reading data and starting integration pipeline with script: ")
L.info(os.path.basename(__file__))
L.info("Running with options: %s", args)

L.info("reading data and starting integration pipeline with script: ")
L.info(os.path.basename(__file__))

adata = sc.read(args.input_anndata)
nnb = int(args.neighbors_within_batch)
# bbknn can't integrate on 2+ variables, so create a fake column with combined information
columns = [x.replace(" ", "") for x in args.integration_col.split(",")]

if len(columns) > 1:
    L.info("using 2 columns to integrate on more variables")
    # comb_columns = "_".join(columns)
    adata.obs["comb_columns"] = adata.obs[columns].apply(lambda x: '|'.join(x), axis=1)

    # make sure that batch is a categorical
    adata.obs["comb_columns"] = adata.obs["comb_columns"].astype("category")
    nn_test = nnb * len(adata.obs["comb_columns"].unique())
    if nn_test > 40:
        nnb = nnb - 1
        if nnb < 1:
            sys.exit("can't work with 0 neighbors_within_batch")
    # run bbknn
    adata = bbknn.bbknn(adata, batch_key="comb_columns", copy=True, neighbors_within_batch=nnb)
else:
    adata.obs[args.integration_col] = adata.obs[args.integration_col].astype("category")
    nn_test = nnb * len(adata.obs[args.integration_col].unique())
    if nn_test > 40:
        nnb = nnb - 1
        if nnb < 1:
            sys.exit("can't work with 0 neighbors_within_batch") 
    # run bbknn
    adata = bbknn.bbknn(adata, batch_key=args.integration_col, copy=True, neighbors_within_batch=nnb ) #calculates the neighbours

L.info("integration run now calculate umap")
sc.tl.umap(adata)
L.info("done umap, saving stuff")
# write out
umap = pd.DataFrame(adata.obsm['X_umap'], adata.obs.index)
umap.to_csv(args.output_csv)

# save the neighbours to put back into the anndata object
# if this is our fave batch correction
# with open('tmp/bbknn_neighbors.pickle', 'wb') as fp:
#     pickle.dump(adata.uns['neighbors'], fp)

# save full bbknn anndata in tmp, cause need more than just neighbors to work 
outfiletmp = "tmp/bbknn_scaled_adata.h5ad" 

L.info("saving full adata")
adata.write(outfiletmp)

L.debug("done")
