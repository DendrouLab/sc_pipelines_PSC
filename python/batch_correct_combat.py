import pandas as pd
import scanpy as sc
import argparse
import h5py
import os
import logging


# parse arguments
parser = argparse.ArgumentParser()

parser.add_argument('--input_anndata',
                    default='adata_log1p.h5ad',
                    help='')
parser.add_argument('--output_csv', default='batch_correction/umap_bc_combat.csv',
                    help='')
parser.add_argument('--integration_col', default='batch', help='')
parser.add_argument('--n_pcs', default=40,
                    help="n_pcs")
parser.add_argument('--n_neighbors', default =40,
                    help="n_comps")

args = parser.parse_args()

logging.basicConfig(filename="logs/2_bc_harmony.log",
                            filemode='a',
                            format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s',
                            datefmt='%H:%M:%S',
                            level=logging.INFO)

L = logging.getLogger("combat_integration")

L.info("reading data and starting integration pipeline with script: ")
L.info(os.path.basename(__file__))
L.info("Running with options: %s", args)

# this should be an object that contains the full log normalised data (adata_log1p.h5ad)
# prior to hvgs and filtering
adata = sc.read(args.input_anndata)

# combat can't integrate on 2+ variables, so create a fake column with combined information
columns = [x.replace(" ", "") for x in args.integration_col.split(",")]
if len(columns) > 1: 
    L.info("using 2 columns to integrate on more variables")
    #comb_columns = "_".join(columns)
    adata.obs["comb_columns"] = adata.obs[columns].apply(lambda x: '|'.join(x), axis=1)

    # make sure that batch is a categorical
    adata.obs["comb_columns"] = adata.obs["comb_columns"].astype("category")
    # run combat
    sc.pp.combat(adata, key="comb_columns")
    
else:

    # run combat
    sc.pp.combat(adata, key=args.integration_col)

L.info("integration run now calculate umap")

sc.pp.highly_variable_genes(adata)
sc.pp.pca(adata, n_comps=int(args.n_pcs), use_highly_variable=True, svd_solver='arpack') #use same args as plot_pca.py
sc.pp.neighbors(adata, n_pcs=int(args.n_pcs), n_neighbors=int(args.n_neighbors))
sc.tl.umap(adata)
L.info("done umap, saving stuff")
#write out
umap = pd.DataFrame(adata.obsm['X_umap'], adata.obs.index)
umap.to_csv(args.output_csv)

adata.write("tmp/combat_scaled_adata.h5ad")

L.debug("done")