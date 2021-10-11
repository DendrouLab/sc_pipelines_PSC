import pandas as pd
import scanpy as sc
import harmonypy as hm
import logging
import argparse
import h5py
import os

# parse arguments
parser = argparse.ArgumentParser()

parser.add_argument('--input_anndata',
                    default='adata_scaled.h5ad',
                    help='')
parser.add_argument('--output_csv', default='batch_correction/umap_bc_harmony.csv',
                    help='')
parser.add_argument('--integration_col', default='batch',
                    help='')
parser.add_argument('--n_pcs', default=40,
                    help="n_pcs")
parser.add_argument('--n_neighbors', default =40,
                    help="n_neighbours")
parser.add_argument('--sigma_val', default =0.1,
                    help="sigma")
parser.add_argument('--plotconvergence', default =False,
                    help="plot convergence")


args = parser.parse_args()
logging.basicConfig(filename="logs/2_bc_harmony.log",
                            filemode='a',
                            format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s',
                            datefmt='%H:%M:%S',
                            level=logging.INFO)

L = logging.getLogger("Harmony_integration")

L.info("reading data and starting integration pipeline with script: ")
L.info(os.path.basename(__file__))
L.info("Running with options: %s", args)

# this should be an object that contains the full log normalised data (adata_log1p.h5ad)
# prior to hvgs and filtering
adata = sc.read(args.input_anndata)

# Harmony can integrate on 2+ variables,
# but for consistency with other approaches create a fake column with combined information
columns = [x.replace (" " ,"") for x in args.integration_col.split(",")]

if len(columns)>1: 
    L.info("using 2 columns to integrate on more variables")
    #comb_columns = "_".join(columns)
    adata.obs["comb_columns"] = adata.obs[columns].apply(lambda x: '|'.join(x), axis=1)

    # make sure that batch is a categorical
    adata.obs["comb_columns"] = adata.obs["comb_columns"].astype("category")
    # run harmony
    ho = hm.run_harmony(adata.obsm['X_pca'][:,0:int(args.n_pcs)], adata.obs, ["comb_columns"], sigma = float(args.sigma_val) , verbose=True,max_iter_kmeans=30, max_iter_harmony=40, plot_convergence= args.plotconvergence )
else:
    # make sure that batch is a categorical
    adata.obs[args.integration_col] = adata.obs[args.integration_col].astype("category")
    # run harmony
    ho = hm.run_harmony(adata.obsm['X_pca'][:,0:int(args.n_pcs)], adata.obs, [args.integration_col], sigma = float(args.sigma_val) , verbose=True,max_iter_kmeans=30, max_iter_harmony=40, plot_convergence= args.plotconvergence )
L.info("integration run now calculate umap")

adjusted_pcs = pd.DataFrame(ho.Z_corr).T
adata.obsm['X_harmony']=adjusted_pcs.values

sc.pp.neighbors(adata, n_pcs=int(args.n_pcs), n_neighbors=int(args.n_neighbors), use_rep='X_harmony')
sc.tl.umap(adata)
L.info("done umap, saving stuff")
# write out
umap = pd.DataFrame(adata.obsm['X_umap'], adata.obs.index)
umap.to_csv(args.output_csv)

# save the harmony dim reduction in case harmony is our favourite
# consider removing these objects as you need the entire obsm and connectivities anyway - 
# no reason to make a fuss just save the entire adata

# hf = h5py.File('tmp/harmony_obsm.h5.gz', 'w')
# hf.create_dataset('harmony_obsm', data=adata.obsm['X_harmony'], compression="gzip")
# hf.close()

adata.write("tmp/harmony_scaled_adata.h5ad")


L.debug("done")