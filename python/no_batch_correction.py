
import numpy as np
import pandas as pd
import scanpy as sc
import os
import argparse


# parse arguments
parser = argparse.ArgumentParser()

parser.add_argument('--input_anndata',
                    default='adata_scaled.h5ad',
                    help='')
parser.add_argument('--output_csv', default='batch_correction/umap_bc_none.csv',
                    help='')
parser.add_argument('--integration_col', default='batch')
parser.add_argument('--n_pcs', default=40,
                    help="n_comps")
parser.add_argument('--n_neighbors', default=40,
                    help="n_neighbours")

args = parser.parse_args()

print("reading data and starting integration pipeline with script: ")
print(os.path.basename(__file__))

# Scanorama is designed to be used in scRNA-seq pipelines downstream of noise-reduction methods,
# including those for imputation and highly-variable gene filtering.

adata = sc.read(args.input_anndata)

columns = [x.replace (" " ,"") for x in args.integration_col.split(",")]
if len(columns)>1: 
    comb_columns = "|".join(columns)
    adata.obs[comb_columns] = adata.obs[columns].apply(lambda x: '|'.join(x), axis=1)
    columns += [comb_columns]

# write out batch
adata.obs[columns].to_csv(os.path.join(os.path.dirname(args.output_csv), 'batch_mtd.csv'))

# run neighbours and umap without batch correction
sc.pp.neighbors(adata, n_pcs=int(args.n_pcs), n_neighbors=int(args.n_neighbors), use_rep='X_pca')
sc.tl.umap(adata)
print("done umap, saving stuff")
#write out
umap = pd.DataFrame(adata.obsm['X_umap'], adata.obs.index)
umap.to_csv(args.output_csv)

print("done")

#
# #produce dataframes
# PCA_LISI = pd.DataFrame(adata.obsm['X_pca'])
# batch = pd.DataFrame(adata.obs['batch'])
# umap = pd.DataFrame(adata.obsm['X_umap'])
#
# label_colnames = batch.columns
# lisi = hm.compute_lisi(umap, batch, label_colnames) #compute LISI - then the aim will be to plot a density plot to assess distribution of LISI plots so this needs more work
# lisi = pd.DataFrame(lisi,adata.obs.index)
#
# a = sns.distplot(lisi, hist = False, kde = True,
#                  kde_kws = {'linewidth': 3},
#                  label = 'no correction')
# a.set(xlabel="LISI Score", ylabel = "Density")
# a.figure.savefig("output/fig/LISI_nocorrection.png", dpi = 150)
#
# #output LISI
# #lisi.to_csv(args.output_dir +'/'+'lisi_no_batch_correction.csv')
# pd.DataFrame(lisi).to_csv(args.output_dir +'/'+'lisi_no_batch_correction.csv')
#
#
# #produce k-BET
# umap = pd.DataFrame(adata.obsm['X_umap'],adata.obs.index)
# batch = adata.obs['batch']
#
# #umap.to_csv(args.output_dir +'/'+'umap_coordinates_no_batch_correction.csv')
# #batch.to_csv(args.output_dir +'/'+'batch_no_batch_correction.csv')
#
# pd.DataFrame(umap).to_csv(args.output_dir +'/'+'umap_coordinates_no_batch_correction.csv')
# pd.DataFrame(batch).to_csv(args.output_dir +'/'+'batch_forbatchcorrectioncomparison.csv')
#
# ###all of the following needs to be done in one 'cell' - conda R environment will need kBET installed
# #%%R -i umap -i batch -o vals
# #library(kBET)
# #batch.estimate <- kBET(umap, batch, do.pca=FALSE)
# #vals = batch.estimate[['summary']][2:4,2]
#
# #nobatchcorrection_ci = vals #save this for later use - this contains kBET score with 95% CIs essentially for no batch correction - we can plot all together later
#
#
# #adata.write(args.output_dir +'/'+'adatanobatchcorrection.h5ad') - UNDO THIS TO WRITE OUT ANNDATA OBJECT
#
