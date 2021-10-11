import numpy as np
import scanpy as sc
import argparse
import h5py
import pickle
import sys
import os
import pandas as pd
# The purpose of this script is to merge the batch correction of choice into the anndata objects before it goes into
# pipeline_scanpy_clustering.py


# parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('--scaled_anndata',
                    default='adata_scaled.h5ad',
                    help='')
parser.add_argument('--correction_choice', default=None,
                    help='')


args = parser.parse_args()

#
bc = args.correction_choice
# bc = "harmony"

# exit if no batch correction
if bc == "None":
    bc = None

if bc is None:
    sys.exit(0)

# check the batch correction is one of the expected options.
if bc not in ["harmony", "scanorama", "combat", "bbknn"]:
    sys.exit("Batch correction not found, choose from 'combat', 'scanorama', 'bbknn', 'harmony'")

scaled_path = args.scaled_anndata
adata = sc.read(scaled_path)

if bc == "harmony":
    # hf2 = h5py.File(os.path.join(tmp_path, 'harmony_obsm.h5.gz'), 'r')
    # harmony_df = hf2.get('harmony_obsm')
    # adata.obsm['X_harmony'] = np.array(harmony_df)
    # hf2.close()
    adata = sc.read("tmp/harmony_scaled_adata.h5ad")  
    adata.write(scaled_path)
elif bc == "scanorama":
    # hf2 = h5py.File(os.path.join(tmp_path, 'scanorama_obsm.h5.gz'), 'r')
    # sc_df = hf2.get('scanorama_obsm')
    # print(np.array(sc_df).shape)
    # adata.obsm['X_scanorama'] = np.array(sc_df)
    # hf2.close()
    adata = sc.read("tmp/scanorama_scaled_adata.h5ad")
    adata.write(scaled_path)
elif bc == "bbknn":
    # bbknn_obj = pickle.load(open(os.path.join(tmp_path, "bbknn_neighbors.pickle"), "rb"))
    # adata.uns['neighbors'] = bbknn_obj
    # re-read anndata from tmp
    adata = sc.read("tmp/bbknn_scaled_adata.h5ad") 
    adata.write(scaled_path)
# elif bc == "combat":
#     adata = sc.read(log1p_path)
#     hf2 = h5py.File(os.path.join(tmp_path, 'combat_X.h5.gz'), 'r')
#     combat_X = hf2.get("X")
#     print(np.array(combat_X).shape)
#     adata.X = np.array(combat_X)
#     hf2.close()
#     # note we need to run normalise_log_hvg_regress_scale.py again, it's a pain.
#     # calculate hvgs
#     if args.flavor == "seurat_v3":
#         if args.n_top_genes is not None:
#             sc.pp.highly_variable_genes(adata, flavor="seurat_v3",
#                                         n_top_genes=args.n_top_genes)
#         else:
#             sc.pp.highly_variable_genes(adata, flavor="seurat_v3")
#     else:
#         sc.pp.highly_variable_genes(adata, flavor=args.flavor,
#                                     min_mean=args.min_mean, max_mean=args.max_mean,
#                                     min_disp=args.min_disp)
#     # exclude hvgs if there is an exclude file
#     if args.exclude_file is not None:
#         excl = pd.read_csv(args.exclude_file, delimiter='\t')
#         excl.status.value_counts()
#
#         print("number of hvgs prior to filtering")
#         print(adata.var.highly_variable.sum())
#
#         adata.var['hvg_exclude'] = [True if gg in excl.gene_id.tolist() else False for gg in adata.var.gene_ids]
#         adata.var.loc[adata.var.hvg_exclude, 'highly_variable']= False
#
#         print("number of hvgs after filtering")
#         print(adata.var.highly_variable.sum())
#     # save object (has to occur before scaling)
#     adata.write(log1p_path)
#     adata.raw = adata # this means that
#     # filter by hvgs
#     adata = adata[:, adata.var.highly_variable]
#     # regress out
#     if args.regress_out is not None:
#         regress_opts = args.regress_out.split(",")
#         sc.pp.regress_out(adata, regress_opts)
#     if args.scale_max_value is not None:
#         sc.pp.scale(adata)
#     else:
#         sc.pp.scale(adata, max_value=args.scale_max_value)
#     # run pca
#     sc.tl.pca(adata, n_comps=250, svd_solver='arpack', random_state=0) #given args above this should work
#     # extract pca coordinates for plotting (in R??)
#     pca_coords = pd.DataFrame(adata.obsm['X_pca'])
#     # add in the rownames
#     pca_coords.index = adata.obs_names
#     # save coordinates to file
#     # (note this saves values values up to 6 significant figures, because why save 20 for a plot
#     pca_coords.to_csv(args.output_prefix + "combat_pca.txt.gz", sep='\t')
#     # save the scaled adata object
#     adata.write(scaled_path)
