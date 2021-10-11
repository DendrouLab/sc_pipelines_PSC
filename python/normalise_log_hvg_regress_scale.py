## Processing anndata based on n_genes, percent mito and min cells/
## originally written by Tom Thomas (https://github.com/tomthomas3000/TAURUS)
## adapted for this pipeline by Charlotte Rich-Griffin 2020-09-30


import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io
import matplotlib.pyplot as plt
import os
import anndata
import argparse

# parse arguments
parser = argparse.ArgumentParser()


parser.add_argument('--input_anndata',
                    default='data/anndata-n5000_filt.h5ad',
                    help='')
parser.add_argument('--output_prefix',
                    default='data/anndata-n5000',
                    help='')
parser.add_argument('--fig_dir', 
                    default="./figures",
                    help="save plots here")
parser.add_argument('--exclude_file', default=None,
                    help='')
# highly variable genes options
parser.add_argument('--flavor', default='seurat')
parser.add_argument('--n_top_genes', default=None)
parser.add_argument('--min_mean', default=0.0125)
parser.add_argument('--max_mean', default=3)
parser.add_argument('--min_disp', default=0.5)
parser.add_argument("--filter_by_hvg", default=False)
# regress out options
parser.add_argument('--regress_out', default=None)
# scale options
parser.add_argument('--scale_max_value', default=None)
args = parser.parse_args()

sc.settings.verbosity = 3
sc.logging.print_header()

figdir = args.fig_dir

if not os.path.exists(figdir):
    os.mkdir(figdir)

sc.settings.figdir = figdir
sc.set_figure_params(scanpy=True, fontsize=14, dpi=300, facecolor='white', figsize=(5,5))

adata = sc.read(args.input_anndata)
# Normalise to depth 10k, store raw data, assess and drop highly variable genes, regress mitochondria and count

# sc.pp.highly variabel genes Expects logarithmized data, except when flavor='seurat_v3' in which count data is expected.
# change the order accordingly
print("run hvgs")
if args.flavor == "seurat_v3":
    if args.n_top_genes is None:
        print("if seurat_v3 is used you must give a n_top_genes value")
        # sc.pp.highly_variable_genes(adata, flavor="seurat_v3",)
    else:
        sc.pp.highly_variable_genes(adata, flavor="seurat_v3",
                                    n_top_genes=int(args.n_top_genes))
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
else:
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, flavor=args.flavor,
                                min_mean=float(args.min_mean), max_mean=float(args.max_mean),
                                min_disp=float(args.min_disp))

sc.pl.highly_variable_genes(adata,show=False, save ="_genes_highlyvar.png")

print("exlcude genes from hvg")
# exclude hvgs if there is an exclude file
if args.exclude_file is not None:
    excl = pd.read_csv(args.exclude_file, delimiter='\t')
    print(excl.status.value_counts())

    print("number of hvgs prior to filtering")
    print(adata.var.highly_variable.sum())

    #adata.var['hvg_exclude'] = [True if gg in excl.gene_id.tolist() else False for gg in adata.var.index]
    adata.var['hvg_exclude'] = [True if gg in excl.gene_name.tolist() else False for gg in adata.var.index]
    adata.var.loc[adata.var.hvg_exclude, 'highly_variable']= False

    print("number of hvgs after filtering")
    print(adata.var.highly_variable.sum())
    sc.pl.highly_variable_genes(adata,show=False, save ="_exclude_genes_highlyvar.png")

genes = adata.var
genes['gene_name'] = adata.var.index
genes.to_csv("filtered_genes.tsv", sep="\t", index=True)
# save object (has to occur before scaling) 
adata.write(args.output_prefix + "_log1p.h5ad")


adata.raw = adata # this means that log normalised counts are saved in raw
# filter by hvgs
if isinstance(args.filter_by_hvg, str):
    # convert string to bool
    filter_by_hvgs = eval(args.filter_by_hvg)
else:
    filter_by_hvgs = args.filter_by_hvg

if filter_by_hvgs is True:
    adata = adata[:, adata.var.highly_variable]

# regress out
if args.regress_out is not None:
    regress_opts = args.regress_out.split(",")
    sc.pp.regress_out(adata, regress_opts)

if args.scale_max_value is not None:
    sc.pp.scale(adata)
else:
    sc.pp.scale(adata, max_value=args.scale_max_value)

# run pca
sc.tl.pca(adata, n_comps=250, svd_solver='arpack', random_state=0) #given args above this should work
# extract pca coordinates for plotting (in R??)
pca_coords = pd.DataFrame(adata.obsm['X_pca'])
# add in the rownames 
pca_coords.index = adata.obs_names
# save coordinates to file
# (note this saves values values up to 6 significant figures, because why save 20 for a plot
pca_coords.to_csv(args.output_prefix + "_pca.txt.gz", sep='\t')

# save the raw counts

# save the scaled adata object
adata.write(args.output_prefix + "_scaled.h5ad")
