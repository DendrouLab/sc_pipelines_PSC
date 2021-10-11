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

parser.add_argument("--input_anndata",
                    default="adata_scaled.h5ad",
                    help="")
parser.add_argument("--fig_dir", default="figures/")
# pca options
parser.add_argument("--n_pcs", default=50)
parser.add_argument("--color_by", default="batch")

args = parser.parse_args()

sc.settings.autoshow=False
sc.settings.figdir=args.fig_dir # change to whatever the output directory should be - provided that "figures" is a folder within the output directory - this should deposit figures there
sc.settings.set_figure_params(dpi_save=300)

adata = sc.read(args.input_anndata)

# do some plots!
sc.pl.pca_variance_ratio(adata, log=True, n_pcs=int(args.n_pcs), save=".png")

col_variables = args.color_by.split(",")
# for cv in col_variables:
#     sc.pl.pca(adata, color=cv, save="_" + cv + ".png")

sc.pl.pca(adata, color=col_variables, save = "_vars.png")

sc.pl.pca_loadings(adata, components="1,2,3,4,5,6", save = ".png")
sc.pl.pca_overview(adata, save = ".png")