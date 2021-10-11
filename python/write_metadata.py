import argparse
import numpy as np
import pandas as pd
import scanpy as sc

parser = argparse.ArgumentParser()
parser.add_argument('--infile',
                    default="../../../data/pipe_test/sc_full.h5ad",
                    help="anndata object that has benn QCd, normalised and logged")
parser.add_argument('--outfile', default='../../../data/pipe_test/sc_preprocess.h5ad',
                    help="anndta object that has been filtered by hvgs, \
                    regressed, scaled and batch corrected with harmony, \
                    file name, format: .h5ad")
args = parser.parse_args()

# read file
adata = sc.read_h5ad(args.infile)

adata.obs.to_csv(args.outfile, sep='\t')