## Filtering anndata based on n_genes, percent mito and min cells/
## originally written by Tom Thomas (https://github.com/tomthomas3000/TAURUS)
## adapted and augmented for this pipeline by Charlotte Rich-Griffin 2020-09-30

import scanpy as sc
import argparse
import logging
import sys
import pandas as pd
logging.basicConfig(format="[ %(asctime)s: %(levelname)s: %(message)s ]", level=logging.INFO, stream=sys.stdout)
L = logging.getLogger(__name__)


sc.settings.verbosity = 2  # verbosity: errors (0), warnings (1), info (2), hints (3)
# comment this because of numba issues
# sc.logging.print_versions()

# parse arguments
parser = argparse.ArgumentParser()

parser.add_argument('--input_anndata',
                    default='data/anndata-n5000-unfilt.h5ad',
                    help='')
parser.add_argument('--output_anndata',
                    default='data/anndata-n5000-filt.h5ad',
                    help='')
parser.add_argument('--min_genes', default=0,
                    help='')
parser.add_argument('--max_genes', default="inf",
                    help='')

parser.add_argument('--min_cells', default=0,
                    help='')

parser.add_argument('--max_counts', default="inf",
                    help='')

parser.add_argument('--percent_mito', default=100,
                    help='exclude any cells with >n% mitochondrial content (pct_counts_mt column) default=100 (no filtering)')
parser.add_argument('--percent_ribo', default=100,
                    help='exclude any cells with >n% ribosomal content (percent.ribo column), default=100 (no filtering)')
parser.add_argument('--percent_hb', default=100,
                    help='exclude any cells with >n% haemoglobin content (percent.ribo column), default=100 (no filtering)')

parser.add_argument('--drop_nas_col', default=None,
                    help='')
L.info("Running filter_anndata")

parser.set_defaults(verbose=True)
args = parser.parse_args()

adata = sc.read(args.input_anndata)

###filter cells, genes, and mitochondria - NOTE: might merit further screening if need to finetune
n_obs = adata.n_obs
L.info("Pre-filter number of cells %d" % adata.n_obs)

sc.pp.filter_cells(adata, min_genes=int(args.min_genes))
sc.pp.filter_genes(adata, min_cells=int(args.min_cells))

# exclude cells with higher total counts or genes than specified (default is inf)

adata = adata[adata.obs['total_counts'] < float(args.max_counts), :]
adata = adata[adata.obs['n_genes'] < float(args.max_genes), :]

if 'pct_counts_mt' in adata.obs.columns:
    if float(args.percent_mito) < 1:
        msg = """percent mito argument is below 1%, 
        this is inadvisable as you will lose the majority of your data!, 
        suggested values are in the range of 5-50%"""
        L.error(msg)
        raise ValueError(msg)
    adata = adata[adata.obs['pct_counts_mt'] < float(args.percent_mito), :]

if 'pct_counts_ribo' in adata.obs.columns:
    if float(args.percent_ribo) < 1:
        msg="""percent ribo argument is below 1%, 
            this is inadvisable as you will lose the majority of your data!, 
            suggested values are in the range of 5-50%"""
        L.error(msg)
        raise ValueError(msg)
    adata = adata[adata.obs['pct_counts_ribo'] < float(args.percent_ribo), :]

if 'pct_counts_hb' in adata.obs.columns:
    if float(args.percent_hb) < 1:
        msg="""percent hb argument is below 1%, 
            this is inadvisable as you will lose the majority of your data!, 
            suggested values are in the range of 5-50%"""
        L.error(msg)
        raise ValueError(msg)
    adata = adata[adata.obs['pct_counts_hb'] < float(args.percent_hb), :]

# drop demultiplexing data without annotation
col_arg=args.drop_nas_col
# print(adata)
if col_arg is not None:
    # check if it refers to multiple cols
    col_choices=col_arg.split(',') # returns list
    # sequentially remove NAs, (dropna not an option because slicing the whole anndata)
    n_obs = adata.n_obs
    for col in col_choices:
        L.debug("Filtering nans from %s column" % col)
        if col in adata.obs.columns:
            L.debug(col)
            adata = adata[adata.obs[col].notna(),:]
        else:
            msg="""demultiplexing column not found in data, check inputs"""
            L.error(msg)
            raise ValueError(msg)
else:
    pass

n_obs = n_obs - adata.n_obs
L.info("No. cells removed for being NA %d" % n_obs)
# print(adata)

L.info("Remaining cells %d" % adata.n_obs)

adata.write(args.output_anndata)

#This stage is the point (i.e. pre-normalisation) where the adata file can be outputted so that we can extract raw matrix for the cellphonedb.
L.info("Completed")
