## Filtering anndata based on n_genes, percent mito and min cells/
## originally written by Tom Thomas (https://github.com/tomthomas3000/TAURUS)
## adapted for this pipeline by Charlotte Rich-Griffin 2020-09-30

import scanpy as sc
import argparse
import logging
import sys
from random import sample


def downsample_adata(adata, nn, cat_col=None):
    # check that the downsample is possible
    if cat_col is not None:
        bcs = []
        for ss in adata.obs[cat_col].unique():
            # only sample if you have more cells than nn in a category
            n_cells=len(list(adata.obs.index[adata.obs[cat_col] == ss]))
            if n_cells > nn:
                bcs = bcs + sample(list(adata.obs.index[adata.obs[cat_col] == ss]), nn)
            else:
                bcs = bcs + list(adata.obs.index[adata.obs[cat_col] == ss])
    else:
        #only sample if you havemore cells htan nn
        n_cells = len(list(adata.obs.index))
        if n_cells > nn:
            bcs = sample(list(adata.obs.index), nn)
        else:
            print("cannot downsample as not enough cells")
            bcs = list(adata.obs.index)

    keep = [True if aa in bcs else False for aa in list(adata.obs.index)]
    adata_sub = adata[keep, ]
    return adata_sub


L = logging.getLogger(__name__)
log_handler = logging.StreamHandler(sys.stdout)
log_handler.setFormatter(logging.Formatter('%(asctime)s %(message)s'))
log_handler.setLevel(logging.INFO)
L.addHandler(log_handler)

sc.settings.verbosity = 0  # verbosity: errors (0), warnings (1), info (2), hints (3)
# comment this because of numba issues
# sc.logging.print_versions()

# parse arguments
parser = argparse.ArgumentParser()

parser.add_argument('--input_anndata',
                    default='data/anndata-n5000-filt.h5ad',
                    help='')
parser.add_argument('--output_anndata',
                    default='data/anndata-n5000-filt.h5ad',
                    help='')
parser.add_argument('--downsample_value', default=None,
                    help='')
parser.add_argument('--batch_col', default=None,
                    help='')

L.warning("Running filter_anndata")

parser.set_defaults(verbose=True)
args = parser.parse_args()

adata = sc.read(args.input_anndata)

L.info("Number of cells per sample")
L.info(print(adata.obs.sample_id.value_counts()))

if args.batch_col == "None":
    args.batch_col = None

adata = downsample_adata(adata, nn=int(args.downsample_value), cat_col=args.batch_col)

L.info("Number of cells per sample")
L.info(print(adata.obs.sample_id.value_counts()))

adata.write(args.output_anndata)

#This stage is the point (i.e. pre-normalisation) where the adata file can be outputted so that we can extract raw matrix for the cellphonedb.
L.warning("Completed")