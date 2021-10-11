
import scanpy as sc
import argparse
import logging
import sys
import pandas as pd
import os


def test_file_or_value(test):
    try:
        float(test)
        out = "value"
    except:
        if os.path.isfile(test):
            out = "file"
        else:
            print("no scrublet thresholds set")
            out = None
            pass
    return out


def merge_with_adata_obs(obs, df, on_col="sample_id"):
    obs['barcode_id'] = obs.index
    new_df = adata.obs.merge(scr_df, how="left", on=on_col)
    new_df.index = new_df['barcode_id']
    new_df.index.name = None
    return new_df


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
parser.add_argument('--output_prefix',
                    default='',
                    help='prefix to prepend to saved files. here just used for metadata file')
parser.add_argument('--scrublet_cutoff', default=None,
                    help='this is either a value or a file')
parser.add_argument('--batch_col', default="sample_id",
                    help='')

L.warning("Running subtract scrublet")

parser.set_defaults(verbose=True)
args = parser.parse_args()
# load data
adata = sc.read(args.input_anndata)

# test what kind of input we are dealing with, value or file?
run = test_file_or_value(args.scrublet_cutoff)


# adata=sc.read_h5ad("./data/run_integration/taurus_test_filt.h5ad")

# run = test_file_or_value("./data/run_integration/scrublet_cutoffs.csv")
# subtract scrublet
if run == 'value':
    # if args.scrublet_cutoff is not None:
    # if it is a value
    adata = adata[adata.obs['doublet_scores'] < float(args.scrublet_cutoff), :]
elif run == 'file':
    # if it is a file
    scr_df = pd.read_csv(args.scrublet_cutoff, names=["sample_id", "scrublet_cutoff"])
    adata.obs = merge_with_adata_obs(adata.obs, scr_df, on_col="sample_id")
    adata = adata[adata.obs['doublet_scores'] < adata.obs['scrublet_cutoff'],]
else:
    pass

L.info("saving anndata and obs in a metadata tsv file")
metafile = adata.obs
metafile["cellbarcode"] = adata.obs.index

savename= args.output_prefix + "_filtered_cell_metadata.tsv"
metafile.to_csv(savename, sep='\t', index=True)  

adata.write(args.output_anndata)

L.info("done")