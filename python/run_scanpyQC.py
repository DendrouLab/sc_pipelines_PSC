'''
scanpy QC script
'''
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import argparse
import sys
import scanpy as sc


parser = argparse.ArgumentParser()
# required option
parser.add_argument("--submissionfile",
                    default="caf_example",
                    help="this file has all the samples of the experiment. it has at least 3 columns: \
                    sample_id, raw_10xpath and filtered_10x path. \
                    it can contain other anonymous columns which will be renamed using the metadatacols information")
parser.add_argument("--sampleprefix",
                    default="",
                    help="prefix to prepend when saving the metadata file")
parser.add_argument("--metadatacols",
                    default="patient,tissue,timepoint",
                    help="names of the metadata to rename the columns of the sample submission file ")
parser.add_argument("--outfile",
                    default="adata_unfilt.h5ad",
                    help="")
parser.add_argument("--figdir",
                    default="./figures/",
                    help="path to save the figures to")
parser.add_argument("--figure_suffix",
                    default="_qc-plot.png",
                    help="figures filename suffix to be appended to figures/umap")
parser.add_argument("--scrubletdir",
                    default="./scrublet/",
                    help="path to save the figures to")
parser.add_argument("--iggenesfile",
                    default=None,
                    help="path to file containing IG genes")
parser.add_argument("--ccgenes",
                    default=None,
                    help="path to file containing cell cycle genes")
parser.add_argument("--demultiplexing",
                    default="False",
                    help="is there demultiplexing data to be incorporated")
parser.add_argument("--demultiplexing_metadatacols",
                    default="antibody",
                    help="column names of the demuliplexing file to incoprate")
args = parser.parse_args()

print("Scanpy qc pipeline")

sc.settings.verbosity = 3
sc.logging.print_header()

print("running with args:")
print(args)
figdir = args.figdir

if not os.path.exists(figdir):
    os.mkdir(figdir)

sc.settings.figdir = figdir
sc.set_figure_params(scanpy=True, fontsize=14, dpi=300, facecolor='white', figsize=(5,5))


sfile = args.submissionfile
if not os.path.exists(sfile):
    sys.exit("Did you forget to create a submission file?")

caf = pd.read_csv(sfile, sep="\t")


metadatacols = args.metadatacols.split(',')
print("reading in all data one by one ...")
adatas = [sc.read_10x_mtx(f, var_names='gene_symbols', cache=True) for f in caf['filtered_path']]
# add metadata to each object
for i in range(len(adatas)):
    print(i)
    # add in the extra meatadata if it exists.
    for col in metadatacols:
        adatas[i].obs[col] = caf[col][i]



# add in the demultiplexing data if applicable
if args.demultiplexing == "True":
    demult_metadatacols = args.demultiplexing_metadatacols.split(',')
    for i in range(len(adatas)):
        # check to see if this sample has demultiplexing data
            # this is false id it is a nan
        if caf['demultiplex_map_file'][i] == caf['demultiplex_map_file'][i]:
            print(caf['demultiplex_map_file'][i])
            demultiplex_df = pd.read_csv(caf['demultiplex_map_file'][i], index_col=None, header=0)
            # for speed subset this df by the barcodes in the adata object
            bcs = list(adatas[i].obs_names)
            demultiplex_df = demultiplex_df[demultiplex_df['barcode_id'].isin(bcs)]
            demultiplex_mtd = pd.read_csv(caf['demultiplex_mtd_file'][i], index_col=None, header=0)
            # turn all the columns into categorys, they mostl likely are and it'll ease plotting
            demultiplex_mtd = demultiplex_mtd.astype("category")
            if "antibody" not in demult_metadatacols:
                demult_metadatacols.insert(0, 'antibody')
            demultiplex_mtd = demultiplex_mtd[demult_metadatacols]
            # need to check the df merge column is compatible, so turn it into a category
            demultiplex_df['antibody']=demultiplex_df['antibody'].astype('category')
            demultiplex_df = demultiplex_df.merge(demultiplex_mtd, how="left", on="antibody", left_index=True).set_index("barcode_id")
            adatas[i].obs = adatas[i].obs.merge(demultiplex_df, how="left", left_index=True, right_index=True)

if len(adatas) != 1:
    print("concatenating ...")
    adata = adatas[0].concatenate(adatas[1:], batch_key="sample_id", batch_categories=caf.sample_id)
else:
    adata = adatas

# remove for mem issues
del adatas

# move sample_id to the front
cols = adata.obs.columns.tolist()
cols.insert(0, cols.pop(cols.index('sample_id')))
adata.obs = adata.obs.reindex(columns=cols)

print("merge in the scrublet scores")
# load the scrublet scores into the anndata (if they have been run)
if args.scrubletdir is not None:
    scrub_dir = args.scrubletdir
    # scrub_dir="./data/run_qc/scrublet"
    doubletscores = [pd.read_csv(scrub_dir + '/' + ss + '_scrublet_scores.txt', sep="\t", header=0) for ss in caf.sample_id]
    doubletscores = pd.concat(doubletscores, keys=caf.sample_id).reset_index(level="sample_id")
    # rename the barcodes to match up qwith what the adata.obs barcodes are
    doubletscores['barcode'] = doubletscores['barcode'] + '-' + doubletscores['sample_id']
    doubletscores = doubletscores.set_index('barcode').drop('sample_id', axis=1)

    # merge with adata.obs
    adata.obs = adata.obs.merge(doubletscores, how="left", left_index=True, right_index=True)

# qc ---------------------
# remove the basic filtering and put it in the integration space
# here just add the QC to the data
# MT proportion
# RP proportion
# Cell Cycle score
# IG genes proportion
# HBB genes proportion




IGgenes = pd.read_csv(args.iggenesfile, sep='\t')

IGgenes = IGgenes['gene_name'].tolist()
IGgenes = [ t for t in IGgenes if t.startswith('IG') ]
ccgenes = pd.read_csv(args.ccgenes, sep='\t')
sgenes = ccgenes[ccgenes["cc_phase"] == "s"]["gene_name"].tolist()
g2mgenes = ccgenes[ccgenes["cc_phase"] == "g2m"]["gene_name"].tolist()
HBgenes = ["HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ"]


adata.var['hb'] = [x in HBgenes for x in adata.var_names] # annotate the group of hb genes as 'hb'
adata.var['ig'] = [x in IGgenes for x in adata.var_names] # annotate the group of ig genes as 'ig'
adata.var['rp'] = adata.var_names.str.startswith('RP') # annotate the group of ribosomal genes as 'rp'
adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt','rp','hb','ig'], percent_top=None,log1p=False, inplace=True)
sc.tl.score_genes_cell_cycle(adata, s_genes=sgenes, g2m_genes=g2mgenes)


print("calculated scores and metrics")
if args.scrubletdir is not None:

    fxd = '_scrublet' + args.figure_suffix
    sc.pl.violin(adata, ['doublet_scores'],
                       multipanel=True,rotation=90,
                       jitter=0.4, show=False, groupby="sample_id",
                       save = fxd )

    fxd = '_ngenes_doubletscores' + args.figure_suffix 
    sc.pl.scatter(adata, x='n_genes_by_counts', y='doublet_scores', show=False, save = fxd)
    

fxd = "_metrics" + args.figure_suffix
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt','pct_counts_rp', 'pct_counts_hb','pct_counts_ig', 'S_score', 'G2M_score'],
                    multipanel=True, jitter=0.4, show=False, groupby="sample_id", rotation=90, save=fxd)

fxd = "_mt_counts" + args.figure_suffix
sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', show=False, save = fxd)

fxd = "_mt_ngenes" + args.figure_suffix
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', show=False, save = fxd)


print("saving anndata and obs in a metadata tsv file")
metafile = adata.obs
metafile["cellbarcode"] = adata.obs.index

savename = args.sampleprefix + "_cell_metadata.tsv"
metafile.to_csv(savename, sep='\t', index=True)  #consider changing to index false

adata.write(args.outfile)

print("done")




