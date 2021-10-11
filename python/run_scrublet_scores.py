'''
run scrublet on single channel
expects sample id and path to input data
'''
import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import argparse
import sys
import scanpy as sc

parser = argparse.ArgumentParser()
parser.add_argument("--sample_id",
                    default="sampleID",
                    help="name of the sample, usually just he name of the cellranger filtered folder")
parser.add_argument("--inputpath",
                    default="/path/to/cellranger/filtered_folder",
                    help="path to the single channel worth of cells for the doublet estimation to be run on")
parser.add_argument("--filetype",
                    default="cellranger",
                    help="type of data, options = ['cellranger', 'csv_matrix','txt_matrix', 'h5ad', 'h5']")
parser.add_argument("--outdir",
                    default="/gpfs3/well/combat/projects/preprocess/citeseq_final/FC_annotation/dotplot_minimal/",
                    help="string, b or t cells?")
parser.add_argument("--expected_doublet_rate",
                    default=0.06,
                    help="the expected fraction of transcriptomes that are doublets, typically 0.05-0.1. Results are not particularly sensitive to this parameter")
parser.add_argument("--sim_doublet_ratio",
                    default=2,
                    help="the number of doublets to simulate, relative to the number of observed transcriptomes. Setting too high is computationally expensive. Min tested 0.5")
parser.add_argument("--n_neighbors",
                    default=20,
                    help="Number of neighbors used to construct the KNN classifier of observed transcriptomes and simulated doublets. The default value of round(0.5*sqrt(n_cells)) generally works well.")
parser.add_argument("--min_counts",
                    default=2,
                    help="Used for gene filtering prior to PCA. Genes expressed at fewer than `min_counts` in fewer than `min_cells` (see below) are excluded")
parser.add_argument("--min_cells",
                    default=3,
                    help="Used for gene filtering prior to PCA. Genes expressed at fewer than `min_counts` (see above) in fewer than `min_cells` are excluded.")
parser.add_argument("--min_gene_variability_pctl",
                    default=85,
                    help="Used for gene filtering prior to PCA. Keep the most highly variable genes (in the top min_gene_variability_pctl percentile), as measured by the v-statistic [Klein et al., Cell 2015]")
parser.add_argument("--n_prin_comps",
                    default=30,
                    help="Number of principal components used to embed the transcriptomes prior to k-nearest-neighbor graph construction")
parser.add_argument("--use_thr",
                    default=True,
                    help="use a user defined thr to define min doublet score to split true from false doublets? if false just use what the software produces")
parser.add_argument("--call_doublets_thr",
                    default=0.25,
                    help="if use_thr is True, this thr will be used to define doublets")
args = parser.parse_args()

if args.filetype == "cellranger":
    if os.path.exists(args.inputpath):
        input_dir = args.inputpath
    else:
        sys.exit("the cells you're trying to load don't exist!")
    counts_matrix = scipy.io.mmread(input_dir + '/matrix.mtx.gz').T.tocsc()
    genes = np.array(scr.load_genes(input_dir + '/features.tsv', delimiter='\t', column=1))
    cellnames = pd.read_csv(input_dir + '/barcodes.tsv.gz', sep='\t', header=None)
elif args.filetype == "h5ad":
    if os.path.exists(args.inputpath):
        input_path = args.inputpath
    else:
        sys.exit("the cells you're trying to load don't exist!")
    adata = sc.read_h5ad(input_path)
    counts_matrix = adata.X.copy()
    genes = list(adata.var_names)
    cellnames = list(adata.obs_names)
elif args.filetype == "csv_matrix":
    if os.path.exists(args.inputpath):
        input_path = args.inputpath
    else:
        sys.exit("the cells you're trying to load don't exist!")
    counts_matrix = pd.read_csv(input_path, sep=',',index_col=0)
    genes = list(counts_matrix.columns)
    cellnames = list(counts_matrix.index)
    counts_matrix = counts_matrix.to_numpy()
elif args.filetype == "txt_matrix":
    if os.path.exists(args.inputpath):
        input_path = args.inputpath
    else:
        sys.exit("the cells you're trying to load don't exist!")
    counts_matrix = pd.read_csv(input_path, sep='\t',index_col=0)
    genes = list(counts_matrix.columns)
    cellnames = list(counts_matrix.index)
    counts_matrix = counts_matrix.to_numpy()
elif args.filetype == "h5":
    if os.path.exists(args.inputpath):
        input_path = args.inputpath
    else:
        sys.exit("the cells you're trying to load don't exist!")
    adata = sc.read_hdf(input_path)
    counts_matrix = adata.X.copy()
    genes = list(adata.var_names)
    cellnames = list(adata.obs_names)
    

print('Counts matrix shape: {} rows, {} columns'.format(counts_matrix.shape[0], counts_matrix.shape[1]))
print('Number of genes in gene list: {}'.format(len(genes)))

print("now initializing the scrublet object with expected_doublet_rate:\n")
print(args.expected_doublet_rate)

scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=float(args.expected_doublet_rate))
print("predicting doublets with params: \nmincells: %s \nmingenes: %s \nmin_gene_variabilty_pctl: %s\nn_prin_comps: %s\n" % (args.min_counts, args.min_cells, args.min_gene_variability_pctl, args.n_prin_comps))
print("predicting using:\n")

doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=int(args.min_counts),
                                                          min_cells=int(args.min_cells),
                                                          min_gene_variability_pctl=float(args.min_gene_variability_pctl),
                                                          n_prin_comps=int(args.n_prin_comps))

# cgat pipelines will probably parse a string
if args.use_thr == "True":
    use_thr = True

if use_thr:
    print("using default threshold to call doublets, instead of predicted %s, using default set to %s" %(round(scrub.threshold_, 2), args.call_doublets_thr))
    predicted_doublets = scrub.call_doublets(threshold=float(args.call_doublets_thr))
    
print("saving plots to outdir")

fig = scrub.plot_histogram()
fig[0].savefig(args.outdir + '/' + args.sample_id + '_' + "doubletScore_histogram.png",
               bbox_inches='tight', dpi=120)

print("cellnames and doublet scores and prediction")

data = pd.DataFrame({'doublet_scores': doublet_scores,
                     'predicted_doublets': predicted_doublets})

data['barcode'] = cellnames
data.to_csv(os.path.join(args.outdir + "/" + args.sample_id + "_scrublet_scores.txt"),
            sep="\t", index=False)
print('done')



