import argparse
import numpy as np
import pandas as pd
import seaborn as sns
import os
import re
import glob

parser = argparse.ArgumentParser()
parser.add_argument("--input_files")
parser.add_argument("--fig_dir")
parser.add_argument("--batch_mtd")
parser.add_argument("--qc_metrics")
parser.add_argument("--groupingvar")
parser.add_argument("--meta_df", default="", help="the full metadata, to save the full umap file at the end")
parser.add_argument('--integration_col', default='batch', help='the column in the adata.obs on which you have integrated')

args = parser.parse_args()

# load files
input_files = args.input_files.split(",")
# input_files = glob.glob("./data/run_scanpy_integration/batch_correction/umap*.csv")
batch_df = pd.read_csv(args.batch_mtd, index_col=0)
# batch_df = pd.read_csv("./data/run_scanpy_integration/batch_correction/batch_mtd.csv", index_col=0)
meta_df = pd.read_csv(args.meta_df, sep="\t", index_col=0) 
fullcolumns = meta_df.columns.to_list()
columns = batch_df.columns.to_list()
missing = np.setdiff1d(columns,fullcolumns)

print("reading in all umaps")
umaps_list = [pd.read_csv(x, index_col=0) for x in input_files]
# get index_names
umaps_index = [re.sub("umap_bc_|.csv", "", os.path.basename(x)) for x in input_files]
# load batch info.

# get dimensions to make umaps list and index it one pandas
nrow = batch_df.shape[0]
n = len(umaps_index)

# add batch info and batch correction methodindex
umaps_list = [np.column_stack([[umaps_index[x]]*nrow, umaps_list[x], meta_df, batch_df[missing.tolist()] ]) for x in range(0, n)]

umaps_df = pd.DataFrame(np.row_stack(umaps_list))

umaps_df.columns = ["method", "umap_1", "umap_2"]+ fullcolumns + missing.tolist()

# convert the umap coord back to numeric
umaps_df["umap_1"] = umaps_df["umap_1"].astype("float")
umaps_df["umap_2"] = umaps_df["umap_2"].astype("float")

# reorder methods so None is at the most left on plot
umaps_df["method"] = umaps_df["method"].astype("category")
# put categories in alphabetical order
umaps_index.sort()
# put none at the top of the list
umaps_index.insert(0, umaps_index.pop(umaps_index.index("none")))
umaps_df['method'] = umaps_df['method'].cat.reorder_categories(umaps_index, ordered=True)

print("plotting with seaborn")


grpvar = [x.replace (" " ,"") for x in args.groupingvar.split(",")]
grpvar = [value for value in grpvar if value in meta_df.columns] 
columns = list(set(columns).difference(set(grpvar))) + grpvar


for col in columns:
    print(col)
    # first plot just facetted on method
    g = sns.FacetGrid(umaps_df, col="method", col_wrap=3, sharex=False, sharey=False)
    g = (g.map(sns.scatterplot, "umap_1", "umap_2", col, s=5, alpha=.7))
    g.add_legend() 
    g.savefig(os.path.join(args.fig_dir, col + "_umap_bc_method.png"))
    # second plot facetted on method and batch
    print("second plot...")
    g = sns.FacetGrid(umaps_df, col="method", row=col, hue=col,sharex=False, sharey=False)
    g = (g.map(sns.scatterplot, "umap_1", "umap_2", s=5, alpha=.7))
    g.add_legend() 
    g.savefig(os.path.join(args.fig_dir,col + "_umap_bc_method_v_batch.png"))

qcmetrics = [x.replace (" " ,"") for x in args.qc_metrics.split(",")]

for col in qcmetrics:
    print(col)
    g = sns.FacetGrid(umaps_df, col="method", col_wrap=3, sharex=False, sharey=False)
    g = (g.map(sns.scatterplot, "umap_1", "umap_2", col, s=5, alpha=.7))
    g.add_legend() #how to add a human readable legend???
    g.savefig(os.path.join(args.fig_dir, col + "_umap_bc_method.png"))
    

# save umap df to file
print("save umap df to file")
umaps_df.to_csv("batch_correction/combined_umaps.tsv", sep="\t")
print("done")