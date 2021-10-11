import argparse
import numpy as np
import pandas as pd
import harmonypy as hm
import seaborn as sns
import os
import re
from itertools import product

parser = argparse.ArgumentParser()
parser.add_argument('--input_files')
parser.add_argument('--batch_mtd')
parser.add_argument('--integration_col', default='batch')
parser.add_argument('--output_file')
parser.add_argument('--fig_dir')
args = parser.parse_args()

# load files
input_files = args.input_files.split(',')
umaps_list = [pd.read_csv(x, index_col=0) for x in input_files]

# get index_names
umaps_index = [re.sub('umap_bc_|.csv', '', os.path.basename(x)) for x in input_files]
# load batch info.
batch_df = pd.read_csv(args.batch_mtd, index_col=0)

# compute LISI for each batch
columns = batch_df.columns.to_list()
list_arrays = [hm.compute_lisi(x, batch_df, columns) for x in umaps_list]
# list_arrays = [hm.compute_lisi(x, batch_df, [args.integration_col]) for x in umaps_list]

# put LISI scores into a pandas df
lisi_df = pd.DataFrame(np.column_stack(list_arrays))
cnames = [(e1 + "_" + e2) for e1,e2 in product(umaps_index, columns)]

# lisi_df.columns = umaps_index
lisi_df.columns = cnames 

# save file
lisi_df.to_csv(args.output_file)

# make a density  plot of LISI scores.
plot_df = lisi_df.melt(var_name="Correction", value_name="LISI score")
sns_plot = sns.displot(plot_df, x='LISI score', hue='Correction', kind="kde")
sns_plot.savefig(os.path.join(args.fig_dir, "LISI_scores.png"))
