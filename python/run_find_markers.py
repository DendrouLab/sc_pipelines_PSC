import scanpy as sc
import argparse
import pandas as pd
import sys
import re
import numpy as np
from itertools import compress
from scipy.sparse import issparse

parser = argparse.ArgumentParser()
parser.add_argument("--adata_object",
                    default="/Users/crg/Documents/Projects/combat/data/scanpy_pipe_test/pbmc3k_nn30_drharmony_metcosine_neighbors.h5ad",
                    help="file name, format: .h5ad")
parser.add_argument("--outputfile",
                    default="/Users/crg/Documents/Projects/combat/data/scanpy_pipe_test/cluster1_markers.txt",
                    help="file name, format: .h5ad")
parser.add_argument("--cluster_file",
                    default="/Users/crg/Documents/Projects/combat/data/scanpy_pipe_test/pbmc3k_nn30_res0.2_algleiden_cluster.txt",
                    help="file name, format: .h5ad")
parser.add_argument("--cluster_value", type=str,
                    default='1', help="the ID of the cluster that we want to look for markers in.")
parser.add_argument("--testuse", default='wilcoxon', help="wilcoxon, t-test or logreg")
parser.add_argument("--pseudo_seurat", default='True', help="apply seurat like filtering before marker test")
parser.add_argument("--minpct", type=str, default='0.1', help="minimum fraction of cells expressing gene")
parser.add_argument("--mindiffpct", type=str, default='-inf',
                    help="minimum fraction difference between cluster and other cells. Setting not recommended.")
parser.add_argument("--threshuse",  type=str, default='0.25',
                    help="testing limited to genes with this (log scale) difference in mean expression level.")
parser.add_argument("--mincells", type=str, default='3',
                    help="minimum number of cells required (applies to cluster and to the other cells)")


args = parser.parse_args()


def exp_mean_sparse(x):
    """
    This returns the log of the mean of the not-logged data
    this version of the function works directly on the sparse matrix but is slow.
    """
    # convert out of compressed sparse matrix
    return np.log((np.sum(x.expm1())/x.shape[1]) + 1)


def exp_mean_dense(x):
    """
    This returns the log of the mean of the not-logged data
    this version of the function requires a dense matrix, so might be memory hungry?
    But it is super fast
    """
    # convert out of compressed sparse matrix
    return np.log((np.sum(np.exp(x)-1)/x.shape[1]) + 1)


def which_ind(bool_list):
    return list((compress(range(len(bool_list)), bool_list)))


def which_val(bool_list, val_list):
    return list((compress(val_list, bool_list)))


def pseudo_seurat(adata, arg_minpct, arg_mindiffpct, arg_logfcdiff):
    """
    alternative method that"s more like seurat (pseudo seurat if you will)
    In that you filter genes before running rank genes
    ---
    1.  define idents
    2. define cluster cells and other cells
    3. check min cells
    4. compute percentages and difference from equiv data slot to s.obj@assays$RNA@data
    5. computer mean expression levels and differences
    6. save stats
    7. define background based on min_pct (pct of cells that have to express marker)
    8. filter stats for testing and save universe
    genes|cluster_mean|other_mean|diff_mean|cluster_pct|other_pct|max_pct|min_pct|diff_pct|background
    9. Filter adata.raw based on these genes?

    9. Find markers based on cluster cells and other cells,
    values from adata.raw ??
    gene: only expressed
    cells: we need to define manually based on above stats.
    10. save results.
    """
    # define cells
    cluster_cells_ind = which_ind(adata.obs["idents"] == "1")
    other_cells_ind = which_ind(adata.obs["idents"] == "0")

    # compute perecentage expressed
    # from normnalised but not scaled data
    # remember cells are rows and genes are columns

    # note: I don't know why norm_counts[cluster_cell_ind:, col_ind] deosn"t work, but it doesn't
    cluster_pct = (adata.X[cluster_cells_ind, :] > 0).sum(axis=0) / len(cluster_cells_ind)
    other_pct = (adata.X[other_cells_ind, :] > 0).sum(axis=0) / len(other_cells_ind)

    pcts = pd.DataFrame(np.vstack((cluster_pct,  other_pct)).transpose())
    max_pct = pcts.max(axis=1)
    min_pct = pcts.min(axis=1)
    diff_pct = max_pct - min_pct
    take_diff_pct = diff_pct > arg_mindiffpct

    # remove genes that are not expressed higher than 0.1 in one of the groups
    take_min_pct = max_pct > arg_minpct


    # KEEP IN CASE NP.ARRAY METHOD USES TOO MUCH MEMORY
    # import time
    # this has the potential to be very slow. Transposeing it speeds it up a bit.
    # I need to understand sparse matrices better to make it work
    # start = time.time()
    # nct = adata.X.T[:,cluster_cells_ind]
    # cluster_mean0 = [exp_mean_sparse(nct[x,:]) for x in range(0,nct.shape[0])]
    # end = time.time()
    # print(end - start)
    #
    # start = time.time()
    # nct = adata.X.T[:, other_cells_ind]
    # other_mean0 = [exp_mean_sparse(nct[x,:]) for x in range(0, nct.shape[0])]
    # end = time.time()
    # print(end - start)

    # extract the counts for cluster cells and calculate exp means on each row
    nct = adata.X.T[:, cluster_cells_ind]
    cluster_mean = np.apply_along_axis(exp_mean_dense, 1, nct.todense())

    # likewise for non-cluster cells
    nct = adata.X.T[:, other_cells_ind]
    other_mean = np.apply_along_axis(exp_mean_dense, 1, nct.todense())
    diff_mean = abs(cluster_mean - other_mean)

    # remove genes with less than threshold difference
    take_thresh = diff_mean > arg_logfcdiff
    # take = if a cell passes all the tests then it is to be kept.
    take = [a and b and c for a, b, c in zip(take_thresh, take_min_pct, take_diff_pct)]
    print("saving universe for fisher test")
    stats_df = pd.DataFrame(np.vstack((adata.var_names, cluster_mean, other_mean, diff_mean,
                                       cluster_pct, other_pct, max_pct, min_pct, diff_pct, take)).transpose(),
                            columns=["gene", "cluster_mean", "other_mean", "diff_mean",
                                     "cluster_pct", "other_pct",
                                     "max_pct", "min_pct", "diff_pct", "background"])
    return stats_df


# script -------------------------------

# check for options

print("Running with options")
print(args)
print('\n')
# read data
adata = sc.read_h5ad(args.adata_object)
# adata = sc.read_h5ad("../../data/scanpy_pipe_test/pbmc3k_nn30_drharmony_metcosine_neighbors.h5ad")

# load clusters
df = pd.read_csv(args.cluster_file, sep="\t", index_col=0)
# df = pd.read_csv("../../data/scanpy_pipe_test/pbmc3k_res0.6_nn30_algleiden_cluster.txt", sep="\t", index_col=0 )
#
adata_shape = adata.shape


# add clusters to adata.obs
adata.obs = pd.merge(adata.obs, df, how="left", left_index=True, right_index=True )
# check check we have not lost cells in merge
if adata.shape != adata_shape:
    print("some cells lost in merge, not all cells have a cluster?")

# convert adata.obs.clusters to category type
# adata.obs["clusters_cat"] = adata.obs["clusters"].astype("category")
# work out mean expression etc.
# work out an appropriate background.
#

# how many unique clusters?

# 1. set idents
adata.obs["idents"] = ["1" if cc==int(args.cluster_value) else "0" for cc in adata.obs["clusters"]]
adata.obs["idents"] = adata.obs["idents"].astype("category")

# check we have enough cells
if args.mincells == "None":
    min_cells = 3
else:
    min_cells = int(args.mincells)

n_cluster = sum([cell is '1' for cell in adata.obs['idents']])
n_other = sum([cell is '0' for cell in adata.obs['idents']])
if n_cluster <= int(args.mincells) | n_other <= int(args.mincells):
    sys.exit("not enough cells in cluster, min cells: %i" % int(args.mincells))

print("Cluster %s number of cells: %i" % (args.cluster_value, n_cluster))
print("Other number of cells: %i\n" % (n_other))

# there is some weird difference between dense and sparse array?
if issparse(adata.X):
    bool_list = (adata.X.sum(axis=0) > 0).tolist()[0]
else:
    bool_list = (adata.X.sum(axis=0) > 0).tolist()
adata = adata[:, bool_list]

if isinstance(args.pseudo_seurat, str):
    if args.pseudo_seurat == 'True':
        run_pseudo_seurat = True
    elif args.pseudo_seurat == 'False':
        run_pseudo_seurat = False
    else:
        print(args.pseudo_seurat)
        sys.exit("spelling_error: please specify --pseudoseurat as True or False")
elif type(args.pseudo_seurat) == 'bool':
    run_pseudo_seurat = args.pseudo_seurat
else:
    sys.exit("type error %s: please specify --pseudo_seurat as True or False" % type(args.pseudo_seurat))


if run_pseudo_seurat:
    print("running pseudo seurat with params: \nminpct: %s \nmindiffpct: %s \nlogFCthreshold: %s\n\n" % (args.minpct, args.mindiffpct, args.threshuse))

    filter_stats = pseudo_seurat(adata,
                                 arg_minpct=float(args.minpct),
                                 arg_mindiffpct=float(args.mindiffpct),
                                 arg_logfcdiff=float(args.threshuse))
    print("number of genes remaining after filtering:  %i\n" % filter_stats['background'].sum())
    # adata = adata.raw.to_adata().copy()
    adata = adata[:, filter_stats['background'].tolist()]
#     adata.raw = adata
    print(adata)
    print('\n')
else:
    # get backgrounds
    pass
    # write out background after rank_gene_groups


# now we can actually run find markers,
# run on a subset of genes
sc.tl.rank_genes_groups(adata, groups=["1"], groupby="idents", reference="0",
                        method="wilcoxon", n_genes=float("inf"), corr_method="bonferroni")

markers = sc.get.rank_genes_groups_df(adata, group="1")

# merge markers with filter_stats_df, so that we can see all the scroes for gsea
#

# filter positive only and adjusted pval < 0.05
# markers = markers[(markers.logfoldchanges > 0) & (markers.pvals_adj < 0.05)]
if sc.__version__ < "1.7.2":
    markers.columns = ['scores', 'gene', 'avg_logFC', 'pvals', 'p.adj.bonferroni']
elif sc.__version__ == "1.7.2"
    markers.columns = ['gene', 'scores', 'avg_logFC', 'pvals', 'p.adj.bonferroni']

markers.head()
print("Number of Markers for cluster %s: %i" % (args.cluster_value, markers.shape[0]))
print("Rank gene groups complete")


print("saving markers to:" + args.outputfile)
markers.to_csv(args.outputfile, sep="\t", index=False)

if run_pseudo_seurat:
    fname = re.sub('_markers.txt', '_background.txt', args.outputfile)
    print("saving background to:" + fname)
    filter_stats.to_csv(fname, sep="\t", index=False)


# quick do some plotting!
