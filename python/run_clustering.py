import pickle
import scanpy as sc
import argparse
import pandas as pd
import sys

parser = argparse.ArgumentParser()
parser.add_argument("--infile", default="/Users/crg/Documents/Projects/combat/data/scanpy_pipe_test/sc_preprocess.h5ad", help="file name, format: .h5ad")
parser.add_argument("--outfile", default="/Users/crg/Documents/Projects/combat/data/scanpy_pipe_test/clusters.txt", help="file name, format: .h5ad")
parser.add_argument("--resolution", default=0.5, help="no. neighbours parameters for sc.pp.neighbors()")
parser.add_argument("--algorithm", default="leiden", help="algortihm choice from louvain and leiden")

args = parser.parse_args()
print(args)
# read data
adata = sc.read_h5ad(args.infile)

# run command
if args.algorithm == "louvain":
    sc.tl.louvain(adata, resolution=float(args.resolution), key_added='clusters')
elif args.algorithm == "leiden":
    sc.tl.leiden(adata, resolution=float(args.resolution), key_added='clusters')
else:
    sys.exit("algorithm not found: please specify 'louvain' or 'leiden'")

## write out clusters as text file
clusters = adata.obs['clusters']


# fname = str(args.sample_prefix) + str(args.neighbors) + ".h5ad"
# save compress_object
clusters = pd.DataFrame(adata.obs['clusters'])
clusters.to_csv(args.outfile, sep='\t')

