import pickle
import argparse
import sys
import scanpy as sc
sc.settings.verbosity = 4

parser = argparse.ArgumentParser()
parser.add_argument('--infile', default='',
                    help="file name, format: anndata_scaled.h5ad")
parser.add_argument('--outfile', default='anndata_wneighbors.h5ad',
                    help="file name, format: .h5ad")
parser.add_argument('--neighbors', default=30, help="no. neighbours parameters for sc.pp.neighbors()")
parser.add_argument('--batch_correction', default='pca', help="dimension reduction for sc.pp.neighbors(), default='pca', \
                                                        alt. 'harmony' or 'scanorama'")
parser.add_argument('--pcs', default=30, help="no. of pcs from dimension reduction to calucalte neighbours from ")
parser.add_argument('--metric', default='euclidean', help="dimension reduction for sc.pp.neighbors(), default='pca'")

args = parser.parse_args()
print("running: run_find_neighbours.py with args")
print(args)

# read data
print("reading adata")
adata = sc.read_h5ad(args.infile)
print("adata read in")
# run command

print("batch correction choice is %s" % args.batch_correction)
if args.batch_correction == 'None':
    print("running neighbours from pca, on data with no batch correction")
    sc.pp.neighbors(adata, n_neighbors=int(args.neighbors), n_pcs=int(args.pcs), metric=args.metric)
elif args.batch_correction == 'combat':
    print("running neighbours from pca on combat corrected data")
    sc.pp.neighbors(adata, n_neighbors=int(args.neighbors), n_pcs=int(args.pcs), metric=args.metric)
elif args.batch_correction == 'harmony':
    print("running neighbours from harmony corrected data")
    sc.pp.neighbors(adata, n_neighbors=int(args.neighbors), n_pcs=int(args.pcs), use_rep='X_harmony', metric=args.metric)
elif args.batch_correction == 'scanorama':
    print("running neighbours from scanorama corrected data")
    sc.pp.neighbors(adata, n_neighbors=int(args.neighbors), n_pcs=int(args.pcs), use_rep="X_scanorama", metric=args.metric)
elif args.batch_correction == 'bbknn':
    print("bbknn neighbours should already be included in the object, not rerunning.")
    try:
        adata.uns['neighbors']
    except KeyError:
        sys.exit("bbknn neighbours not found, make sure they are in the scaled object")
    pass
else:
    sys.exit("dimension_reduction not recognised")
# save with n_neighbours included in output file.


print("saving data")
adata.write(args.outfile)
