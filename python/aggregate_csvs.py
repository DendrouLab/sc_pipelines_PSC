import argparse
import pandas as pd
import re
from itertools import chain
import os


parser = argparse.ArgumentParser()
parser.add_argument('--input_files_str',
                    default='',
                    help='')
parser.add_argument('--output_file',
                    default='',
                    help='')
parser.add_argument('--sample_prefix',
                    default='',
                    help='')
parser.add_argument('--clusters_or_markers',
                    default='',
                    help='')
parser.set_defaults(verbose=True)
args = parser.parse_args()


def is_float_try(string):
    try:
        float(string)
        return True
    except ValueError:
        return False


def splitall(path):
    allparts = []
    while 1:
        parts = os.path.split(path)
        if parts[0] == path:  # sentinel for absolute paths
            allparts.insert(0, parts[0])
            break
        elif parts[1] == path: # sentinel for relative paths
            allparts.insert(0, parts[1])
            break
        else:
            path = parts[0]
            allparts.insert(0, parts[1])
    return allparts


def extract_parameter_from_fname(fname, parameter, prefix):
    """
    Extract parameters from filename under the following circumstances:
    INPUT fname: res0.6_cluster.txt.gz or <parametervalue>_<analysis_type><file suffix>
    INPUT parameters: "res"
    OUTPUT: "0.6"
    Currently it doesn't work if the paramtervalue combo is next to the file suffix. (e.g. res0.6.txt.gz)
    """
    # first remove the sample prefix from the fname
    # so that any character strings that match the parameters are not picked up instead
    fname = re.sub(prefix, "", fname)
    # split by underscores
    fname_split0 = [re.split('_', x) for x in splitall(fname)]  # split by path and _
    fname_split1 = list(chain.from_iterable(fname_split0))
    # extract the string that matches the requested parameter
    fname_grep = list(filter(lambda x: parameter in x, fname_split1))[0]
    # extract the specific value
    value = fname_grep.replace(parameter, "")
    # if the value is a number, convert to float or integer as appropriate
    if is_float_try(value):
        value = float(value)
        # convert integer values to integers e.g. 1.0 to 1 to match seurat columns
        if int(value) == float(value):
            value = int(value)
    return value


infiles = custom_markers = re.split(',', args.input_files_str)
if args.clusters_or_markers == "clusters":
    combined_csv = pd.concat([pd.read_csv(f, sep='\t', index_col=0) for f in infiles], axis=1)
    # get colnames
    cnames = []
    for f in infiles:
        alg = extract_parameter_from_fname(f, 'alg', prefix=args.sample_prefix)
        res = extract_parameter_from_fname(f, 'res', prefix=args.sample_prefix)
        cnames.append(alg + '_res_' + str(res))
    combined_csv.to_csv(args.output_file, sep='\t', header=cnames, index=True)


if args.clusters_or_markers == "markers":
    li = []
    top_markers_file = re.sub("_top", "_all", args.output_file)
    excel_file = re.sub("_top.txt.gz", "_top.xlsx", args.output_file)
    with pd.ExcelWriter(excel_file) as writer:
        for ff in infiles:
            # print(os.path.join(in_path, ff))
            df = pd.read_csv(ff, sep='\t')
            # add a cluster column
            clust_val = extract_parameter_from_fname(ff, "cluster", prefix=args.sample_prefix)
            sname = "cluster" + str(clust_val)
            df['cluster'] = clust_val
            li.append(df)
            df.to_excel(writer, sheet_name=sname, index=False)
        frame = pd.concat(li, axis=0, ignore_index=True)
        frame.to_csv(args.output_file, sep='\t', header=True, index=False)
        frame_sub = frame[frame['p.adj.bonferroni'] < 0.05]
        frame_sub = frame_sub[frame_sub['avg_logFC'] > 0]
        frame_sub.to_csv(top_markers_file, sep='\t', header=True, index=False)

