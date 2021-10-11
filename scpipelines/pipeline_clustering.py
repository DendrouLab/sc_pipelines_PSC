"""
CGAT pipeline for clustring single cell data with Scanpy.
# ASSUMED INPUT: 2 anndata objects, one containing all the data logn normalised but unscaled.
# the second containing scaled data and subset by highly variable genes.
# dimension reduction such as PCA or equivalent from harmobny or scanorama
# must have has been computed and saved within object

# This pipeline is designed to follow on from pipeline_integration.py

"""

from ruffus import *
import sys
import os
from cgatcore import pipeline as P
import pandas as pd
import cgatcore.iotools as IOTools
import re
from itertools import chain
import glob

PARAMS = P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
     "../pipeline.yml",
     "pipeline.yml"])

PARAMS['py_path'] =  os.path.join(os.path.dirname(os.path.dirname(__file__)), 'python')
PARAMS['r_path'] = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'R')


# ------------------------------------
# Generic functions
# ------------------------------------
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
    fname_split0 = [re.split('_',x) for x in splitall(fname)]  # split by path and _
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


# ------------------------------------
# Find neighbors
# ------------------------------------
def gen_neighbor_jobs():
    """
    Generate find neighbor jobs with all parameter combinations.
    """
    # same infile for all jobs
    # define files based on jobs
    infile = PARAMS['scaled_obj']
    batch_correct_str = str(PARAMS["batch_correction"])
    batch_correct = batch_correct_str.strip().replace(" ", "").split(",")
    nbrs_str = str(PARAMS["neighbors_n"])
    nbrs = nbrs_str.strip().replace(" ", "").split(",")
    metrics_str = str(PARAMS['neighbors_metric'])
    metrics = metrics_str.strip().replace(" ", "").split(",")
    pcs_str = str(PARAMS['neighbors_pcs'])
    pcs = pcs_str.strip().replace(" ", "").split(",")
    for mm in metrics:
        for bc in batch_correct:
            for pc in pcs:
                for rr in nbrs:
                    outfile = os.path.join(PARAMS['sample_prefix'] + "_nneigh" + str(rr) + "_pcs" + str(pc) + "_dir",
                                           PARAMS['sample_prefix'] + "_bc" + bc + "_met" + mm + "_neighbors.h5ad")
                    yield [infile, outfile]


# @follows(create_anndata)
@files(gen_neighbor_jobs)
def calc_neighbors(infile, outfile):
    print(outfile)
    dr = extract_parameter_from_fname(fname=outfile, parameter='bc', prefix=PARAMS['sample_prefix'])
    nn = extract_parameter_from_fname(fname=outfile, parameter='nneigh', prefix=PARAMS['sample_prefix'])
    met = extract_parameter_from_fname(fname=outfile, parameter="met", prefix=PARAMS['sample_prefix'])
    pc = extract_parameter_from_fname(fname=outfile, parameter="pcs", prefix=PARAMS['sample_prefix'])
    logfile = os.path.join(os.path.dirname(outfile),str(nn) + "_" + str(pc) + "_" + met + "_" + "calcneighbours.log")
    cmd = "python %(py_path)s/run_find_neighbors.py --infile %(infile)s \
    --outfile %(outfile)s \
    --batch_correction %(dr)s \
    --neighbors %(nn)s \
    --pcs %(pc)s \
    --metric %(met)s > %(logfile)s"
    P.run(cmd, job_threads=PARAMS["resources_threads_high"])


# ------------------------------------
# UMAP
# ------------------------------------
@subdivide(calc_neighbors, regex(r'(.*)/(.*)_(.*)_(.*)_neighbors.h5ad'), r"\1/\2*_umap.sentinel", "_umap.sentinel")
def gen_umap_jobs(infile, outfiles, outfile_suffix):
    """
    Generate cluster jobs with all parameter combinations.
    """
    # prefix = re.sub('_neighbors.h5ad', '', infile)
    prefix = os.path.split(infile)[0]
    mindist_str = str(PARAMS["umap_mindist"])
    mindist_vals = mindist_str.strip().replace(" ", "").split(",")
    for md in mindist_vals:
        output_file = os.path.join(prefix, PARAMS['sample_prefix'] + '_md' + str(md) + outfile_suffix)
        # yield [infile, output_file]
        IOTools.touch_file(output_file)


@transform(gen_umap_jobs, regex(r"(.*).sentinel"), r"\1.txt.gz")
def calc_umap(infile, outfile):
    print(infile)
    print(outfile)
    prefix = os.path.split(infile)[0]
    neighbour_file = glob.glob(os.path.dirname(infile) + '/*.h5ad')[0]
    md = extract_parameter_from_fname(infile, 'md', prefix=PARAMS['sample_prefix'])
    print(md)
    print(neighbour_file)
    logfile = prefix + "umapcalc.log" 
    cmd = "python %(py_path)s/run_umap.py \
            --infile %(neighbour_file)s \
            --outfile %(outfile)s \
            --min_dist %(md)s > %(logfile)s"
    P.run(cmd, job_threads=PARAMS["resources_threads_high"])


# ------------------------------------
# Clustering
# ------------------------------------
@subdivide(calc_neighbors,
           regex(r'(.*)/(.*)_(.*)_(.*)_neighbors.h5ad'),
           r"\1/\2_*/*_clusters.sentinel",
           "_clusters.sentinel")
def gen_cluster_jobs(infile, outfiles, outfile_suffix):
    """
    Generate cluster jobs with all parameter combinations.
    """
    # same infile for all jobs
    # infile = PARAMS['scaled_obj']
    # define files based on parameters
    infile_path = os.path.split(infile)[0]
    print(infile_path)
    sample_prefix = PARAMS['sample_prefix']
    resolutions_str = str(PARAMS["clusterspecs_cluster_resolutions"])
    resolutions = resolutions_str.strip().replace(" ", "").split(",")
    algorithm_str = str(PARAMS["clusterspecs_algorithm"])
    algorithms = algorithm_str.strip().replace(" ", "").split(",")

    for alg in algorithms:
        # for nn in neighbors_vals:
        for rr in resolutions:
            output_folder = os.path.join(infile_path, sample_prefix + "_alg" + alg + "_res" + str(rr))
            output_file = os.path.join(output_folder, "alg" + alg + "_clusters.sentinel")
            # yield [infile, output_file]
            if not os.path.exists(output_folder):
                os.makedirs(output_folder)
            IOTools.touch_file(output_file)


@transform(gen_cluster_jobs, regex(r"(.*)/(.*)/(.*)_clusters.sentinel"),
           r"\1/\2/\3_clusters.txt.gz")
def calc_cluster(infile, outfile):
    # pull resolution from the gen_cluster_jobs file names
    sentinel_fname = infile
    print(sentinel_fname)
    res = extract_parameter_from_fname(sentinel_fname, 'res', prefix=PARAMS['sample_prefix'])
    alg = extract_parameter_from_fname(sentinel_fname, 'alg', prefix=PARAMS['sample_prefix'])
    # select adata in file withe cooreect nneighbours
    nn = extract_parameter_from_fname(sentinel_fname, 'nneigh', prefix=PARAMS['sample_prefix'])
    pcs = extract_parameter_from_fname(sentinel_fname, 'pcs', prefix=PARAMS['sample_prefix'])
    # get the correct neighbour file
    nneigh_dir = os.path.dirname(os.path.dirname(sentinel_fname))
    nneigh_file = glob.glob(nneigh_dir + "/*.h5ad")[0]
    logfile = os.path.join(nneigh_dir, "calcclustering.log")
    cmd = "python %(py_path)s/run_clustering.py \
            --infile %(nneigh_file)s \
            --outfile %(outfile)s \
            --resolution %(res)s \
            --algorithm %(alg)s > %(logfile)s"
    P.run(cmd, job_threads=PARAMS["resources_threads_medium"])


# @follows(mkdir("test"))

@collate(calc_cluster,
         regex("(.*)/(.*)/(.*)clusters.txt.gz"),
         r"\1/all_res_clusters_list.txt.gz")
def aggregate_clusters(infiles, outfile):
    print(infiles)
    print(outfile)
    infiles_str = ','.join(infiles)
    cmd = "python %(py_path)s/aggregate_csvs.py \
               --input_files_str %(infiles_str)s \
               --output_file %(outfile)s \
               --sample_prefix %(sample_prefix)s \
               --clusters_or_markers clusters"
    P.run(cmd, job_threads=PARAMS["resources_threads_low"])


@follows(aggregate_clusters, calc_umap)
@transform(calc_neighbors, regex("(.*)/(.*).h5ad"), r'logs/\1_plot_clusters_umaps.sentinel')
def plot_cluster_umaps(infile, outfile):
    # print(infile)
    # print(outfile)
    # get associated umap
    nneigh_dir = os.path.dirname(infile)
    umaps = glob.glob(nneigh_dir + "/*_umap.txt.gz")
    print(umaps)
    clusters = os.path.join(nneigh_dir, "all_res_clusters_list.txt.gz")
    print(clusters)
    for uf in umaps:
        md = extract_parameter_from_fname(uf, "md", PARAMS['sample_prefix'])
        fig_dir = os.path.join(nneigh_dir, "figures")
        fig_suffix = "_all_res_md" + str(md) + ".png"
        umap_file = uf
        cmd = """python %(py_path)s/plot_cluster_umaps.py \
        --adata_object %(infile)s \
        --clusters_file %(clusters)s \
        --umap_file %(uf)s \
        --figure_dir %(fig_dir)s \
        --figure_suffix %(fig_suffix)s 
        """
        # cmd += " > %(outfile)s 2>&1"
        P.run(cmd, job_threads=PARAMS['resources_threads_low'], jobs_limit=1)

# all the defs
# def clustering_umbrella
@follows(plot_cluster_umaps)
@originate("cluster_analysis.sentinel")
def cluster_analysis(fname):
    IOTools.touch_file(fname)
    pass

# ------------------------------------
# Markers
# ------------------------------------
# The function below doesn't use the outfile parameter, but ruffus uses it to
# track files, the function uses outfile_suffix to build the filename
# Generate sentinel files in order to parallise the findMarkers step.
# I kind of hate having to use the sentinel files

#
@subdivide(calc_cluster,
           regex("(.*)/(.*)/(.*)clusters.txt.gz"),
           r"\1/\2/markers/cluster*_markers.sentinel",
           r"\1/\2/markers", "_markers.sentinel"
           )
def gen_marker_jobs(infile, outfiles, outfile_path,  outfile_suffix):
    """
    Generate cluster jobs for each resolution and cluster combo.
    Note that this might have to go into a separate script, in order to submit to server!
    depends how big the pickle files get.
    """
    print(infile)
    print(outfiles)
    # ignore unwanted input
    outfile_prefix = re.sub(".txt.gz", "", infile)
    # extract resolution from filename
    # load clusters
    df = pd.read_csv(infile, sep="\t")
    # how many clusters are we doing this for? get unique values in second column
    uniq_clust = list(set(df['clusters']))
    # generate sentinel files which will get used as input for FindMarkers
    if not os.path.exists(outfile_path):
        os.makedirs(outfile_path)
    for uc in uniq_clust:
        out = os.path.join(outfile_path, "cluster" + str(uc) + outfile_suffix)
        IOTools.touch_file(out)
#
#
# # for some reason this is submitting two jobs to the sever for each regex match.
# # I think it has semething to do with the regex in the previous call.
@active_if(PARAMS['resources_fewer_jobs']==False)
@transform(gen_marker_jobs, suffix("_markers.sentinel"),
           add_inputs(PARAMS['full_obj'], calc_cluster),
           "_markers.txt.gz")
def find_markers(infiles, outfile):
    """
    Runs scanpy.tl.rank_gene_groups in parallel for each cluster
    """
    print(infiles)
    print(outfile[0])
    # pull resolution from the gen_cluster_jobs file names
    # get seurat object
    sentinel_fname = os.path.basename(infiles[0])
    cluster_dname = os.path.dirname(os.path.dirname(infiles[0]))

    # log_file = re.sub(".sentinel", ".log", sentinel_fname)
    adata_infile = infiles[1]
    # # # get parameters
    cluster_val = str(extract_parameter_from_fname(sentinel_fname, 'cluster', prefix=PARAMS['sample_prefix']))
    cluster_file = os.path.join(cluster_dname, "alg" + PARAMS["clusterspecs_algorithm"] + "_clusters.txt.gz")
    cmd = """python %(py_path)s/run_find_markers.py \
    --adata_object %(adata_infile)s \
    --outputfile %(outfile)s \
    --cluster_file %(cluster_file)s \
    --cluster_value %(cluster_val)s \
    --pseudo_seurat %(markerspecs_pseudo_seurat)s \
    --minpct %(markerspecs_minpct)s \
    --threshuse %(markerspecs_threshuse)s \
    --mincells %(markerspecs_mincells)s
    """
    P.run(cmd, job_threads=PARAMS["resources_threads_high"])

# @collate([find_markers,find_markers_multi],
#          regex(r"(.*)/(.*)/markers/(.*).txt.gz"),
#          r"\1/\2/markers/markers_list.txt",
#          r"\1/\2/markers/")


@active_if(PARAMS['resources_fewer_jobs']==True)
# @follows(mkdir("out"))
@transform(calc_cluster,
           regex("(.*)/(.*)/(.*)clusters.txt.gz"),
           r"\1/\2/markers/cluster0_markers.txt.gz",
           r"\1/\2/markers/cluster", "_markers.txt.gz"
           )
def find_markers_multi(infile, outfile, out_path, suffix):
    """
    Runs scanpy.tl.rank_gene_groups in parallel for each cluster
    """
    print(infile)
    print(outfile)
    print(out_path)
    print(suffix)
    cmd = """python %(py_path)s/run_find_markers_multi.py \
    --adata_object %(full_obj)s \
    --output_prefix %(out_path)s \
    --output_suffix %(suffix)s \
    --cluster_file %(infile)s \
    --pseudo_seurat False \
    --minpct %(markerspecs_minpct)s \
    --threshuse %(markerspecs_threshuse)s \
    --mincells %(markerspecs_mincells)s
    """
    P.run(cmd, job_threads=PARAMS["resources_threads_high"])


# AIM make the regex flow be the same for find markers and find markers multi
# The problem with this is that "transform will make lots of jobs"
# how to make it pass mutliple jobs through the regex but only run one?
# 1. we could cheat and put the aggregate step in the find_markers_multi step
# 2. we could just remove find_markers_multi.
# 3. OR we could get aggregate markers to work it out by itself.

# @follows(mkdir("test"))
# @jobs_limit(10)
@collate([find_markers, find_markers_multi],
         regex(r"(.*)/(.*)/markers/cluster(.*).txt.gz"),
         r"\1/\2/markers/markers_list_top.txt.gz",
         r"\1/\2/markers/")
def aggregate_markers(infiles, outfile, in_path):
    # if using find_markers_multi
    if len(infiles) == 1:
        top_path=top_path=os.path.dirname(infiles[0])
        marker_files = glob.glob(top_path + "/cluster*.txt.gz")
        infiles_str = ','.join(marker_files)
    # else using find markers
    else:
        marker_files = infiles
        infiles_str = ','.join(marker_files)
    print(infiles_str)
    cmd = "python %(py_path)s/aggregate_csvs.py \
                  --input_files_str %(infiles_str)s \
                  --output_file %(outfile)s \
                  --sample_prefix %(sample_prefix)s \
                  --clusters_or_markers markers"
    P.run(cmd, job_threads=PARAMS["resources_threads_low"])


@originate(PARAMS['sample_prefix'] + "_cell_mtd.txt")
def write_metadata(outfile):
    cmd = """
    python %(py_path)s/write_metadata.py \
    --infile %(full_obj)s \
    --outfile %(outfile)s
    """
    P.run(cmd, job_threads=PARAMS['resources_threads_low'])


@follows(aggregate_markers)
@originate("marker_analysis.sentinel")
def marker_analysis(fname):
    IOTools.touch_file(fname)
    pass


# ------------------------------------
# Marker/cluster visualisation
# ------------------------------------

# limit jobs because reading the h5 file simultaneouly is problematic
@jobs_limit(1)
@follows(write_metadata)
@transform(aggregate_markers,
           regex(r"(.*)/(.*)/markers/markers_list_top.txt.gz"),
           add_inputs(output_from(write_metadata)),
           r"\1/\2/figures/\2_heatmaps.png", r"\1/\2/")
def heatmap_markers(infiles, outfile, data_path):
    """
    Plots a standard Seurat Heatmap, or a Complex Heatmap
    """
    # get cluster data
    fig_path = os.path.dirname(infiles[0])
    if not os.path.exists(fig_path):
        os.makedirs(fig_path)
    marker_file = infiles[0]
    metadata_file = infiles[1]
    cluster_file = glob.glob(data_path + "*clusters.txt.gz")[0]
    logfile = os.path.dirname(outfile) + "/marker_heatmap.log"
    cmd = """
    Rscript %(r_path)s/marker_heatmaps.R \
    --anndata %(scaled_obj)s \
    --marker_file %(marker_file)s \
    --output_file %(outfile)s \
    --cluster_file %(cluster_file)s \
    --metadata_file %(metadata_file)s \
    --docomplex %(plotspecs_docomplex)s \
    --subgroup %(plotspecs_discrete_variables)s > %(logfile)s"""
    P.run(cmd, job_threads=PARAMS['resources_threads_low'])


# limit jobs because reading the h5 file simultaneouly is problematic
# @follows(heatmap_markers)
@jobs_limit(1)
@transform(aggregate_markers,
           regex(r"(.*)/(.*)/markers/markers_list_top.txt.gz"),
           r"\1/\2/figures/dotplot_top_markers.png", r"\1/\2/")
def dotplot_markers(infile, outfile, data_path):
    """
    Plots some additional marker plots
    """
    # check there is a figures directory
    fig_path = os.path.dirname(outfile)
    if not os.path.exists(fig_path):
        os.makedirs(fig_path)
    # get cluster data for this level
    cluster_file = glob.glob(data_path + "*clusters.txt.gz")[0]
    cmd = """
    python %(py_path)s/plot_scanpy_markers.py \
    --adata_object %(full_obj)s \
    --marker_file %(infile)s \
    --cluster_file %(cluster_file)s \
    --figure_prefix %(fig_path)s"""
    P.run(cmd, job_threads=PARAMS['resources_threads_low'])




# Plot cluster metrics
@follows(write_metadata)
@transform(calc_cluster,
    regex("(.*)/(.*)/(.*)clusters.txt.gz"),
    add_inputs(output_from(write_metadata)),
    r"\1/\2/figures/metrics_plot.sentinel", r"\1", r"\1/\2/figures/\2")
def plot_cluster_metrics(infile, outfile, neighbor_dir, out_prefix):
    print(neighbor_dir)
    # print(out_prefix)
    dimred = "umap"
    # find umaps in the neighbours dir
    clusterid = infile[0]
    mtd = infile[1]
    fig_path = os.path.dirname(outfile)
    if not os.path.exists(fig_path):
        os.makedirs(fig_path)

    umaps = glob.glob(neighbor_dir + "/*umap.txt.gz")
    # print(umaps)
    for coords in umaps:
        print(coords)
        coords_uniq = "md" + str(extract_parameter_from_fname(coords, "md", prefix=PARAMS['sample_prefix']))
        print(coords_uniq)
        logfile = fig_path + "/cluster_metrics.log"
        cmd = """
        Rscript %(r_path)s/plot_cluster_metrics.R \
            --mtd_object %(mtd)s \
            --coords %(coords)s \
            --coords_suffix %(coords_uniq)s \
            --clusterid %(clusterid)s \
            --continuous_variables %(plotspecs_continuous_variables)s \
            --discrete_variables %(plotspecs_discrete_variables)s \
            --dimred %(dimred)s \
            --outfile_prefix %(out_prefix)s > %(logfile)s
            """
        P.run(cmd, job_threads=PARAMS['resources_threads_low'])
    IOTools.touch_file(outfile)



@collate(calc_cluster,
        regex("(.*)/(.*)/(.*)clusters.txt.gz"),
         # regex("(.*)alg(.*)_(.*)_cluster.txt.gz"),
         r"\1/figures/clustree.pdf", r"\1")
def plot_clustree(infiles, outfile, prefix):
    # convert infiles to comma sep. string
    infiles_str = ','.join(infiles)
    prefix = re.sub('_dir', '', prefix)
    # call R
    cmd = "Rscript %(r_path)s/plotclustree.R \
        --infiles %(infiles_str)s  \
        --infile_prefix %(prefix)s \
        --outfile %(outfile)s"
    P.run(cmd, job_threads=PARAMS['resources_threads_low'])


@active_if(PARAMS['custom_markers_files'] is not None)
@transform(calc_cluster,
           regex("(.*)/(.*)/(.*)clusters.txt.gz"),
           r'\1/\2/figures/custom_markers.sentinel', r"\1/\2/figures/")
def plot_custom_markers(infile, outfile, fig_path):
    cluster_file = infile
    # add the full path to each of the custom marker files then recombine
    custom_markers = re.split(',', PARAMS['custom_markers_files'])
    custom_markers_list = [os.path.join(PARAMS['custom_markers_path'], cm) for cm in custom_markers]
    marker_file_str = ",".join(custom_markers_list)
    # run plotting
    cmd = """
    python %(py_path)s/plot_custom_markers.py \
    --adata_object %(full_obj)s \
    --adata_layer X
    --marker_files %(marker_file_str)s \
    --clusters_file %(cluster_file)s \
    --figure_dir %(fig_path)s
    """
    P.run(cmd, job_threads=PARAMS['resources_threads_low'])
    IOTools.touch_file(outfile)


@active_if(PARAMS['custom_markers_files'] is not None)
@transform(calc_umap, regex('(.*)/(.*)_umap.txt.gz'),
           r'\1/figures/umap_\2_markers.png',
           r'\1/figures/', r'_\2_markers.png')
def plot_custom_markers_umap(infile, outfile, fig_path, figure_suffix):
    print(infile)
    print(outfile)
    marker_file = os.path.join(PARAMS['custom_markers_path'], PARAMS['custom_markers_minimal'])
    logfile = os.path.dirname(outfile) + "/custom_markers.log"
    cmd = """
    python %(py_path)s/plot_custom_markers_umap.py \
    --adata_object %(full_obj)s \
    --adata_layer X
    --variables_file %(marker_file)s \
    --umap_coords %(infile)s \
    --figure_dir %(fig_path)s \
    --figure_suffix %(figure_suffix)s > %(logfile)s
    """
    P.run(cmd, job_threads=PARAMS['resources_threads_low'])


@follows(plot_clustree, plot_custom_markers,
         plot_custom_markers_umap, plot_cluster_umaps,
         heatmap_markers,
         plot_cluster_metrics, dotplot_markers)
@originate("visualisation.sentinel")
def visualisation(fname):
    IOTools.touch_file(fname)
    pass

#  add plot custom markers heatmap here


@follows(cluster_analysis, calc_umap, marker_analysis, visualisation)
def full():
    """
    All cgat pipelines should end with a full() function which updates,
    if needed, all branches of the pipeline.
    The @follows statement should ensure that all functions are covered,
    either directly or as prerequisites.
    """
    pass


def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
