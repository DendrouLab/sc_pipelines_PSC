from ruffus import *
import sys
import os
# import re
from cgatcore import pipeline as P
import pandas as pd
import cgatcore.iotools as IOTools

PARAMS = P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
     "../pipeline.yml",
     "pipeline.yml"])

unfilt_file = PARAMS['sample_prefix'] + '_unfilt.h5ad'
sprefix = PARAMS['sample_prefix']


@follows(mkdir("logs"))
@follows(mkdir("figures"))
@originate("logs/0_plot10x.log")
def plot_tenx_metrics(outfile):
    r_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), "R")
    cmd = """
        Rscript %(r_path)s/produce_barplot_10xmetric.v3.R
            --csvpaths %(submission_file)s
            --outdir ./
            --figdir ./figures/
            --project %(project)s
            --kneeplot %(kneeplot)s > %(outfile)s
            """
    P.run(cmd, job_threads=PARAMS['resources_threads_low'])



def gen_scrublet_jobs():
    """
    Generate a run_scrublet job for each line in submission.txt
    """
    # infile: cellranger output
    # outfile: args.outdir + "/" + args.sample_id +"_scrublet_scores.txt"
    caf = pd.read_csv(PARAMS['submission_file'], sep='\t')
    for nn in range(0, caf.shape[0]):
        infile = caf['filtered_path'][nn]
        outfile = "./scrublet/" + caf['sample_id'][nn] + "_scrublet_scores.txt"
        out_dir = "./scrublet"
        sample_id = caf['sample_id'][nn]
        yield infile, outfile, out_dir, sample_id



@follows(mkdir("logs"))
@follows(mkdir("figures"))
@files(gen_scrublet_jobs)
def run_scrublet(infolder, outfile, outdir, sample_id):
    py_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), "python")
    cmd = """
        gunzip -c %(infolder)s/features.tsv.gz > %(infolder)s/features.tsv;
        python %(py_path)s/run_scrublet_scores.py
        --sample_id %(sample_id)s
        --inputpath %(infolder)s
        --outdir %(outdir)s
        """
    if PARAMS['scr_expected_doublet_rate'] is not None:
        cmd += " --expected_doublet_rate %(scr_expected_doublet_rate)s"
    if PARAMS['scr_sim_doublet_ratio'] is not None:
        cmd += " --sim_doublet_ratio %(scr_sim_doublet_ratio)s"
    if PARAMS['scr_n_neighbours'] is not None:
        cmd += " --n_neighbors %(scr_n_neighbours)s"
    if PARAMS['scr_min_counts'] is not None:
        cmd += " --min_counts %(scr_min_counts)s"
    if PARAMS['scr_min_cells'] is not None:
        cmd += " --min_cells %(scr_min_cells)s"
    if PARAMS['scr_min_gene_variability_pctl'] is not None:
        cmd += " --min_gene_variability_pctl %(scr_min_gene_variability_pctl)s"
    if PARAMS['scr_n_prin_comps'] is not None:
        cmd += " --n_prin_comps %(scr_n_prin_comps)s"
    if PARAMS['scr_use_thr'] is not None:
        cmd += " --use_thr %(scr_use_thr)s"
    if PARAMS['scr_call_doublets_thr'] is not None:
        cmd += " --call_doublets_thr %(scr_call_doublets_thr)s"
    cmd += " > logs/1_scrublet_" + sample_id + ".log;"
    cmd += "rm %(infolder)s/features.tsv"
    P.run(cmd, job_threads=PARAMS['resources_threads_medium'])
    IOTools.touch_file(outfile)


orfile=sprefix + "_cell_metadata.tsv"
#@originate("sample_metadata.sentinel")


@follows(run_scrublet)
@originate(orfile)
def run_scanpy_qc(metadata):
    # infile = submission file
    py_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), "python")
    resources_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), "resources")
    anndata_out = PARAMS['sample_prefix'] + "_unfilt.h5ad"
    sampleprefix = PARAMS['sample_prefix']
    cmd = """
        python %(py_path)s/run_scanpyQC.py
          --submissionfile %(submission_file)s
          --sampleprefix %(sampleprefix)s
          --metadatacols  %(metadatacols)s
          --outfile %(anndata_out)s
          --figdir figures
          --scrubletdir scrublet
          """
    # if params is specified otherwise default values will be used.
    if PARAMS['iggenesfile'] is not None:
        if PARAMS['iggenesfile'] == "default":
            iggenesfile = resources_path + "/exclude_genes_HLAIGTR_v1.txt"
        else:
            iggenesfile = PARAMS['iggenesfile']
        cmd += " --iggenesfile %(iggenesfile)s"
    if PARAMS['ccgenes'] is not None:
        if PARAMS['ccgenes'] == "default":
            ccgenesfile = resources_path + "/cell_cycle_genes.tsv"
        else:
            ccgenesfile = PARAMS['ccgenes']
        cmd += " --ccgenes %(ccgenesfile)s"
    if PARAMS['demultiplex_include'] is not None:
        cmd += " --demultiplexing %(demultiplex_include)s"
        cmd += " --demultiplexing_metadatacols %(demultiplex_metadatacols)s"
    # add log file
    cmd += " > logs/2_qc.log;"
    cmd += " rm -r cache"
    P.run(cmd, job_threads=PARAMS['resources_threads_high'])
    #IOTools.touch_file(metadata)

# ----------
# plotting script
# ----------

@transform(run_scanpy_qc,
            regex(rf"{sprefix}_cell_metadata.tsv"),
            r"logs/3_plot.log")
def plot_qc(infile, log_file):
    R_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), "R")
    sampleprefix = PARAMS['sample_prefix']
    #metadatafile = sampleprefix + "_metadata.tsv"
    cmd = """
    Rscript %(R_path)s/plotQC.R 
    --prefilter TRUE
    --metadata %(infile)s 
    --sampleprefix %(sampleprefix)s
    --groupingvar %(plotqc_grouping_var)s
    --qcmetrics %(plotqc_metrics)s > %(log_file)s 2>&1
    """
    P.run(cmd, job_threads=PARAMS['resources_threads_medium'])


# ------------
@follows(plot_qc, plot_tenx_metrics)
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
