from ruffus import *
import sys
import os
from cgatcore import pipeline as P
import pandas as pd
import cgatcore.iotools as IOTools
import re
from itertools import chain
import glob
# from itertools import chain
# import glob

# import pandas as pd

PARAMS = P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
     "../pipeline.yml",
     "pipeline.yml"])

filt_file = PARAMS['sample_prefix'] + '_filt.h5ad'
PARAMS['py_path'] = os.path.join(os.path.dirname(os.path.dirname(__file__)), "python")
PARAMS['R_path'] = os.path.join(os.path.dirname(os.path.dirname(__file__)), "R")
PARAMS['resources_path'] = os.path.join(os.path.dirname(os.path.dirname(__file__)), "resources")


@mkdir("logs")
@originate(filt_file)
def filter_anndata(outfile):
    if PARAMS['filtering_run']:
        cmd = """
        python %(py_path)s/run_filter.py
        --input_anndata %(unfiltered_obj)s
        --output_anndata %(outfile)s
        """
        if PARAMS['filtering_min_genes'] is not None:
            cmd += " --min_genes %(filtering_min_genes)s"
        if PARAMS['filtering_max_genes'] is not None:
            cmd += " --max_genes %(filtering_max_genes)s"
        if PARAMS['filtering_min_cells'] is not None:
            cmd += " --min_cells %(filtering_min_cells)s"
        if PARAMS['filtering_max_counts'] is not None:
            cmd += " --max_counts %(filtering_max_counts)s"
        if PARAMS['filtering_percent_mito'] is not None:
            cmd += " --percent_mito %(filtering_percent_mito)s"
        if PARAMS['filtering_percent_ribo'] is not None:
            cmd += " --percent_ribo %(filtering_percent_ribo)s"
        if PARAMS['filtering_percent_hb'] is not None:
            cmd += " --percent_hb %(filtering_percent_hb)s"
        if PARAMS['filtering_drop_nas_col'] is not None:
            cmd += " --drop_nas_col %(filtering_drop_nas_col)s"
        cmd += " > logs/0_filtering.log 2>&1"
        P.run(cmd,  job_threads=PARAMS['resources_threads_low'])
    else:
        try:
            f = open(outfile)
            f.close(outfile)
        except IOError:
            print("filter is not set to True, but there is no %s file" % outfile)
            sys.exit(1)


#@follows(filter_anndata)
@transform(filter_anndata, regex(r'(.*)_filt.h5ad'), r"logs/\1_subtractscrublet.log")
def subtract_scrublet(filt_obj, log_file):
    print(PARAMS['scrublet_treshold'])
    sprefix = PARAMS['sample_prefix']
    # we need subtract scrublet to run even if no threshold set, to write out cell_metadata
    cmd = """
         python %(py_path)s/subtract_scrublet.py
          --input_anndata %(filt_obj)s \
          --output_anndata %(filt_obj)s \
          --output_prefix %(sprefix)s \
          --scrublet_cutoff %(scrublet_threshold)s \
          --batch_col sample_id 
         """
    cmd += " > %(log_file)s 2>&1"
    P.run(cmd, job_threads=PARAMS['resources_threads_low'])


# @transform(subtract_scrublet, suffix("subtract_scrublet.log"), "postfilterplot.log")
@active_if(PARAMS['plotqc_metrics'] is not None)
@follows(subtract_scrublet)
@originate("logs/postfilterplot.log")
def postfilterplot(filter_log):
    sampleprefix = PARAMS['sample_prefix']
    metadatafile = sampleprefix + "_filtered_cell_metadata.tsv"
    cmd = """
    Rscript %(R_path)s/plotQC.R 
    --prefilter FALSE
    --metadata %(metadatafile)s 
    --sampleprefix %(sampleprefix)s
    --groupingvar %(plotqc_grouping_var)s
    --qcmetrics %(plotqc_metrics)s > %(filter_log)s 2>&1
    """
    P.run(cmd, job_threads=PARAMS['resources_threads_low']) 


@follows(subtract_scrublet)
@transform(filter_anndata, formatter("_filt.h5ad"), "logs/downsample.log")
def downsample(filt_obj, log_file):
    downsample_threshold = PARAMS['downsample']
    print(downsample_threshold)
    # downsample all
    if downsample_threshold is not None:
        cmd="""
        python %(py_path)s/downsample.py
         --input_anndata %(filt_obj)s \
         --output_anndata %(filt_obj)s \
         --downsample_value %(downsample)s \
         --batch_col sample_id 
        """
        cmd += " > %(log_file)s  2>&1"
        P.run(cmd, job_threads=PARAMS['resources_threads_low'])
    IOTools.touch_file(log_file)
    pass

# setting these follows means that it still works if filter is False
@follows(downsample, postfilterplot, filter_anndata, subtract_scrublet)
@transform(filter_anndata, regex(r'(.*)_filt.h5ad'), r"\1_log1p.h5ad")
def normalise_log_hvg_regress_scale(infile, outfile):
    cmd = """python %(py_path)s/normalise_log_hvg_regress_scale.py 
            --input_anndata %(infile)s 
            --output_prefix %(sample_prefix)s
            --fig_dir figures/"""
    # add in the options if specified in pipeline.yml
    if PARAMS['hvg_exclude'] is not None:
        if PARAMS['hvg_exclude'] == "default":
            hvg_exclude = PARAMS['resources_path'] + "/exclude_genes_HLAIGTR_v1.txt"
        cmd += " --exclude_file %(hvg_exclude)s"
    if PARAMS['hvg_flavor'] is not None:
        cmd += " --flavor %(hvg_flavor)s"
    if PARAMS['hvg_n_top_genes'] is not None:
        cmd += " --n_top_genes %(hvg_n_top_genes)s"
    if PARAMS['hvg_min_mean'] is not None:
        cmd += " --min_mean %(hvg_min_mean)s"
    if PARAMS['hvg_max_mean'] is not None:
        cmd += " --max_mean %(hvg_max_mean)s"
    if PARAMS['hvg_min_disp'] is not None:
        cmd += " --min_disp %(hvg_min_disp)s"
    if PARAMS['hvg_filter'] is not None:
        cmd += " --filter_by_hvg %(hvg_filter)s"
    if PARAMS['regress_variables'] is not None:
        cmd += " --regress_out %(regress_variables)s"
    if PARAMS['scale_max_value'] is not None:
        cmd += " --scale_max_value %(scale_max_values)s"
    cmd += " > logs/1_processing.log 2>&1"
    P.run(cmd, job_threads=PARAMS['resources_threads_medium'])


@follows(normalise_log_hvg_regress_scale)
@originate("logs/plot_pcas.log")
@mkdir("figures")
def plot_pcas(outfile):
    anndata_obj = PARAMS['sample_prefix'] + "_scaled.h5ad"
    cmd = "python %(py_path)s/plot_pca.py \
            --input_anndata %(anndata_obj)s \
            --fig_dir figures/"
    if PARAMS['pca_scree_n_pcs'] is not None:
        cmd += " --n_pcs %(pca_scree_n_pcs)s"
    if PARAMS['pca_color_by'] is not None:
        cmd += " --color_by %(pca_color_by)s"
    cmd += " > %(outfile)s 2>&1"
    P.run(cmd, job_threads=PARAMS['resources_threads_low'])


# ------------------------------------------------------------------------
# Integration
# ------------------------------------------------------------------------
# No batch correction
@follows(normalise_log_hvg_regress_scale)
@mkdir("batch_correction")
@originate("batch_correction/umap_bc_none.csv")
def run_no_batch(outfile):
    anndata_obj = PARAMS['sample_prefix'] + "_scaled.h5ad"
    cmd = """python %(py_path)s/no_batch_correction.py 
     --input_anndata %(anndata_obj)s
     --output_csv %(outfile)s
     --integration_col %(bc_column)s
     """
    if PARAMS['bc_npcs'] is not None:
        cmd += " --n_pcs %(bc_npcs)s"
    if PARAMS['bc_npcs'] is not None:
        cmd += " --n_neighbors %(bc_neighbors)s"
    cmd += " > logs/2_bc_no_correct.log 2>&1"
    P.run(cmd, job_threads=PARAMS['resources_threads_high'])


# BBKNN
@active_if(PARAMS['bc_run'])
@active_if('bbknn' in PARAMS['bc_tools'])
@follows(normalise_log_hvg_regress_scale)
@mkdir("batch_correction")
@mkdir("tmp")
@originate("batch_correction/umap_bc_bbknn.csv")
def run_bbknn(outfile):
    anndata_obj = PARAMS['sample_prefix'] + "_scaled.h5ad"
    cmd = """python %(py_path)s/batch_correct_bbknn.py 
     --input_anndata %(anndata_obj)s
     --output_csv %(outfile)s
     --integration_col %(bc_column)s
     """
    cmd += " > logs/2_bc_bbknn.log 2>&1"
    P.run(cmd, job_threads=PARAMS['resources_threads_high'])


# COMBAT
@active_if(PARAMS['bc_run'])
@active_if('combat' in PARAMS['bc_tools'])
@follows(normalise_log_hvg_regress_scale)
@mkdir("tmp")
@mkdir("batch_correction")
@originate("batch_correction/umap_bc_combat.csv")
def run_combat(outfile):
    anndata_obj = PARAMS['sample_prefix'] + "_scaled.h5ad"
    cmd = """python %(py_path)s/batch_correct_combat.py 
     --input_anndata %(anndata_obj)s
     --output_csv %(outfile)s
     --integration_col %(bc_column)s
     """
    if PARAMS['bc_npcs'] is not None:
        cmd += " --n_pcs %(bc_npcs)s"
    if PARAMS['bc_npcs'] is not None:
        cmd += " --n_neighbors %(bc_neighbors)s"
    cmd += " > logs/2_bc_combat.log 2>&1"
    P.run(cmd, job_threads=PARAMS['resources_threads_high'])


# HARMONY
@active_if(PARAMS['bc_run'])
@active_if('harmony' in PARAMS['bc_tools'])
@follows(normalise_log_hvg_regress_scale)
@mkdir("batch_correction")
@mkdir("tmp/")
@originate("batch_correction/umap_bc_harmony.csv")
def run_harmony(outfile):
    anndata_obj = PARAMS['sample_prefix'] + "_scaled.h5ad"
    cmd = """python %(py_path)s/batch_correct_harmony.py 
     --input_anndata %(anndata_obj)s
     --output_csv %(outfile)s 
     --integration_col %(bc_column)s
     """
    if PARAMS['bc_npcs'] is not None:
        cmd += " --n_pcs %(bc_npcs)s"
    if PARAMS['bc_neighbors'] is not None:
        cmd += " --n_neighbors %(bc_neighbors)s"
    if PARAMS['bc_sigma'] is not None:
        cmd += " --sigma_val %(bc_sigma)s"
    if PARAMS['bc_plotconvergence'] is not None:
        cmd += " --plotconvergence %(bc_plotconvergence)s"
    
    # cmd += " > logs/2_bc_harmony.log" the log is now created inside the script to allow capturin harmony logging
    P.run(cmd, job_threads=PARAMS['resources_threads_high'])

# SCANORAMA
@active_if(PARAMS['bc_run'])
@active_if('scanorama' in PARAMS['bc_tools'])
@follows(normalise_log_hvg_regress_scale)
@mkdir("batch_correction")
@mkdir("tmp/")
@originate("batch_correction/umap_bc_scanorama.csv")
def run_scanorama(outfile):
    anndata_obj = PARAMS['sample_prefix'] + "_scaled.h5ad"
    cmd = """python %(py_path)s/batch_correct_scanorama.py 
     --input_anndata %(anndata_obj)s
     --output_csv %(outfile)s 
     --integration_col %(bc_column)s
     """
    if PARAMS['bc_npcs'] is not None:
        cmd += " --n_pcs %(bc_npcs)s"
    if PARAMS['bc_neighbors'] is not None:
        cmd += " --n_neighbors %(bc_neighbors)s"
    cmd += " > logs/2_bc_scanorama.log 2>&1"
    P.run(cmd, job_threads=PARAMS['resources_threads_high'])


@collate([run_no_batch, run_scanorama, run_bbknn, run_harmony, run_combat], regex(r"(.*)/(.*)"),  r'\1/lisi.csv')
def run_lisi(infiles, outfile):
    print(infiles, outfile)
    infiles_string = ','.join(infiles)
    cmd = """python %(py_path)s/run_lisi.py 
    --input_files %(infiles_string)s 
    --batch_mtd batch_correction/batch_mtd.csv 
    --integration_col %(bc_column)s
    --output_file %(outfile)s 
    --fig_dir figures  > logs/3_lisi.log 2>&1
    """
    P.run(cmd, job_threads=PARAMS['resources_threads_low'])


@collate([run_no_batch, run_scanorama, run_bbknn, run_harmony, run_combat], regex(r"(.*)/(.*)"), r'batch_correction/combined_umaps.tsv')
def plot_umaps(infiles, outfile):
    print(infiles, outfile)
    sampleprefix = PARAMS['sample_prefix']
    metadatafile = sampleprefix + "_filtered_cell_metadata.tsv"
    infiles_string = ','.join(infiles)
    cmd = """python %(py_path)s/plot_umaps_batch_correct.py 
    --input_files %(infiles_string)s 
    --batch_mtd batch_correction/batch_mtd.csv
    --meta_df %(metadatafile)s
    --qc_metrics %(plotqc_metrics)s
    --groupingvar %(plotqc_grouping_var)s 
    --fig_dir figures/  > logs/4_plot_umaps.log
    """
    P.run(cmd, job_threads=PARAMS['resources_threads_low'])


@follows(run_lisi, plot_umaps)
@originate("logs/batch_correction_complete.log")
def batch_correction(outfile):
    IOTools.touch_file(outfile)


# Final choices: Leave blank until you have reviewed the results from running
# sc_pipelines integration make full
# choose the parameters to integrate into the final anndata object
# then run
# `sc_pipelines integration make merge_batch_correction`
@originate("logs/merge_final_obj.log")
def merge_batch_correction(outfile):
    scaled_obj = PARAMS['sample_prefix'] + "_scaled.h5ad"
    log1p_obj = PARAMS['sample_prefix'] + "_log1p.h5ad"
    cmd = """
    python %(py_path)s/batch_correct_merge.py 
    --scaled_anndata %(scaled_obj)s
    --correction_choice %(bc_choice)s
    """
    cmd += " > %(outfile)s"
    P.run(cmd, job_threads=PARAMS['resources_threads_medium'])
    # clear up tmp since it is no longer required
    P.run("rm -r tmp")


@follows(plot_pcas,postfilterplot, batch_correction)
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
