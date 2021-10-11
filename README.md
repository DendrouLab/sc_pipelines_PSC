# Dendrou group single cell pipelines

- Maintained by Charlotte Rich-Griffin
- Contributors: Charlotte Rich-Griffin, Tom Thomas and Fabiola Curion

### Available pipelines:
- qc_cellranger
- qc_general_start
- integration
- clustering

Coming soon:
- demultiplexing

## Introduction
These pipelines use cgat-core pipeline software


## Installation:
It is advisable to run everything in a virtual environment.

e.g.
```
mkdir my-project
cd my-project
python3 -m venv --prompt=sc_pipelines python3-venv-scpipelines/
# This will create a my-project/venv folder
```

activate the environment

```
cd my-project
source python3-venv/bin/activate
```

Download and install this repo

```
git clone https://github.com/DendrouLab/sc_pipelines
cd sc_pipelines
pip install --upgrade pip
pip install -r requirements.txt # installs required python packages
python setup.py develop
```
The pipelines are now installed as a local python package.

The pipelines use R (mostly for ggplot visualisations). The pipeline will call a local R installation (as opposed to requireing a specific build within the virtual environment)
Install required R packages by copying the following code into R
```
install.packages(c())
```

Create a yml file for the cgat core pipeline software to read

```
vim ~/.cgat.yml
```
containing the following information
```
cluster:
    queue_manager: sge
```



To check the installation was successful run the following line
```
sc_pipelines --help
```
A list of available pipelines should appear!

## General principles for running pipelines:
Run the pipeline for the login node on your server, it will use in built the job submission system to submit jobs.

Navigate to the directory where you want to run your analysis (this should not be within the dendrou_pipelines folder)
```
mkdir data_dir/
cd data_dir/
sc_pipelines qc_cellranger config
```

This will produce two files, `pipeline.log` and `pipeline.yml`

Edit `pipeline.yml` as appropriate for your data.

Then check which jobs will run with the command
```
sc_pipelines qc_cellranger show full
```
The output of this will show a list of tasks that will be run as part of the pipeline.

To run use the command
```
sc_pipelines qc_cellranger make full
```
or
to run it in the background, and prevent the jobs from hanging up when you log off the server
```
nohup sc_pipelines qc_cellranger make full &
```

Occasionally you might want to run tasks individually (e.g. to debug)
In order to do this you can run any task in the `show full` list such as:
```

sc_pipelines clustering make find_markers
```

## Running the complete pipeline

Run each of pipeline qc, integration and clustering in separate folders.

1. Run `sc_pipelines qc_cellranger make full `
2. Use outputs to decide filtering thresholds. Note that the actual filtering occurs in the first step of integration pipeline
3. Run `sc_pipelines integration make full`
4. Use outputs to decide on the best batch correction method
5. Edit the pipeline yml with your preferred batch correction 
6. Run `sc_pipelines integration make merge_batch_correction`
7. Run the clustering pipeline  `sc_pipelines clustering make full`

## Inputs to QC pipeline
There are two qc pipelines
- qc_cellranger


For qc_cellranger the minimum required columns are 

sample_id  |     raw_path   |     filtered_path   
---|---|---
Example at resources/sample_file_cellranger.txt


#### demultiplexing data
If you have demultiplexing data you can should include two extra columns in your samples 
file (e.g. resources/sample_file_inc_demultiplexing.txt); demultiplex_map_file  and   
demultiplex_mtd_file examples at resource/demult_map.csv and resources/demult_mtd.csv, 
espectively. demultiplexing_map_file should contain all barcodes, "antibody" 
demultimplexing_mtd_file should contain "antibody" which corresponds to the map file and any other metadata associated to your samples
Other metadata that you want to add to the amnndata object can be specified in 
the pipeline.yml.

Note that if you are combining multiple datasets from different sources the final anndata object will only contain the intersection of the genes
from all the data sets. For example if the mitochondrial genes have been excluded from one of the inputs, they will be excluded from the final data set.
In this case it might be wise to run qc separately on each dataset, and them merge them together to create on h5ad file to use as input for
integration pipeline.

## Running pipeline modules separately
For circumstances where you arelady have a qc'd anndata object
### Integration

It is possible to run the integration pipeline starting from one combined anndata object containing all your samples containing raw counts in the X slot, 
either with or without running filtering first.
If your data is set to run simply call your anndata object [PROJECT_PREFIX]_filt.h5ad and set filtering_run: False.
You must have a column called sample_id which groups the data in some way (otherwise the plotting script will break) TODO: Fix this

If you have not filtered your data then you can set run filtering_run: True, and set the remaining parameters, 
BUT you must ensure that your obs columns names which you want to use for filtering match the column names in [resources/qc_column_names.txt](https://github.com/DendrouLab/sc_pipelines/blob/matrix_start/resources/qc_column_names.txt)

### Clustering
To run clustering_scanpy without the prior steps, you will need to produce 2 anndata objects
[PROJECT_PREFIX]_log1p.h5ad and [PROJECT_PREFIX]_scaled.h5ad 

[PROJECT_PREFIX]_log1p.h5ad  
- log normalised data in the adata.X slot
- highly variabel genes calculated

[PROJECT_PREFIX]_scaled.h5ad
- log normalised data saved in adata.raw.X
- scaled data (optionally regress) in adata.X
- pca

Minimal code:
```
sc.pp.normalize_total(adata, target_sum=1e4);
sc.pp.log1p(adata))
sc.pp.highly_variable_genes(adata)
adata.write("anndata_log1p.h5ad")

# optional steps
# adata = adata[:, adata.var.highly_variable]
# sc.pp.regress_out()

sc.pp.scale(adata)
sc.tl.pca(adata)
adata.write("anndata_scaled.h5ad")
``` 

If you want to use a specific batch correction then fit it into the above minimal code as appropriate

(It's probably easier to just run integration again with `filtering_run: False`)

## Plans:
- adapt demultiplexing addition to qc to work for any per barcode metadata?

CRG 2021-02-25

