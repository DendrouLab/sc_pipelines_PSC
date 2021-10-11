# if(interactive()){reticulate::use_python("/Users/crg/Documents/Projects/pipelines/bin/pipeline_venv/bin/python")} #why this?
library(optparse)
library(Seurat)
library(dplyr)
library(magrittr)
library(ComplexHeatmap)
library(RColorBrewer)
library(ggplot2)
library(reticulate)
library(Matrix)
set.seed(42)

get_parameter_from_fname <- function(fname, val){
  # extract reoltuion from filename
  # INPUT example: "res0.2_clustering.txt"
  # OUTPUT example: "0.2"
  # ------------------------
  # split to get parameters
  f_split = strsplit(fname, '_')[[1]]
  # get resolution string
  res_str = f_split[grep(val, f_split)]
  # extract number
  return(gsub(val, "", res_str))
}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

option_list <- list(
  make_option(c("--anndata"), default="./data/run_scanpy_integration/pbmc3k_log1p.h5ad",
              help=""),
  make_option(c("--scanpy_layer"), default="X",
              help=""),
  make_option(c("--marker_file"),
              help=""),
  make_option(c("--output_file"), default="./heatmap.pdf",
              help="name of txt file for output"),
  make_option(c("--cluster_file"),
              help="An rds object containing the cluster identities"),
  make_option(c("--metadata_file"),
              help="An rds object containing the cluster identities"),
  make_option(c("--docomplex"), default=FALSE,
              help="Whether to plot complex heatmap. Default (if FALSE is set) to Seurat Doheatmap.If true, specify subgroups"),
  make_option(c("--subgroup"), default="batch",
              help="comma separated list of factors to group by if doing complex heatmap.")
)

opt <- parse_args(OptionParser(option_list=option_list))
#
# if(interactive()){
#   setwd("/Users/crg/Documents/Projects/pipelines/data/run_clustering")
#   opt$anndata="../run_integration/alt_start_scaled.h5ad"
#   opt$marker_file="alt_start_nneigh30_pcs50_dir/alt_start_algleiden_res0.2/markers/markers_list_top.txt.gz"
#   opt$cluster_file="alt_start_nneigh30_pcs50_dir/alt_start_algleiden_res0.2/algleiden_clusters.txt.gz"
#   opt$output_file="alt_start_nneigh30_pcs50_dir/alt_start_algleiden_res0.2/figures/heatmap.png"
#   opt$metadata_file="alt_start_cell_mtd.txt"
#   opt$docomplex=TRUE
#   opt$subgroup="sample_id,sex"
# }
message("Running with options:")
print(opt)
message("Saving to:")

out_path = dirname(opt$output_file)
print(out_path)


message("Reading cluster info...")
cl.id<- read.table(opt$cluster_file, sep='\t', header=T, row.names=1)


# use reticulate to read the scaled data and metadata

# read scaled data
message("Using reticulate to import scanpy")
sc <- reticulate::import("scanpy")
# use scanpy to read in anndata
adata= sc$read_h5ad(opt$anndata)
if(opt$scanpy_layer=="X"){
  # exprs = as(adata$X,"matrix")
  exprs <- as(t(adata$X), "CsparseMatrix")
}else{
  exprs = t(adata$layers[opt$scanpy_layer])
}
colnames(exprs) <- adata$obs_names$to_list()
rownames(exprs) <- adata$var_names$to_list()



message("Reading metadata...")
mtd = data.table::fread(opt$metadata_file) %>% as.data.frame
row.names(mtd) = mtd$cellbarcode


mtd$clusters = cl.id$clusters


# opT$docomplex is a string, convert to Bool
docomplex <- as.logical(opt$docomplex)


message("Reading markers...")
gde.all <- read.table(opt$marker_file, header=T)

message("Reading cluster info...")
table(cl.id$clusters) -> sk
data.frame(sk) -> tmp
colnames(tmp) <- c("cluster","ncells")

summary_stats <- gde.all %>% group_by(cluster) %>% dplyr::summarize(n_pos_markers=sum(avg_logFC>0), n_neg_markers=sum(avg_logFC<0), n_markers=n())
summary_stats <- merge(summary_stats,tmp, by="cluster")
summary_stats<- summary_stats[,c("cluster", "ncells","n_pos_markers" ,"n_neg_markers" ,"n_markers")]


summary_stats_fname=gsub(".txt", ".summary.stats.txt", opt$cluster_file)
write.table(summary_stats,
            file=summary_stats_fname,
            quote=F,sep="\t",row.names=F)

top_n <- gde.all %>% filter(avg_logFC>0) %>% group_by(cluster) %>% top_n(10,scores)

nclust <-length(sk)

if(nclust>=12) {
  ww<-20
  hh<-18
} else {
  ww<-10
  hh<-12
}
min_max_clip<- function(x, a, b) {
  ifelse(x <= a,  a, ifelse(x >= b, b, x))
}

# which genes and cells
genes_use <-as.character(top_n$gene)
genes_use <- genes_use[genes_use %in% rownames(exprs)]


cell_use <- row.names(cl.id)


# define colours for heatmap (match seurat colors)
exprs_cols <- circlize::colorRamp2(c(-2.5,0,2.5), c("#FF00FF", "#000000", "#FFFF00"))
cell_clusters <- cl.id$clusters #reset the order FIx

if(docomplex){
  
  sub_group <- strsplit(opt$subgroup,",")[[1]]
  
  sub_group[sub_group %in% colnames(mtd)] -> plotvar
  if(!is.null(plotvar)) {
    message("Plotting with")
    print(plotvar)
  }else{
    stop("specified sub group not found in the metadata")
  }
  
  for(sg in plotvar){
    message("Plotting complex heatmap")
    message(sg)
    
    # get the expression data
    x <- as.matrix(exprs[genes_use,cell_use])
    x <- min_max_clip(x, -2.5, 2.5) # seurat default values
    n_rows <- nrow(x)
    row.cex <- min(0.5, 60/n_rows)
    
    #clusters <- exprs[,cell_use]
    # set up the subgroup colour palette
    cell_sub_groups <- mtd[cell_use, sg]
    # account for missing data
    cell_sub_groups[is.na(cell_sub_groups)] = "missing"
    cell_sub_groups[cell_sub_groups==""] = "missing"
    sub_groups <- gtools::mixedsort(unique(cell_sub_groups))
    
    # sub_group_cols <- brewer.pal(length(sub_groups),"Set1")
    sub_group_cols <- gg_color_hue(length(sub_groups))
    sub_group_cols <- sub_group_cols[1:length(sub_groups)]
    names(sub_group_cols) <- sub_groups
    
    print(sub_group_cols)
    
    cell_order <- order(cell_clusters, cell_sub_groups)
    
    subgroupAnnotation = HeatmapAnnotation(df=data.frame(subgroup=cell_sub_groups[cell_order]),
                                           col = list(subgroup=sub_group_cols),
                                           show_annotation_name = FALSE)
    
    x <- x[,cell_order]
    
    # draw the heatmap
    ComplexHeatmap::Heatmap(x,
                            cluster_rows = FALSE,
                            col = exprs_cols,
                            row_names_gp = gpar(#fontsize = row_names_gp,
                              cex = row.cex),
                            cluster_columns = FALSE,
                            show_column_names = FALSE,
                            row_title = NULL,
                            row_split=top_n %>% filter(gene  %in% genes_use) %>% .$cluster,
                            column_split = cell_clusters[cell_order],
                            top_annotation = subgroupAnnotation,
                            row_gap=unit(0.7,"mm"),
                            column_title_side="top",
                            column_gap=unit(0.75, "mm"),
                            use_raster=TRUE,
                            raster_device="CairoPNG",
                            raster_quality=4,
                            show_heatmap_legend = TRUE,
                            name="Scaled\nExprs") ->gg
    
    png(paste0(gsub(".png","_",opt$output_file), sg,".png"), width=ww, height=hh,units = 'in', res = 120, type="cairo")
    print(gg)
    dev.off()
    
  }
}
# also make a heamtap with no variables

# get the expression data
x <- as.matrix(exprs[genes_use,cell_use])
x <- min_max_clip(x, -2.5, 2.5) # seurat default values
n_rows <- nrow(x)
row.cex <- min(0.5, 60/n_rows)

sub_group = 'clusters'
print(sub_group)
# set up the subgroup colour palette
cell_sub_groups <- mtd[cell_use, sub_group]
sub_groups <- gtools::mixedsort(unique(cell_sub_groups))

# sub_group_cols <- brewer.pal(length(sub_groups),"Set1")
sub_group_cols <- gg_color_hue(length(sub_groups))
sub_group_cols <- sub_group_cols[1:length(sub_groups)]
names(sub_group_cols) <- sub_groups

print(sub_group_cols)

cell_order <- order(cell_clusters, cell_sub_groups)

subgroupAnnotation = HeatmapAnnotation(df=data.frame(cluster=cell_sub_groups[cell_order]),
                                       col = list(cluster=sub_group_cols),
                                       show_annotation_name = FALSE)

x <- x[,cell_order]
# match the Seurat colors
ComplexHeatmap::Heatmap(x,
                        cluster_rows = FALSE,
                        col = exprs_cols,
                        row_names_gp = gpar(#fontsize = row_names_gp,
                          cex = row.cex),
                        cluster_columns = FALSE,
                        show_column_names = FALSE,
                        row_title = NULL,
                        row_split=top_n %>% filter(gene  %in% genes_use) %>% .$cluster,
                        column_split = cell_clusters[cell_order],
                        top_annotation = subgroupAnnotation,
                        row_gap=unit(0.7,"mm"),
                        column_title_side="top",
                        column_gap=unit(0.75, "mm"),
                        use_raster=TRUE,
                        raster_device="CairoPNG",
                        raster_quality=4,
                        show_heatmap_legend = TRUE,
                        name="Scaled\nExprs") ->gg



png(opt$output_file, width=ww, height=hh,units = 'in', res = 120, type="cairo")
print(gg)
dev.off()


message("done")