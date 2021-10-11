# Plot cluster metrics based on clusters from python/run_clustering.py
# Script originally written by Fabiola Curion, adapted to take python outputs by Charlotte Rich-Griffin
####


library(optparse)
library(Seurat)
library(dplyr)
library(tidyselect)
library(ggplot2)
library(foreign)
library(xtable)
library(stringr)
library(reticulate)
library(ggrepel)

set.seed(42)
option_list <- list(
	make_option(c("--mtd_object"), default="none",help="seurat or scanpy object to take as input"),
	make_option(c("--mtd_type"), default=NULL, help="specify whether a seurat or scanpy object \
	              is  the input type, leave blank for txt file, default=NULL"),
	make_option(c("--coords"), default="none",
							help="A txt file containing coordinates froma dimension reduction"),
	make_option(c("--coords_suffix"), default="", help="string to add to umap plots filenames, optional"),
	make_option(c("--clusterid"), default="none",
							help="An rds object containing the cluster identities"),
	make_option(c("--continuous_variables"), default='percent.mito',
				help="continuous variables to plot, csv list of variables stored in adata.obs"),
	make_option(c("--discrete_variables"), default="condition",
							help="variables to use to plot. Should be a comma separated list"),
	make_option(c("--dimred"), default='umap',
							help="What dim reduction to use for the cell dimred space. Default to umap. Needs a tsne.txt table to be precomputed (pca,umap) "),
	make_option(c("--outfile_prefix"), default=NULL,
							help="prefix for output files, (should contain parameters)")
	)

opt <- parse_args(OptionParser(option_list=option_list))
message("Running with options:")
print(opt)

message("Saving to:")
if(grepl('/', opt$outfile_prefix)){
	out_path = paste0(dirname(opt$outfile_prefix), '/')
}else{
	out_path='./'
}
print(out_path)



# read in metadata
if(is.null(opt$mtd_type)){
	mtd <- read.table(opt$mtd_object, sep='\t', header=TRUE, row.names=1)
}else if(opt$mtd_type == "scanpy"){
	sc <- reticulate::import("scanpy")
	# use scanpy to read in anndata
	adata= sc$read_h5ad(opt$mtd_object)
	mtd = adata$obs
} else if(opt$mtd_type == "seurat"){
	seu = readRDS(opt$mtd_object)
	mtd = seu@meta.data
} else{
	mtd <- read.table(opt$mtd_object, sep='\t', header=TRUE, row.names=1)
}


# read in coordinates
tt<-read.table(opt$coords, header=T, sep='\t', row.names=1)
colnames(tt) <- c('UMAP_1', 'UMAP_2')


cont_plot <- strsplit(opt$continuous_variables,",")[[1]]
print("continuous variables to plotted")
print(cont_plot[which(unlist(cont_plot) %in% colnames(mtd))])
print("continuous variables not found in metadata")
print(cont_plot[!(which(unlist(cont_plot) %in% colnames(mtd)))])


# merge mtd and tt
tt <- cbind(tt, mtd)

# read in cluster id
cl.id <- read.table(opt$clusterid, sep='\t', header=TRUE, row.names = 1)
print(head(cl.id))

tt$cluster <- factor(cl.id$cluster)

##dim reduction plots 

if(opt$dimred=="umap"){
	plout = opt$outfile_prefix
	print("plotting prefix")
	print(plout)

  # plot continuous variables on umap 
  print("---------------")
  print("plot continuous variables" )
  print("---------------")

	for(cp in cont_plot){
		print(cp)
		tt %>% ggplot(aes_string("UMAP_1","UMAP_2", color=cp)) +
			geom_point(size=0.3) + scale_color_viridis_c() + theme_bw() +
			theme(legend.key.height=unit(0.5,'cm'), legend.key.width=unit(0.1,'cm'), axis.title = element_text(size=15), axis.text = element_text(size=14)) -> g
		ggsave(g, filename=paste0(opt$outfile_prefix,"_",cp,"_", opt$coords_suffix, ".png" ), width=7, height=6, type="cairo-png", dpi=300)
	}
	
	# plot discrete variables on umap 
	
	varplot1 <- strsplit(opt$discrete_variables,",")[[1]]
	if(length(varplot1) ==1){
		use.var1 = varplot1
	}else{
		use.var1 <- varplot1[varplot1 %in% colnames(tt)]
	}
	print("Plotting with")
	print(use.var1)
	# check vars are not integers in tt
	tt[,use.var1] <- tt %>% select_at(.vars=use.var1) %>% mutate_all(as.factor)
	

	print("---------------")
	print("plot discrete variables" )
    print("---------------") 
	for (vv in use.var1){
		print(vv)
		tt %>% ggplot(aes_string("UMAP_1", "UMAP_2", color=vv)) + geom_point(size=0.3) + labs(color=vv) +
			theme_bw() +
			theme(legend.key.height=unit(0.5,'cm'), legend.key.width=unit(0.2,'cm'), axis.title = element_text(size=15), axis.text = element_text(size=14)) -> g

		if (vv=="cluster"){
          tt %>% select(UMAP_1,UMAP_2, cluster) %>% group_by(cluster) %>% summarise_all( .funs = median) ->dm
          g <- g + 
               geom_text_repel(data = dm, aes(UMAP_1,UMAP_2, label=cluster), color="black", size=3) +
               guides(colour = guide_legend(override.aes = list(shape=16,size=2)))  
        }

		ggsave(g, filename=paste0(opt$outfile_prefix,"_umap_", vv, "_", opt$coords_suffix, ".png" ), width=7, height=6, type="cairo-png", dpi=300)
	}


}



##barplots
print("---------------")
print("plot barplots")
print("---------------")
use.var1 <- use.var1[!use.var1=="cluster"]
for (vv in use.var1){
	gtools::mixedsort(as.character(unique(tt[,vv]))) ->setlim
	print(vv)
	#print(tt %>% group_by_at(.vars=c("cluster",vv)) %>% summarize(tot=n()))
	tt %>% group_by_at(.vars=c("cluster",vv)) %>% summarize(tot=n()) %>%
		ggplot(aes_string("cluster", "tot", fill=vv ))+ geom_bar(stat="identity", position="dodge", color="black") +
		labs(fill=vv) +
		theme_bw() +
		theme(legend.key.height=unit(0.5,'cm'), legend.key.width=unit(0.1,'cm'),
					axis.text = element_text(size=12), axis.title = element_text(size=13, face="bold")) +
		ylab("Number of cells") -> g
	if(length(unique(tt[,"cluster"])) >12) { ww=23}else{ww=12}

	ggsave(g, filename=paste0(opt$outfile_prefix,"_barplot_cluster_",vv,".pdf" ), width=ww, height=6, dev=cairo_pdf)
}

for (vv in use.var1){
	gtools::mixedsort(as.character(unique(tt[,vv]))) ->setlim
	tt %>% group_by_at(.vars=c(vv,"cluster")) %>% summarize(tot.per.var=n()) %>%
		group_by_at("cluster") %>%
		mutate(tot.per.cluster= sum(tot.per.var)) %>%
		mutate(percent=100*tot.per.var/tot.per.cluster)%>%
		ggplot(aes_string("cluster","percent", fill=vv)) + geom_bar(stat="identity", position="stack", color="black") +
		labs(x="cluster") + theme_bw() +
		theme(axis.text.x = element_text(angle=30, hjust=1, size=12),
					axis.text.y = element_text(size=12),
					axis.title = element_text(size=13, face="bold"),legend.key.height=unit(0.5,'cm'), legend.key.width=unit(0.1,'cm')) +
		ylab("Percent of cells") -> g
	if(length(unique(tt[,"cluster"])) >12) { ww=23}else{ww=12}
	ggsave(g, filename=paste0(opt$outfile_prefix,"_proportion_cluster_",vv,".pdf" ), width=ww, height=6, dev=cairo_pdf)
}

##violin plots
print("---------------")
print("plot violins")
print("---------------")

for (vv in use.var1){
	gtools::mixedsort(as.character(unique(tt[,vv]))) ->setlim
	g <- tt %>%  reshape2::melt(id.vars=c("cluster",vv), measure.vars=cont_plot) %>%
		ggplot( aes_string("cluster","value", fill=vv)) +geom_violin()  + facet_wrap(~variable, scales='free', ncol=min(length(cont_plot),2)) +
		theme_bw() +
		theme(legend.key.height=unit(0.5,'cm'), legend.key.width=unit(0.1,'cm'),
			  strip.background.x = element_rect(fill="white", color="white"),
			  strip.text.x = element_text(size=10, face='bold.italic'),
              axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
		ylab("QC value") 

	if(length(unique(tt[,"cluster"])) >12) { ww=23}else{ww=12}
	ggsave(g, filename=paste0(opt$outfile_prefix,"_violin_cluster_",vv,".pdf" ), width=ww, height=10, dev=cairo_pdf)
}

print("done")