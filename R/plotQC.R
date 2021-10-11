# plotting utilities for qc 
suppressPackageStartupMessages({library(tidyverse)
library(magrittr)
library(circlize)
library(ComplexHeatmap)
library(R.utils)
library(optparse)
library(Matrix)
library(matrixStats)})
options(stringsAsFactors = F)


option_list <- list(
  make_option(c("--prefilter"), default=TRUE,
              help="am i parsing data before or after filtering?"),
  make_option(c("--metadata"), default=NULL,
              help="the path to the anndata object"),
  make_option(c("--sampleprefix"), default=NULL,
              help="the prefix to prepend to save summary filtering plots"),
  make_option(c("--groupingvar"), default="sample_id,tissue,patient,channel",
              help="names of grouping variables"),
  make_option(c("--qcmetrics"), default="pct_counts_mt,pct_counts_rp",
              help="names of qcmetrics"),
  make_option(c("--outdir"), default="./figures/",
              help="the name of the output folder")
)


message("Plot QC data")
opt <- parse_args(OptionParser(option_list=option_list))
if(is.null(opt$outdir)) { opt$outdir <- paste0(getwd(),"/")}
if(!grepl("\\/$", opt$outdir)){opt$outdir <- paste(opt$outdir, "/", sep = "")}
if (!file.exists(opt$outdir)){ dir.create(opt$outdir)}
print("Running with options:")
print(opt)


run<- opt$outdir
opt[which(opt=="NULL")] <- NULL
opt[which(opt=="None")] <- NULL
opt$prefilter <- as.logical(opt$prefilter)

data_plot <- read.delim(opt$metadata, header=T, sep="\t", row.names = 1)
sprefix <- opt$sampleprefix

if (!is.null(opt$groupingvar)){
  source_facet <- strsplit(opt$groupingvar,",")[[1]]
  source_facet[source_facet %in% colnames(data_plot)] -> source_facet
  if(length(source_facet)>0){
  message("Plotting with")
  print(source_facet)
  }
    # add sample_id as a minimum requirement if it's not there already
  source_facet = c("sample_id", source_facet)

}else{
  stop("i don't have the minimum variable _sampleid_ to facet on, will stop here")
}


if (!is.null(opt$qcmetrics)){
  qcmetrics <- strsplit(opt$qcmetrics,",")[[1]]
  qcmetrics[qcmetrics %in% colnames(data_plot)] ->qcmetrics
  if(length(qcmetrics)>0){
  message("Qc metrics i will use are")
  print(qcmetrics)
  }




for (qc in qcmetrics){
  print(qc)
  for (sc in source_facet[!source_facet=="sample_id"]){
   for (ss in source_facet[!source_facet=="sample_id"]){

     
     
       g <- data_plot %>%
        ggplot(aes_string(x="sample_id",y=qc)) +
        {if(ss!=sc) geom_violin(aes(fill=interaction(get(ss),get(sc))))}+
        {if(ss==sc) geom_violin(aes_string(fill=ss))}+ 
        {if(qc=="doublet_scores") geom_hline(yintercept=0.25) }+
        {if(any(qc %in% c('pct_counts_mt','pct_counts_rp', 'pct_counts_hb','pct_counts_ig'))) geom_hline(yintercept=c(5, 10, 20 ,70), color="seagreen1")}+
        {if(any(qc %in% c('pct_counts_mt','pct_counts_rp', 'pct_counts_hb','pct_counts_ig'))) coord_cartesian(ylim=c(0,100))}+
        theme_bw()+
        theme(axis.text.x=element_text(size=8,angle=45, hjust=1, vjust=1),
              axis.text.y=element_text(size=13, face="bold"))
        
        
        strname <- ifelse(sc==ss, sc,paste(sc,ss,sep=".")) 
        ggsave(g, filename=paste0(run, "sample_", strname, "_", qc,"_violin.pdf"), width= 8, height=7)
    }
  }
  
}


data_plot %>% 
  ggplot(aes( n_genes_by_counts,doublet_scores)) +
  {if("pct_counts_hb" %in% colnames(data_plot)) geom_point(aes(size=pct_counts_hb, color=pct_counts_mt))}+
  {if("pct_counts_hb" %in% colnames(data_plot)) scale_size(range=c(0,1))}+
  {if(!("pct_counts_hb" %in% colnames(data_plot))) geom_point(aes(color=pct_counts_mt)) }+
  scale_color_viridis_c() +
  geom_hline(yintercept=0.25, color="red", linetype="dashed") +
  geom_vline(xintercept=c(100,2500,3000), color="red", linetype="dashed") + 
  theme_bw() + theme(axis.text = element_text(size=8, face="bold", color="black"))+
  facet_wrap(~sample_id) ->g 
ggsave(g, filename=paste0(run, "sample_scatter_multifilter_ptcmito.png"), width= 12, height=10, type="cairo-png", dpi=300)

data_plot %>% 
  ggplot(aes( n_genes_by_counts,doublet_scores)) +
  {if("pct_counts_hb" %in% colnames(data_plot)) geom_point(aes(size=pct_counts_hb, color=total_counts))}+
  {if("pct_counts_hb" %in% colnames(data_plot)) scale_size(range=c(0,1))}+
  {if(!("pct_counts_hb" %in% colnames(data_plot))) geom_point(aes(color=total_counts)) }+
  scale_color_viridis_c() +
  geom_hline(yintercept=0.25, color="red", linetype="dashed") +
  geom_vline(xintercept=c(100,2500,3000), color="red", linetype="dashed") + 
  theme_bw() + theme(axis.text = element_text(size=8, face="bold", color="black"))+
  facet_wrap(~sample_id) ->g 
ggsave(g, filename=paste0(run, "sample_scatter_multifilter_numi.png"), width= 12, height=10, type="cairo-png", dpi=300)

if(opt$prefilter){

  message ("saving some counts tables for references")

    data_plot %>% 
    filter(pct_counts_mt<=10 & pct_counts_hb<=70 &n_genes_by_counts>=100 & doublet_scores<=0.25) %>%
    group_by_at(.vars=c(source_facet)) %>%
    group_by_at(.vars=c(source_facet)) %>%
    summarise(cell.count= n()) %>%
    group_by_at(.vars="sample_id") %>% 
    rename(cell.count_f1=cell.count) ->f1
    
    data_plot %>% 
    filter(pct_counts_mt<=5 & pct_counts_hb<=50 &n_genes_by_counts>=100 & doublet_scores<=0.25) %>%
    group_by_at(.vars=c(source_facet)) %>%
    summarise(cell.count= n()) %>% 
    group_by_at(.vars="sample_id") %>% 
    rename(cell.count_f2=cell.count) ->f2
    
  data_plot %>% 
    filter(pct_counts_mt<=5 & pct_counts_hb<=50 &n_genes_by_counts>=100 & n_genes_by_counts<=3000) %>%
    group_by_at(.vars=c(source_facet)) %>%
    summarise(cell.count= n()) %>% 
    group_by_at(.vars="sample_id") %>% 
    rename(cell.count_f3=cell.count) ->f3

  data_plot %>% 
    group_by_at(.vars=c(source_facet)) %>%
    summarise(cell.count= n()) %>% 
    group_by_at(.vars="sample_id") %>% 
    rename(baseline.counts=cell.count) ->baseline 

  merge(merge(merge(f1,f2,by=source_facet,all=TRUE),f3, by=source_facet, all=TRUE), baseline, by=source_facet, all=T) %>%
  mutate(percent_retain_f1 = 100*cell.count_f1/baseline.counts,
        percent_retain_f2 = 100*cell.count_f2/baseline.counts,
        percent_retain_f3 = 100*cell.count_f3/baseline.counts) ->info

  

  write.table(info, file=paste0( sprefix,"_threshold_filter.tsv"), col.names=T, row.names=F, sep="\t", quote=F)

  data.frame(qcmetric=c("pct_counts_mt_max","pct_counts_hb_max","n_genes_by_counts_min","doublet_scores_max","n_genes_by_counts_max"),
  f1=c(10,70,100,0.25,NA),
  f2=c(5,50,100,0.25,NA),
  f3=c(5,50,100,NA,3000)) ->explain

  write.table(explain, file=paste0(sprefix,"_threshold_filter_explained.tsv"), col.names=T, row.names=F, sep="\t", quote=F)

  lab <- c("%Mt <=10 & %HB<70 & minGenes>=100 & scrublet<=0.25",
          "%Mt <=5 & %HB<50 & minGenes>=100 & scrublet<=0.25",
          "%Mt <=5 & %HB<50 & minGenes>=100 & maxGenes<=3000")
  names(lab) <- c("percent_retain_f1","percent_retain_f2","percent_retain_f3") 

  info %>% 
    pivot_longer(cols=starts_with("percent"), names_to="filterclass", values_to="percent")%>%
    ggplot(aes(sample_id,percent, fill=filterclass)) +
    geom_bar(stat="identity", position="dodge", color="black") +
    facet_wrap(~filterclass, ncol=1, labeller = labeller(filterclass=lab)) +
    theme_bw()+
    theme(axis.text.x=element_text(size=8,angle=45, hjust=1, vjust=1),
          axis.text.y=element_text(size=13, face="bold")) +
    scale_fill_manual(values = c("red", "yellow","blue"), 
    breaks=c("percent_retain_f1","percent_retain_f2","percent_retain_f3"),
    limits=c("percent_retain_f1","percent_retain_f2","percent_retain_f3")) + 
    coord_cartesian(ylim=c(0,100)) -> g

  ggsave(g, file = paste0(sprefix,"_thresholds_filter.pdf"), width=9, height=9)
} else{
  message("producing filed with final counts for cells after filtering")

  data_plot %>% 
  group_by_at(.vars=c(source_facet)) %>% 
  summarise(cell.count= n()) %>% 
  group_by_at(.vars="sample_id") ->baseline 
  
  write.table(baseline, file=paste0( sprefix,"_filtered_data.tsv"), col.names=T, row.names=F, sep="\t", quote=F)
  
  baseline %>% 
  ggplot(aes_string(x="sample_id", y="cell.count")) +
  geom_bar(stat="identity", position="dodge", color="black", fill="red") +
  theme_bw()+
  theme(axis.text.x=element_text(size=8,angle=45, hjust=1, vjust=1),
        axis.text.y=element_text(size=13, face="bold")) -> g
  ggsave(g, file = paste0(run,sprefix,"_filtered_data.pdf"), width=9, height=6)
  

}
}else{
  message("i don't have the QC metrics to use!")
}

message("done")
