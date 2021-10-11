# Clustree diagram from aggregated cluster data.
# Script originally written by Fabiola Curion, adapted by Charlotte Rich-Griffin
####


library(clustree)
library(dplyr)
library(optparse)
library(ggplot2)
library(readr)

get_parameter_from_fname <- function(fname, val){
	# extract reoltuion from filename
	# INPUT example: "res0.2_clustering.txt"
	# OUTPUT example: "0.2"
	# ------------------------
	# split to get parameters
	f_split = strsplit(fname, '/|_')[[1]]
	# get resolution string
	res_str = f_split[grep(val, f_split)]
	# extract number
	return(gsub(val, "", res_str))
}

option_list <- list(
    make_option(c("--infiles"), default="none",
                help="a comma separated list of clustering files"),
	make_option(c("--infile_prefix"), default="none",
                help="file prefix for clustree to use as input"),
  	make_option(c("--outfile"), default="clustree.pdf",help = "save file")
    )

opt <- parse_args(OptionParser(option_list=option_list))

# 1. mload and merge infiles
# 2. name columns
# 3. run clustree

message("Running with options:")

print(opt)


# read in files
#lf <- list.files(path="../../data/scanpy_pipe_test", pattern = "cluster.txt", full.names = T)
lf <- strsplit(opt$infiles,",")[[1]]
alg <- get_parameter_from_fname(lf[1], 'alg')

# load data
m <- list()
print(lf)
for(fname in lf){
	print(fname)
	resolution = get_parameter_from_fname(fname, 'res')
	clu =read.table(fname, header=T, row.names=1)
	m[[paste0('res', resolution)]]<-clu[,'clusters']
}
# convert to dataframe
m<-do.call("cbind",m)
rownames(m) <- rownames(clu)
head(m)
# run clustree
head_name = opt$infile_prefix
clustree(m, prefix =paste0('res')) + ggtitle(head_name) -> gg


if (!(dir.exists(dirname(opt$outfile)))){
	dir.create(dirname(opt$outfile))
}

# save
ggsave(gg, filename=opt$outfile, height=12,width=16,dev=cairo_pdf)

message("clustree done")
