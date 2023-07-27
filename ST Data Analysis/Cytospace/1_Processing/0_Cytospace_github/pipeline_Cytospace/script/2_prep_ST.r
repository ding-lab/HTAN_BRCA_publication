# conda activate seurat4.3
# Prepare ST files for CytoSpace

library(optparse)
opt_list = list(
    make_option(c("-n", "--sample_name"), type="character", default=NULL, help="sample name"),
    make_option(c("-f", "--ST_file_path"), type="character", default=NULL, help="ST file path"),
    make_option(c("-o", "--st_out_dir"), type="character", default=NULL, help="ST out dir"),
    make_option(c("-s", "--script_root"), type="character", default=NULL, help="script root")
)

opt_parser = OptionParser(option_list=opt_list)
opt=parse_args(opt_parser)

st_out_dir = opt$st_out_dir
ST_file_path = opt$ST_file_path
sample_name = opt$sample_name
script_root = opt$script_root

# Generate required file from seurat objects
library(tidyverse)
library(googlesheets4)
library(qs)
library(Seurat)

# Process function
source(str_glue('{script_root}/tools/generate_cytospace_from_seurat_object.R'))

# Load ST
# If ST_Sex_split_obj not NA, load from there. otherwise load from ST_obj
st = if(str_detect(ST_file_path, '.rds')){
    message(str_c('Loading ST object using readRDS: ', ST_file_path))
    readRDS(ST_file_path)
}else{
    message(str_c('Loading ST object using qread: ', ST_file_path))
    qs::qread(ST_file_path, nthreads = 30)
}


# Run function to load generate_cytospace_from_ST_seurat_object
message(str_c('Processing ST sample: ', sample_name))
dir.create(str_c(st_out_dir,"/"))
generate_cytospace_from_ST_seurat_object(
    st_seurat = st, 
    dir_out = str_c(st_out_dir,"/"),
    slice = Images(st)[[1]] # First slice 
)
