# conda activate seurat4.3
library(optparse)

# Options
opt_list = list(
    make_option(c("-n", "--sn_sample_name"), type="character", default=NULL, help="sn sample name"),
    make_option(c("-f", "--sn_file_path"), type="character", default=NULL, help="sn file path"),
    make_option(c("-c", "--cell_type_column_name"), type="character", default=NULL, help="cell type column name"),
    make_option(c("-o", "--out_path"), type="character", default=NULL, help="out path"),
    make_option(c("-s", "--script_root"), type="character", default=NULL, help="script root")
)

opt_parser = OptionParser(option_list=opt_list)
opt = parse_args(opt_parser)

# Assign paramters
out_path = opt$out_path
cell_type_column_name = opt$cell_type_column_name
sn_sample_name = opt$sn_sample_name
sn_file_path = opt$sn_file_path
script_root = opt$script_root

# Main library
library(tidyverse)
library(googlesheets4)
library(qs)
library(Seurat)

# process function 
source(str_glue('{script_root}/tools/generate_cytospace_from_seurat_object.R'))

# # load sn object
sn = if(str_detect(sn_file_path, '.rds')){
    message(str_c('Loading snRNA object: ', sn_file_path))
    readRDS(sn_file_path)
}else{
    message(str_c('Loading snRNA object: ', sn_file_path))
    qs::qread(sn_file_path, nthreads = 30)
}

## New process - use my own code
message(str_c('Processing snRNA sample: ', sn_sample_name))

# 1. Subset out cell_type_Abbr = NA spot
message(str_c('Number of cells original: ', ncol(sn)))
cell_keep = sn@meta.data %>% filter(!is.na(.data[[cell_type_column_name]])) %>% rownames()
message(str_c('Number of cells keep: ', length(cell_keep)))
sn_narm = subset(sn, cells = cell_keep)

# 2. Save cell_type_labels.txt
cell_type_labels = sn_narm@meta.data[,cell_type_column_name, drop=F] %>% 
    rownames_to_column('Cell IDs') %>%
    setNames(c('Cell IDs', 'CellType')) %>%
    select(c('Cell IDs', 'CellType')) 

write_tsv(cell_type_labels, str_glue('{out_path}/cell_type_labels.txt'))

# 3. Save scRNA_data.txt
count_mtx = sn_narm@assays$RNA %>% GetAssayData(slot = 'counts')
message(str_c('Matrix size: ', toString(dim(count_mtx))))
count_df = count_mtx %>% as_tibble %>%
    mutate(GENES = rownames(count_mtx)) %>%
    select(GENES, everything()) 

write_tsv(count_df, str_glue('{out_path}/scRNA_data.txt'))
