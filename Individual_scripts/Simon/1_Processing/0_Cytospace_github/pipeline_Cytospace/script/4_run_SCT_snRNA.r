library(optparse)

opt_list = list(
    make_option(c("-s", "--sn_file_path"), type="character", default=NULL, help="snRNA object"),
    make_option(c("-n", "--sn_sample_name"), type="character", default=NULL, help="snRNA name"),
    make_option(c("-o", "--out_path"), type="character", default=NULL, help="out dir")
)

opt_parser = OptionParser(option_list=opt_list)
opt=parse_args(opt_parser)

# Assign Arguments
sn_object=opt$sn_file_path
sn_name=opt$sn_sample_name
out_dir=opt$out_path

# Expression : snRNA matrix
# Cells : assigned_locations.csv, OriginalCID
# conda activate seurat4.3                                                        
library(tidyverse)
library(qs)
library(Seurat)

#############################################################
# Load snRNA object
sn = qs::qread(sn_object)

################################################################################
# IF sn object missing SCT assay, calculate it. else just save it
################################################################################
# first ch if has SCT assay 
if('SCT' %in% Assays(sn)){
    message(str_glue("SCT assay already exists in {sn_object}. Saving as qs and exit..."))
    qs::qsave(sn, str_glue("{out_dir}/{sn_name}.qs"), nthreads = 30)
}else{
    message(str_glue("SCT assay does not exist in {sn_object}. Calculating SCT..."))
    # SCT
    sn = SCTransform(sn, assay = 'RNA', verbose = FALSE)
    # Run Dim reduction - PCA and UMAP
    sn <- RunPCA(sn, assay = "SCT", verbose = FALSE)
    sn <- RunUMAP(sn, reduction = "pca", dims = 1:30)
    # Save snRNA object
    print(str_glue("saving to location: {out_dir}/{sn_name}.qs"))
    qs::qsave(sn, str_glue("{out_dir}/{sn_name}.qs"), nthreads = 30)
}