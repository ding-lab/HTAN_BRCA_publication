# conda activate seurat4.3
library(optparse)

option_list = list(
    make_option(c("-c", "--cytospace_path"), type="character", default=NULL, help="cytospace result path"),
    make_option(c("-s", "--snRNA_obj"), type="character", default=NULL, help="snRNA object"),
    make_option(c("-S", "--ST_obj"), type="character", default=NULL, help="ST object"),
    make_option(c("-o", "--out_path"), type="character", default=NULL, help="out dir"),
    make_option(c("-n", "--ST_sample_name"), type="character", default=NULL, help="ST sample name"),
    make_option(c("-r", "--script_root"), type="character", default=NULL, help="script root")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# input
path_cytospace_path=opt$cytospace_path
path_snRNA_obj=opt$snRNA_obj
path_ST_obj=opt$ST_obj
script_root=opt$script_root

# output
out_path=opt$out_path
ST_sample_name=opt$ST_sample_name


library(tidyverse)
library(googlesheets4)
library(qs)
library(Seurat)

################################################################################# 
#    Functions                                                            
################################################################################
# Load Cytospace results                                                 
LoadCytospaceResult = function(CytospaceResultPath){
    list(                                      
        assigned_locations = read_csv(str_glue("{CytospaceResultPath}/assigned_locations.csv")),
        cell_counts = read_csv(str_glue("{CytospaceResultPath}/cell_type_assignments_by_spot.csv")),
        cell_fractions = read_csv(str_glue("{CytospaceResultPath}/fractional_abundances_by_spot.csv"))
    )
}

LoadFunction = function(path){
    if(str_detect(path, 'qs')) return(qs::qread(path, nthreads = 30))
    if(str_detect(path, 'rds')) return(readRDS(path))
}

#################################################################################
# Load results
#################################################################################
cyto_result = LoadCytospaceResult(path_cytospace_path)
sn_obj = LoadFunction(path_snRNA_obj)
st_obj = LoadFunction(path_ST_obj)

################################################################################# 
# Create cytospace seurat object
################################################################################
## FUNCTION
source(str_glue('{script_root}/tools/function_generate_cytospace_seurat.r'))
source(str_glue('{script_root}/tools/function_plot_cytospace_spatial.r'))
source(str_glue('{script_root}/tools/function_calculate_spatial_distribution.r'))

## Create cytospace seurat object
message(str_glue("Processing {ST_sample_name}"))
cytoseurat_obj = combine_sn_st_cyto(
    sn_obj,
    st_obj,
    cyto_result
)

## Calculate spatial distribution
cytoseurat_obj = CalculateSpatialDistribution(cytoseurat_obj)

## Replace Sphere with larger number = 90
cytoseurat_obj = CalculateSpatialDistribution(cytoseurat_obj, method = 'Sphere', radius = 90)

## Run SCT 
message("Running SCT on the cytospace seurat object")
cytoseurat_obj = SCTransform(cytoseurat_obj, assay = 'RNA', verbose = FALSE)

## Save cytospace seurat object
#qs::qsave(cytoseurat_obj, str_glue("{out_path}/Cytospace_seurat_{ST_sample_name}.qs"))
# Switch to not use the sample name for now
qs::qsave(cytoseurat_obj, str_glue("{out_path}/Cytospace_seurat_object.qs"))

