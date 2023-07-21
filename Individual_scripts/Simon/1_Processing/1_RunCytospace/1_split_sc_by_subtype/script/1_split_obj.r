library(tidyverse)
library(Seurat)
library(googlesheets4)
library(qs)

# load object
obj = readRDS("/PATH/TO/BRCA/scRNA/all_merged/scRNA_merged.rds")

# Add molecular subtype. Load PAM50
sheet_brca = read_tsv('../0_PAM50_table/HTAN_BRCA_PAM50.txt')

# save a copy
todaydata = format(Sys.Date(), "%Y%m%d")
analysis_root = '/ROOT/FOLDER/'
write_tsv(sheet_brca, str_glue('{analysis_root}/table/googlesheet_metadata_backup_{todaydata}.tsv'))

# 1. Filer out unknown
celltype_column = 'old_cell_type_specific'
obj_filtered = subset(obj, cells = obj@meta.data %>% filter(!str_detect(.data[[celltype_column]], 'Unknown')) %>% rownames)
# check 
obj_filtered@meta.data %>% count(.data[[celltype_column]])

# 2. Add PAM50 subtype
subtype_meta = sheet_brca[, c('Case', 'PAM50_bulk','PAM50_sn')] %>%
    mutate(PAM50_final = ifelse(is.na(PAM50_sn), PAM50_bulk, PAM50_sn)) %>%
	# filter out NA and 'Normal' since they doen't contribute to determining PAM50 subtype
	filter(!is.na(PAM50_final) & !PAM50_final %in% c('NA', 'Normal' )) %>% 
    select(Case, PAM50_final) %>% 
	distinct()

new_meta = left_join(obj_filtered@meta.data %>% rownames_to_column(), subtype_meta, by = c('Sample' = 'Case')) %>% column_to_rownames()

# 3. Split and save object based on the rownames for each subtype
# Split intrinsicly ignore the NA 
new_meta_list_by_subtype = new_meta %>% split(., .$PAM50_final) %>% glimpse

# 4. Save 
dir.create(str_glue("{analysis_root}/out/{todaydata}/"), showWarnings = FALSE, recursive = TRUE)
iwalk(new_meta_list_by_subtype, function(df, subtype){
	message("Processing: ", subtype)
	newobj = subset(obj_filtered, cells = rownames(df))
	qs::qsave(newobj, str_glue("{analysis_root}/out/{todaydata}/scRNA_merged_{subtype}.qs"), nthreads = 2)
})
