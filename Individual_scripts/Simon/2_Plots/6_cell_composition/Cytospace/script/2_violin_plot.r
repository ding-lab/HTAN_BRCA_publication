# conda acitvate seurat4.3
library(tidyverse)
library(Seurat)
library(qs)
library(ggplot2)
library(patchwork)
library(googlesheets4)

### Preprocess metadata -----------------------------------------------------
# Add molecular subtype
# load google sheet
gs4_deauth()
sheet_brca = read_sheet('https://docs.google.com/spreadsheets/d/10wyXuZaGAhKx0EWDI4ipCw8qgjTWS3uN9U_mScTS2DA/edit#gid=0')
todaydata = format(Sys.Date(), "%Y%m%d")
analysis_root = '/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/BRCA_revision/1_RCTD/3_scRNA_by_subtype/6_cell_composition/Cytospace/'

# 1. Filer out unknown
celltype_column = 'old_cell_type_specific'

# 2. Add PAM50 subtype
subtype_meta = sheet_brca[, c('Case', 'FINAL_CALL')] %>% 
	distinct()



# Run all the samples available 
object_root = '/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/BRCA_revision/1_RCTD/3_scRNA_by_subtype/4_runCytospace/Cytospace/'


celltype_interests = c('exhausted|Tumor|Macro|NK|DC|B|Mast')
celltype_interests_vector = c("Mono_Macro", "CD4_T exhausted", 
"cDC1", "CD8_T exhausted", "CD8_T preexhausted GZMK+", "DC")
# Load data ---------------------------------------------------------------
celltype_list = map(list.files(object_root), function(sample_use){
    object_path = str_glue('{object_root}/{sample_use}/output/5_create_cytospace_seurat/Cytospace_seurat_object.qs')
    # skip if object not exisits
    if (!file.exists(object_path)){
        message(str_glue('Skip {sample_use}'))
        return(NULL)
    }
    message(str_glue('Processing {sample_use}'))
    st_use = qread(object_path, nthreads = 50)

    # Extract cell type compostion ---------------------------------------------
    cell_type_composition = st_use@misc$cytospace_result$cell_fractions %>% as.data.frame %>%
        mutate(SampleID = sample_use) 

    return(cell_type_composition)
}) %>% setNames(list.files(object_root))

# Remove NULL
celltype_list_new = celltype_list[!map(celltype_list, is.null) %>% unlist]
# Get sample level cell type composotion 
celltype_composition_df = imap(celltype_list_new, function(df, sample_id){
    data.frame(
        celltype_ratio = df %>% select(-SampleID, -SpotID) %>% as.matrix %>% apply(2, function(x) mean(x,na.rm = T)),
        sample_id = sample_id
        ) %>% 
    rownames_to_column('celltype') 
}) %>% bind_rows()

# Add Cancer type 
celltype_composition_df = celltype_composition_df %>% 
    # Extract Case ID
    mutate(Case = str_extract(sample_id, 'HT[0-9]{3}[A-Z][0-9]')) %>%
    # Left join the subtype_meta by Case ID
    left_join(subtype_meta, by = 'Case')

# Save the cell type composition
celltype_composition_df %>% write_csv(str_glue('{analysis_root}/out/cell_type_composition.csv'))

# Color palette
color_pam50 = c("Basal"="#FF0000",'Her2'='#FF69B4','LumA'='#00008B','LumB'='#ADD8E6')

# Plot cell type composition violin plot ------------------------------------
# stat
my_comparisons <- list( c("Basal", "LumA"), c("Basal", "LumB"))

############################################################################################################
# Make violin plots
############################################################################################################
# 1. all cell types ---------------------------------------------------------
# Per cancer type
library(ggpubr)
pdf(str_glue('{analysis_root}/out/cell_type_composition_violin.pdf'), w = 15, h = 15)
p_cell_comp = ggplot(celltype_composition_df, aes(x = FINAL_CALL, y = celltype_ratio, fill = FINAL_CALL)) +
    geom_violin(alpha = 0.6, scale = 'width')  +
    geom_jitter(aes(color = FINAL_CALL)) +
    geom_boxplot(width = 0.1, fill = 'gray90') + 
    stat_compare_means(comparisons = my_comparisons, method = 't.test', label = 'p.format', hide.ns = T, ref.group = '.all.') + 
    facet_wrap(~celltype, scales = 'free_y') +
    theme_classic2() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1)#,
        #legend.position = 'none'
    ) +
    labs(
        x = 'PAM50 subtype',
        y = 'Cell type composition',
        title = 'Cell type composition by PAM50 subtype'
    ) + 
    scale_y_continuous(expand = expansion(mult = c(0,0.2))) + 
    scale_fill_manual(values = color_pam50) +
    scale_color_manual(values = color_pam50)
p_cell_comp
dev.off()


# 2. all cell types ---------------------------------------------------------
# But plot them seperatedly and combine them together
# Per cancer type - V2 
library(ggpubr)
library(patchwork)
celltype_vect = celltype_composition_df$celltype %>% unique %>% sort
p_vlin_list = map(celltype_composition_df$celltype %>% unique, function(CellType){
    df_use = celltype_composition_df %>% filter(celltype == CellType)
    ggplot(df_use, aes(x = FINAL_CALL, y = celltype_ratio, fill = FINAL_CALL)) +
        geom_violin(alpha = 0.6, scale = 'width') +
        geom_jitter(aes(color = FINAL_CALL)) +
        geom_boxplot(width = 0.1, fill = 'gray90') + 
        stat_compare_means(comparisons = my_comparisons, method = 't.test', label = 'p.format', hide.ns = T, ref.group = '.all.') + 
        theme_classic2() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1)#,
            #legend.position = 'none'
        ) +
        labs(
            x = 'PAM50 subtype',
            y = 'Cell type composition',
            title = CellType
        ) + 
    scale_y_continuous(expand = expansion(mult = c(0,0.2)))  + 
    scale_fill_manual(values = color_pam50) +
    scale_color_manual(values = color_pam50)
}) %>% setNames(celltype_vect)

pdf(str_glue('{analysis_root}/out/cell_type_composition_violin_v2.pdf'), w = 15, h = 15)
wrap_plots(p_vlin_list, guides = 'collect')
dev.off()


# 3. Basal related cell type violin plot ---------------------------------------------------------
# Select Basal related cell type violin plot
selected_cell_basal = c('CD8_T proliferating', 'CD8_T exhausted', 'CD4_T exhausted', 'CD4_T activated', 'cDC2','Treg','CD4_T')
library(ggpubr)
pdf(str_glue('{analysis_root}/out/cell_type_composition_violin_Basal_selected.pdf'), w = 12, h = 8)
celltype_composition_df %>% 
    filter(celltype %in% selected_cell_basal) %>%
    mutate(celltype = fct_reorder(celltype, celltype_ratio)) %>%
ggplot(aes(x = FINAL_CALL, y = celltype_ratio, fill = FINAL_CALL)) +
    geom_violin(alpha = 0.6, scale = 'width')  +
    geom_jitter(aes(color = FINAL_CALL)) +
    geom_boxplot(width = 0.1, fill = 'gray90') + 
    stat_compare_means(comparisons = my_comparisons, method = 't.test', label = 'p.format', hide.ns = T, ref.group = '.all.') + 
    facet_wrap(~celltype, scales = 'free_y', ncol = 4) +
    theme_classic2() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1)#,
        #legend.position = 'none'
    ) +
    labs(
        x = 'PAM50 subtype',
        y = 'Cell type composition',
        title = 'Cell type composition by PAM50 subtype - Basal Enriched'
    ) + 
    scale_y_continuous(expand = expansion(mult = c(0,0.2))) +
    scale_fill_manual(values = color_pam50) +
    scale_color_manual(values = color_pam50)
dev.off()
