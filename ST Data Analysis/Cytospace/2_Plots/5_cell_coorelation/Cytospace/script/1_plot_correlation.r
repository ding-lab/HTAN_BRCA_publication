# conda acitvate SPOTlight
library(tidyverse)
library(Seurat)
library(SPOTlight)
library(qs)
library(ggplot2)
library(patchwork)

# Run all the samples available 
object_root = '/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/BRCA_revision/1_RCTD/3_scRNA_by_subtype/4_runCytospace/Cytospace/'


celltype_interests = c('exhausted|Tumor|Macro|NK|DC|B|Mast')
celltype_interests_vector = c("Mono_Macro", "CD4_T exhausted", 
"cDC1", "CD8_T exhausted", "CD8_T preexhausted GZMK+", "DC")
# Load data ---------------------------------------------------------------
for(sample_use in list.files(object_root)){
    object_path = str_glue('{object_root}/{sample_use}/output/5_create_cytospace_seurat/Cytospace_seurat_object.qs')
    # skip if object not exisits
    if (!file.exists(object_path)){
        message(str_glue('Skip {sample_use}'))
        next
    }
    message(str_glue('Processing {sample_use}'))
    st_use = qread(object_path)

    # Extract cell type compostion ---------------------------------------------
    cell_type_composition = st_use@misc$cytospace_result$cell_fractions %>% as.data.frame %>% column_to_rownames('SpotID') %>% as.matrix
    
    # subset cell type to cell type of interests
    cell_keep = intersect(colnames(cell_type_composition), celltype_interests_vector)
    print(cell_keep)
    cell_type_composition = cell_type_composition[, cell_keep]
    message(str_glue('cell_type_composition: {str_c(colnames(cell_type_composition), collapse = ",")}'))


    # Plot correlation analysis -------------------------------------------------
    out_root = '/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/BRCA_revision/1_RCTD/3_scRNA_by_subtype/5_cell_coorelation/Cytospace/out/'
    out_dir = str_glue('{out_root}/{sample_use}/')
    print(out_dir)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    # https://marcelosua.github.io/SPOTlight/articles/SPOTlight_kidney.html
    pdf(str_glue('{out_dir}/1_plotCorrelationMatrix.pdf'))
        plotCorrelationMatrix(cell_type_composition) %>% print()
    dev.off()

    # 2. Co-localization -------------------------------------------------------
    pdf(str_glue('{out_dir}/2a_plotInteractions_heatmap.pdf'))
        p_heat_prop = plotInteractions(cell_type_composition, which = 'heatmap', metric = 'prop') + labs(title =sample_use, subtitle = 'Proportion')
        p_heat_jaccard = plotInteractions(cell_type_composition, which = 'heatmap', metric = 'jaccard') + labs(title = sample_use, subtitle = 'Jaccard')
        print(p_heat_prop)
        print(p_heat_jaccard)
    dev.off()

    # network version
    pdf(str_glue('{out_dir}/2b_plotInteractions_network.pdf'))
        plotInteractions(cell_type_composition, which = 'network', metric = 'prop') 
        plotInteractions(cell_type_composition, which = 'network', metric = 'jaccard')
    dev.off()
}

