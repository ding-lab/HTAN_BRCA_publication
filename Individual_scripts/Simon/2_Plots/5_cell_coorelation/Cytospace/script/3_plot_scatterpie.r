# conda acitvate SPOTlight
library(tidyverse)
library(Seurat)
library(SPOTlight)
library(qs)
library(ggplot2)
library(patchwork)

# Run all the samples available 
object_root = '/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/BRCA_revision/1_RCTD/3_scRNA_by_subtype/4_runCytospace/Cytospace/'
rctd_seruat_root = '/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/BRCA_revision/1_RCTD/3_scRNA_by_subtype/3_runRCTD/4_Add_to_Seurat/out/'

# Parameters
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
    st_use = qread(object_path, nthreads = 30)
    # skip if the object is not at the path
    if(!file.exists(str_glue('{rctd_seruat_root}/{sample_use}/ST_Seurat_RCTD_old_cell_type_specific.qs'))){
        message(str_glue('Skip {sample_use}. seruat object not exists'))
        next
    }
    st_rctd_seurat = qread(str_glue('{rctd_seruat_root}/{sample_use}/ST_Seurat_RCTD_old_cell_type_specific.qs', nthreads = 30))
    
    print(dim(st_rctd_seurat))
    print(dim( st_use@misc$cytospace_result$cell_fractions))
    # Extract cell type compostion ---------------------------------------------
    cell_type_composition = st_use@misc$cytospace_result$cell_fractions %>% as.data.frame %>% column_to_rownames('SpotID') %>% as.matrix


   # subset cell type to cell type of interests
    #celltype_keep = str_subset(colnames(cell_type_composition), celltype_interests)
    #cell_type_composition = cell_type_composition[, celltype_keep]
    cell_type_composition = cell_type_composition[, intersect(colnames(cell_type_composition), celltype_interests_vector)]
    message(str_glue('cell_type_composition: {str_c(colnames(cell_type_composition), collapse = ",")}'))

    # skip sample if the ncol of st_rctd_seruat != nrow of cell_type_composition
    if(ncol(st_rctd_seurat) != nrow(cell_type_composition)){
        message(str_glue('Skip {sample_use}. ncol of st_rctd_seruat != nrow of cell_type_composition'))
        next
    }

    # Plot scatterpie -------------------------------------------------
    mat = cell_type_composition
    ct <- colnames(mat)
    mat[mat < 0.1] <- 0

    # Define color palette
    # (here we use 'paletteMartin' from the 'colorBlindness' package)
    paletteMartin <- c(
        "#000000", "#004949", "#009292", "#ff6db6", "#ffb6db", 
        "#490092", "#006ddb", "#b66dff", "#6db6ff", "#b6dbff", 
        "#920000", "#924900", "#db6d00", "#24ff24", "#ffff6d")

    pal <- colorRampPalette(paletteMartin)(length(ct))
    names(pal) <- ct

    # Plot
    out_root = '/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/BRCA_revision/1_RCTD/3_scRNA_by_subtype/5_cell_coorelation/Cytospace/out/'
    out_dir = str_glue('{out_root}/{sample_use}/')
    print(out_dir)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    pdf(str_glue('{out_dir}/3_plotSpatialScatterpie.pdf'))
    p_scatter = plotSpatialScatterpie(
            x = st_rctd_seurat, # need ot use this since Cytospace uses different coordinates
            y = mat,
            cell_types = colnames(mat),
            img = FALSE,
            scatterpie_alpha = 1,
            pie_scale = 0.4) +
            scale_fill_manual(
                values = pal,
                breaks = names(pal)) 
    print(p_scatter)
    dev.off()
}
