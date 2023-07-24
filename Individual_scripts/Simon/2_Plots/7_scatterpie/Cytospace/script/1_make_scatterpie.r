# conda acitvate SPOTlight
library(tidyverse)
library(Seurat)
library(SPOTlight)
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
analysis_root = '/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/BRCA_revision/1_RCTD/3_scRNA_by_subtype/7_scatterpie/Cytospace/'

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

# Get all cell types and assign colors
col_celltype = c(
    # Lymphoid https://www.khanacademy.org/science/biology/developmental-biology/development-and-differentiation/a/introduction-to-development
    # Innate
    "NK cells"="#ffcc24",
    "NKT_cells"="#ffff6d",
    setNames(c("#920000", "#924900", "#db6d00", "#ff6d24"), c("cDC1", "cDC2", "DC", "pDC")), 
    # adaptive
    setNames(c("#ff6db6", "#ffa6ab"), c("B", "Plasma")),
    "Treg"="#8899cc",
    setNames(c("#001166","#006ddb","#6db6ff", "#b6dbff"), c("CD4_T", "CD4_T activated", "CD4_T exhausted", "CD4_T RPhigh")),
    setNames(c("#38774a", "#43bd62", "#c3f266", "#85c267", "#46801b"),c("CD8_CTL", "CD8_T", "CD8_T proliferating", "CD8_T preexhausted GZMK+","CD8_T exhausted")),
    # Myeloid
    "Mast"="#ff00ee",
    "Mono_Macro"="#ff3324",
    setNames(c("#ffccff", "#ff99ff", "#ff66ff", "#ff33ff"),c("CAF", "cCAF", "mCAF", "vCAF")),
    "Endothelial"="#ffd3bd",
    "Normal_duct"="#224949",
    "Tumor"="#112211"
)
############################################################################################################
## Function
############################################################################################################
# Plot function: Plot whole pie chart plot add to each spatilpie
MakeSampleLevelPieFromMatrix = function(mtx, percent_cutoff = 0.01){
    sample_pie_df = mtx %>% apply(2, mean) %>% as.data.frame %>% rownames_to_column() %>% setNames(c('celltype','percent'))
    sample_pie_df2 <- sample_pie_df %>% 
        mutate(csum = rev(cumsum(rev(percent))), 
         pos = percent/2 + lead(csum, 1),
         pos = if_else(is.na(pos), percent/2, pos))
    
    sample_pie_df2 =  sample_pie_df2 %>% mutate(label_celltype = ifelse(percent > percent_cutoff, celltype, ""))
    
    p_pie_sample = ggplot(sample_pie_df, aes(x="", y = percent, fill = celltype)) + 
        geom_bar(stat="identity", width=1) + 
        coord_polar(theta = "y") +
        # label
        geom_text(aes(label = str_c(round(percent, 2) * 100, '%')), position = position_stack(vjust = 0.5), color = "white") +
        #theme_void() + 
        scale_y_continuous(breaks = sample_pie_df2$pos, labels = sample_pie_df2$label_celltype) +
        theme(axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(size = 15), 
        legend.position = "none", # Removes the legend
        panel.background = element_rect(fill = "white")) +
        scale_fill_manual(
                        values = pal,
                        breaks = names(pal)) 
    return(p_pie_sample)
}

# Plot function: Plot whole pie chart plot add to each spatilpie
MakeSampleLevelPieFromMatrix = function(mtx, percent_cutoff = 0.01){
    sample_pie_df = mtx %>% apply(2, mean) %>% as.data.frame %>% rownames_to_column() %>% setNames(c('celltype','percent'))
    sample_pie_df2 <- sample_pie_df %>% 
        mutate(csum = rev(cumsum(rev(percent))), 
         pos = percent/2 + lead(csum, 1),
         pos = if_else(is.na(pos), percent/2, pos))
    sample_pie_df2 =  sample_pie_df2 %>% mutate(label_celltype = ifelse(percent > percent_cutoff, celltype, ""))
    
    p_pie_sample = ggplot(sample_pie_df2, aes(x="", y = percent, fill = fct_inorder(celltype))) + 
        geom_bar(stat="identity", width=1) + 
        coord_polar(theta = "y") +
        # label
        ggrepel::geom_label_repel(data = sample_pie_df2,
                   aes(y = pos, label = paste0(round(percent,2)*100, "%")),
                   size = 4.5, nudge_x = 1, show.legend = FALSE) + 
        #theme_void() + 
        scale_y_continuous(breaks = sample_pie_df2$pos, labels = sample_pie_df2$label_celltype) +
        theme(axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(size = 15), 
        legend.position = "none", # Removes the legend
        panel.background = element_rect(fill = "white")) +
        scale_fill_manual(
                        values = pal,
                        breaks = names(pal)) 
    return(p_pie_sample)
}
    


############################################################################################################
# Plot cell type composition scatterpie ------------------------------------
############################################################################################################
# Load data ---------------------------------------------------------------
# Load and Plot data ---------------------------------------------------------------
out_root = '/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/BRCA_revision/1_RCTD/3_scRNA_by_subtype/7_scatterpie/Cytospace/out/'
rctd_seruat_root = '/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/BRCA_revision/1_RCTD/3_scRNA_by_subtype/3_runRCTD/4_Add_to_Seurat/out/'
selected_cell_list = list(
    basal_enriched = c('CD8_T proliferating', 'CD8_T exhausted', 'CD4_T exhausted', 'CD4_T activated', 'cDC2','Treg','CD4_T'),
    cd8T = c('CD8_T proliferating', 'CD8_T exhausted', 'CD8_T preexhausted GZMK+', 'CD8_T'),
    cd8Texhausted = c('CD8_T','CD8_T exhausted')
)


# parameter
min_cell_cutoff= 0.1
sample_uses = list.files(object_root)

for(sample_use in samples_use){
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
    
    
    # Extract cell type compostion ---------------------------------------------
    cell_type_composition = st_use@misc$cytospace_result$cell_fractions %>% as.data.frame %>% column_to_rownames('SpotID') %>% as.matrix

    # Make sure size of the st_rctd_seruat == size of the cell_type_composition
    st_rctd_seurat = subset(st_rctd_seurat, cells = rownames(cell_type_composition))
    
    # skip sample if the ncol of st_rctd_seruat != nrow of cell_type_composition
    if(ncol(st_rctd_seurat) != nrow(cell_type_composition)){
        message(str_glue('Skip {sample_use}. ncol of st_rctd_seruat != nrow of cell_type_composition'))
        next
    }

    # Plot scatterpie -------------------------------------------------
    mat = cell_type_composition
    ct <- colnames(mat) %>% sort
    mat[mat < min_cell_cutoff] <- 0

    # Define color palette
    pal = col_celltype[intersect(names(col_celltype), ct)]

    # Plot
    out_dir = str_glue('{out_root}/{sample_use}/')
    print(out_dir)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    p_scatter = plotSpatialScatterpie(
            x = st_rctd_seurat, # need ot use this since Cytospace uses different coordinates
            y = mat,
            cell_types = colnames(mat),
            img = FALSE,
            scatterpie_alpha = 1,
            pie_scale = 0.35) +
            scale_fill_manual(
                values = pal,
                breaks = names(pal)) 
    
    # Adding sample level pie
    p_sample_pie = MakeSampleLevelPieFromMatrix(mat, percent_cutoff = 0.05) 

    pdf(str_glue('{out_dir}/1_plotSpatialScatterpie.pdf'))
    print(p_scatter)
    print(p_sample_pie )
    dev.off()

    ## -----------------------------------------------------------------------
    # Plot cell type composition scatterpie ------------------------------------
    # Of select cells 
    # Plot scatterpie -------------------------------------------------
    iwalk(selected_cell_list, function(celltype_interests_vector, plt_list_name){
        celltype_interests_vector = intersect(colnames(cell_type_composition), celltype_interests_vector)
        message(str_glue('Processing {sample_use} {plt_list_name}'))
        mat = cell_type_composition
        mat[mat < min_cell_cutoff] <- 0

        # Define color palette
        pal = col_celltype[intersect(names(col_celltype), celltype_interests_vector)]

        # Plot
        p_scatter = plotSpatialScatterpie(
                x = st_rctd_seurat, # need ot use this since Cytospace uses different coordinates
                y = mat,
                cell_types = celltype_interests_vector,
                img = TRUE,
                scatterpie_alpha = 1,
                pie_scale = 0.35) +
                scale_fill_manual(
                    values = pal,
                    breaks = names(pal)) 

        # Adding sample level pie
        p_sample_pie = MakeSampleLevelPieFromMatrix(mat[,celltype_interests_vector], percent_cutoff = 0.01) 

        out_dir = str_glue('{out_root}/{sample_use}/')
        print(out_dir)
        dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
        pdf(str_glue('{out_dir}/2_plotSpatialScatterpie_selected_{plt_list_name}.pdf'))
        print(p_scatter)
        print(p_sample_pie)
        dev.off()
    })
}

