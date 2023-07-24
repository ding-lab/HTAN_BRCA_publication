library(optparse)
################################################################################
option_list = list(
    make_option(c("-s", "--sample_name"), type="character", default=NULL, help="sample name"),
    make_option(c("-c", "--cytospace_seurat_path"), type="character", default=NULL, help="cytospace seurat path"),
    make_option(c("-o", "--out_path"), type="character", default=NULL, help="out dir"),
    make_option(c("-l", "--cell_type_column_name"), type="character", default='cell_type_Abbr', help="cell type column name"),
    make_option(c("-p", "--pt_size_factor"), type="numeric", default=1, help="point size factor"),
    make_option(c("-a", "--image_bg_alpha"), type="numeric", default=0.1, help="image background alpha"),
    make_option(c("-t", "--label_st"), type="logical", default=F, help="label ST"),
    make_option(c("-r", "--script_root"), type="character", default=NULL, help="script root")
)

opt_parser = OptionParser(option_list=option_list)
opt=parse_args(opt_parser)

# Assign Arguments
sample_name=opt$sample_name
cytospace_seurat_path=opt$cytospace_seurat_path
out_path=opt$out_path
cell_type_column_name=opt$cell_type_column_name
script_root=opt$script_root
# Plotting parameters
pt_size_factor = opt$pt_size_factor
image_bg_alpha = opt$image_bg_alpha
label_st = opt$label_st


# Expression : snRNA matrix
# Cells : assigned_locations.csv, OriginalCID
# conda activate seurat4.3                                                        
library(tidyverse)
library(googlesheets4)
library(qs)
library(Seurat)
library(patchwork)

source(str_glue('{script_root}/tools/function_plot_cytospace_spatial.r'))

################################################################################
# Load CytospaceSeurat object
################################################################################
## Load cytospace seurat object
cytoobj = qs::qread(cytospace_seurat_path, nthreads = 20)

################################################################################
# Plotting constants
################################################################################
# Get all the cell types
ident_all = cytoobj@meta.data[[cell_type_column_name]] %>% unique()
length(ident_all) ; ident_all

color_use = c(
    RColorBrewer::brewer.pal(8, "Set1"), 
    #RColorBrewer::brewer.pal(8, "Set2"), 
    RColorBrewer::brewer.pal(8, "Set3"),
    RColorBrewer::brewer.pal(8, "Dark2"),
    RColorBrewer::brewer.pal(8, "Paired"),
    RColorBrewer::brewer.pal(8, "Pastel1"),
    RColorBrewer::brewer.pal(8, "Pastel2"),
    RColorBrewer::brewer.pal(8, "Accent"),
    RColorBrewer::brewer.pal(8, "Spectral")
) %>% setNames(ident_all) 


# Select color
ident_exists = cytoobj@meta.data[[cell_type_column_name]] %>% unique() 
ident_color_use = color_use[ident_exists]

################################################################################
# Plot
################################################################################
################################################################################
# PDF version
################################################################################
# Spatial Dim plot
# constant
stroke =  NA

# Make out dir
dir.create(str_glue('{out_path}/'), showWarnings = F)

message(str_glue("Processing {sample_name}"))
message(str_glue("pt_size_factor: {pt_size_factor}"))

# Spatial Dim plot
pdf(str_glue('{out_path}/1_STDim_{cell_type_column_name}.pdf'), w=15, h =7)
    pst_c = SpatialCytoPlot(cytoobj, group.by = cell_type_column_name, label = label_st, label.size = 3, pt.size.factor = pt_size_factor , stroke = stroke, image.alpha = image_bg_alpha, repel = T,
        method = 'Circle') +
        labs(title = str_glue("SpatialDimPlot: Circle")) + 
        scale_fill_manual(values = ident_color_use) + 
        NoLegend()
    pst_s = SpatialCytoPlot(cytoobj, group.by = cell_type_column_name, label = label_st, label.size = 3, pt.size.factor = pt_size_factor , stroke = stroke, image.alpha = image_bg_alpha, repel = T,
        method = 'Sphere') +
        labs(title = str_glue("SpatialDimPlot: Distributed")) + 
        scale_fill_manual(values = ident_color_use) + 
        # Make legend dot larger
        guides(fill = guide_legend(override.aes = list(size=6))) 
    # psn = DimPlot(cytoobj, group.by = cell_type_column_name, label = TRUE, label.size = 3) + 
    #     labs(title = str_glue("UMAP")) 
    #     #scale_color_manual(values = ident_color_use)

    # Also just histology
    p_hist = SpatialCytoPlot(cytoobj, group.by = cell_type_column_name, label = F, label.size = 3, pt.size.factor = 0, stroke = stroke, image.alpha = 1,
        method = 'Circle') +
        labs(title = str_glue("Histology")) + NoLegend()

    ### @@@@ Page 2 #### Plot Umap and other info
    p_umap = DimPlot(cytoobj, group.by = cell_type_column_name, label = TRUE, label.size = 3,  repel = T) + 
        labs(title = str_glue("UMAP")) + 
        scale_color_manual(values = ident_color_use) + 
        coord_fixed()

    ## Plot ALL - With Circle
    p_all = wrap_plots(p_hist, pst_c, p_umap, guides = 'collect', nrow = 1) + #, psn) + 
        plot_annotation(title = str_glue("{sample_name}"), 
            theme = theme(title = element_text(size = 15))) &
            theme(legend.position = 'bottom')
    print(p_all)

    ## Plot ALL - With Sphere
    p_all2 = wrap_plots(p_hist, pst_s, p_umap, guides = 'collect', nrow = 1) + #, psn) + 
        plot_annotation(title = str_glue("{sample_name}"), 
            theme = theme(title = element_text(size = 15))) &
            theme(legend.position = 'bottom')
    print(p_all2)

dev.off()

# Part 2 Make Cell type proportion bar plot
# 2.1 get the proportin count
cell_type_count_df = cytoobj@meta.data %>% count(.data[[cell_type_column_name]]) %>% arrange(desc(n))
# Get cell type proportion order
cell_type_order = cell_type_count_df[[cell_type_column_name]] %>% unique() 

# 2.2 Plot
pdf(str_glue('{out_path}/2_CellTypeProportion_{cell_type_column_name}.pdf'), w=7, h =7)
    p = ggplot(cell_type_count_df, aes(x = factor(.data[[cell_type_column_name]], levels = rev(cell_type_order)), y = n)) + 
        geom_bar(stat = 'identity', fill = ident_color_use) + 
        labs(title = str_glue("Cell type proportion"), x=cell_type_column_name) + 
        cowplot::theme_cowplot() + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
        coord_flip()
    print(p)
dev.off()


message(str_glue("Done plotting {sample_name}") )
