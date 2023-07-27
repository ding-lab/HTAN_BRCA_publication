################################################
# Function to create Cytospace Seurat object
################################################

## Combine snRNA and ST image with cytospace result
combine_sn_st_cyto <- function(sn_use, st_use, cyto_result) {
    # 0. extract relevant objects from inputs
    loc_use <- cyto_result$assigned_locations # Location files 
    img_name = Images(st_use)[[1]]# Image name, default use the first one
    st_img_use <- st_use@images[[img_name]] # Image
    st_img_use_coords <- st_img_use@coordinates # Image Coordinate
    
    
    # 1. Generate new coordinate table
    coor_df_new_all <- loc_use %>% 
        left_join(
        st_img_use_coords, by = c('row', 'col')
        ) %>% 
        column_to_rownames('UniqueCID') 
    
    coor_df_new <- coor_df_new_all %>%     
        select(c('tissue','row','col','imagerow','imagecol'))
    
    # 1.1 Replace with new image coordinate
    cyto_img <- st_img_use
    cyto_img@coordinates = coor_df_new
    
    # 1.2 Create a OriginalCID -> UniqueCID map
    loc_use_map <- loc_use %>% 
        group_by(OriginalCID) %>% 
        slice(1) %>% 
        select(OriginalCID, UniqueCID) %>% 
        ungroup()
    
    snRNAid_UniqudCID_mapping <- setNames(loc_use_map$UniqueCID, loc_use_map$OriginalCID)
    
    # 2. Replace image coordinate with new cell coordinate
    st_img_new <- st_img_use
    st_img_new@coordinates = coor_df_new
    
    # 3. Make NEW snRNA object Manually
    sn_use_new_RNA <- sn_use %>% 
        GetAssay('RNA') %>% 
        subset(cells = loc_use$OriginalCID) %>%
        RenameCells(loc_use$UniqueCID)
    
    sn_use_new_Meta <- sn_use@meta.data[loc_use$OriginalCID, ] %>% 
        rownames_to_column('OriginalCellID') %>%
        `rownames<-`(loc_use$UniqueCID) 
    
    # # 3.1 Create new Seurat object
    message('Creating new Seurat object')
    st_new2 <- CreateSeuratObject(counts = sn_use_new_RNA, assay = 'RNA', meta.data = sn_use_new_Meta)
    # Add SCT
    if('SCT' %in% Assays(sn_use)) {
      sn_use_new_SCT <- sn_use %>% 
        GetAssay('SCT') %>% 
        subset(cells = loc_use$OriginalCID) %>%
        RenameCells(loc_use$UniqueCID)
      st_new2[['SCT']] <- sn_use_new_SCT
    }

    # Add image object to new object
    message('Adding image to the new object')
    DefaultAssay(cyto_img) <- 'RNA'
    Key(cyto_img) = img_name
    st_new2[[img_name]] <- cyto_img

    # 3.3 Extract relevant information from snRNA object
    # 3.4 Replace new embedding into the reduction slot
    reductions_new <- sn_use@reductions
    if('pca' %in% Reductions(sn_use)) {
      # PCA
      pca_embadding_new <- sn_use@reductions$pca@cell.embeddings[loc_use$OriginalCID, ] %>% `rownames<-`(loc_use$UniqueCID)
      reductions_new$pca@cell.embeddings = pca_embadding_new
    }
    # UMAP
    if('umap' %in% Reductions(sn_use)) {
      umap_embadding_new <- sn_use@reductions$umap@cell.embeddings[loc_use$OriginalCID, ] %>% `rownames<-`(loc_use$UniqueCID)
      reductions_new$umap@cell.embeddings = umap_embadding_new
    }
    st_new2@reductions = reductions_new
    
    
    # 3.5 save Cytospace result to misc slot
    message('Saving cytospace result to misc slot @misc$cytospace_result')
    st_new2@misc$cytospace_result = cyto_result
  
  return(st_new2)
}
