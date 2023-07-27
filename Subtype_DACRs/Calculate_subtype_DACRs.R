library(Seurat)
library(Signac)
library(plyr)
library(dplyr)

library(future)

plan("multiprocess", workers =5)
options(future.globals.maxSize = 100 * 1024^3) # for 200 Gb RAM


path_to_ATAC_merged_obj=''
Piece_ID_subtype_annotation=''

#Add annotation of normal duct cells for these cluters:
#13 - Lum mature
#15 - Basal prog
#30 - Lum prog

ATAC=readRDS(path_to_ATAC_merged_obj)

#Use cell_type_manual for cell type annotation

ATAC$cell_type_manual=ifelse(ATAC$seurat_clusters==13,"Lum mature",ATAC$cell_type_manual)
ATAC$cell_type_manual=ifelse(ATAC$seurat_clusters==15,"Basal progenitors",ATAC$cell_type_manual)
ATAC$cell_type_manual=ifelse(ATAC$seurat_clusters==30,"Lum progenitors",ATAC$cell_type_manual)
atac=ATAC

ATAC=subset(atac, cell_type_manual %in% c('Tumor'))


annot=read.table(Piece_ID_subtype_annotation,sep='\t',header=T)
lum=unique(annot$Sample_ID[annot$PAM50_sn %in% c('LumA','LumB')])
basal=unique(annot$Sample_ID[annot$PAM50_sn %in% c('Basal')])
her2=unique(annot$Sample_ID[annot$PAM50_sn %in% c('Her2')])


#Add fragment counts in peaks for the latest peaks list
peak.data <- GetAssayData(object = ATAC, assay = 'peaksMACS2', slot = "counts")
peak.counts <- colSums(x = peak.data)
ATAC <- AddMetaData(object = ATAC, metadata = peak.counts, col.name = 'peak_RF_500MACS2')
DefaultAssay(ATAC)='peaksMACS2'


ATAC$test='Other'
ATAC$test=ifelse(ATAC$Piece_ID %in% basal, 'Basal', ATAC$test)
ATAC$test=ifelse(ATAC$Piece_ID %in% lum, 'Lum', ATAC$test)
ATAC$test=ifelse(ATAC$Piece_ID %in% her2, 'Her2', ATAC$test)
atac=ATAC

#Remove cells annotated as tumor from NATs:
ATAC=subset(atac, test!='Other')

Idents(ATAC)=ATAC$test

cell_groups=unique(ATAC$test)
all_da_peaks=NULL
for (cell_group in cell_groups){

da_peaks <- FindMarkers(
  object = ATAC,
  ident.1 = cell_group,
  only.pos = FALSE,
  min.pct = 0.1,
  min.diff.pct=0,
  logfc.threshold=0,
  test.use = 'LR',
  latent.vars = 'peak_RF_500MACS2'
)

da_peaks$cell_type=cell_group
da_peaks$peak=rownames(da_peaks)
all_da_peaks=rbind(all_da_peaks,da_peaks)
print(cell_group)
}

write.table(all_da_peaks, paste("out/Subtype_DACRs_TumorCells.minPct0.1.tsv", sep=""),sep="\t",quote=FALSE,row.names=FALSE)
