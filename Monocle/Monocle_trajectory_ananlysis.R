library(monocle)
library(Signac)
library(Seurat)
library(RColorBrewer)
library(ggplot2)


n_components=10
sample='MergedObj'
n_cells=1000

cols=c(brewer.pal(n=3, name='Dark2'),'#FF69B4','#FF0000','#00008B','#ADD8E6')
names(cols)=c('Lum progenitors','Basal progenitors','Lum mature','Her2','Basal','LumA','LumB')

path_to_ATAC_merged_epi_obj=''
Piece_ID_subtype_annotation=''

#Use already subsetted object with epithelial cells:
ATAC=readRDS(path_to_ATAC_merged_epi_obj)


ATAC$cell_type_manual=ifelse(ATAC$seurat_clusters==13,"Lum mature",ATAC$cell_type_manual)
ATAC$cell_type_manual=ifelse(ATAC$seurat_clusters==15,"Basal progenitors",ATAC$cell_type_manual)
ATAC$cell_type_manual=ifelse(ATAC$seurat_clusters==30,"Lum progenitors",ATAC$cell_type_manual)

annot=read.table(Piece_ID_subtype_annotation, sep='\t',header=T)
luma=unique(annot$Piece_ID[annot$PAM50_sn %in% c('LumA')])
lumb=unique(annot$Piece_ID[annot$PAM50_sn %in% c('LumB')])
basal=unique(annot$Piece_ID[annot$PAM50_sn %in% c('Basal')])
her2=unique(annot$Piece_ID[annot$PAM50_sn %in% c('Her2')])

ATAC$test=ifelse(ATAC$cell_type_manual=='Tumor',ATAC$Piece_ID, ATAC$cell_type_manual)

Idents(ATAC)=ATAC$test

#Downsample object:
ATAC_downs=subset(x = ATAC, downsample = n_cells)

all_ids=c(luma, lumb, basal, her2)

for (id in all_ids){

print(id)

sel_list=c(id,'Basal progenitors','Lum mature','Lum progenitors')

atac_s=subset(ATAC_downs, test %in% sel_list)


#Use slot counts:
x=atac_s@assays$peaksMACS2@counts

peaks_e=rowSums(x)

peaks_e=peaks_e[order(-peaks_e)]

#Select the most covered across cells peaks: 
peaks=names(peaks_e[1:50000])

x1=x[rownames(x) %in% peaks,]
atac_s[['for_m']]<-CreateAssayObject(x1)


#Load Seurat object
seurat_object <- atac_s

#Extract data, phenotype data, and feature data from the SeuratObject

data <- as(as.matrix(seurat_object@assays$for_m@counts), 'sparseMatrix')

pd <- new('AnnotatedDataFrame', data = seurat_object@meta.data)

fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))

fd <- new('AnnotatedDataFrame', data = fData)



#Construct monocle cds
monocle_cds <- newCellDataSet(data,
                         phenoData = pd,
                         featureData = fd,
                         lowerDetectionLimit = 0.5,
                         expressionFamily = negbinomial.size())


monocle_cds<-estimateSizeFactors(monocle_cds)

monocle_cds<-estimateDispersions(monocle_cds)

monocle_cds <- reduceDimension(monocle_cds, method = 'DDRTree', max_components =n_components)

monocle_cds <- orderCells(monocle_cds)


monocle_cds$ID=monocle_cds$cell_type_manual

subtype='NA'
subtype=ifelse(id %in% luma, 'LumA', subtype)
subtype=ifelse(id %in% lumb, 'LumB', subtype)
subtype=ifelse(id %in% basal, 'Basal', subtype)
subtype=ifelse(id %in% her2, 'Her2', subtype)

monocle_cds$ID=ifelse(monocle_cds$cell_type_manual=='Tumor',subtype,monocle_cds$cell_type_manual)
monocle_cds$ID=factor(monocle_cds$ID)

#Plot cells along the trajectory colored by cell group:
p1=plot_cell_trajectory(monocle_cds, color_by = "ID", cell_size = 1) +
scale_color_manual(values=cols)

#Save the plot:
pdf(paste('plots/',n_cells,'_MergedObj','_',id,'_',subtype,'_', n_components,
'_maxComp_subsetted_monocle.pdf',sep=''),width=6,height=5,useDingbats=F)
print(p1)
dev.off()

#Save RDS file:
saveRDS(monocle_cds,paste('Monocle_RDS/',n_cells,'_MergedObj','_',id,'_',subtype,'_', n_components,
'_monocle.rds',sep=''))
}