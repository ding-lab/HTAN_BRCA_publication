library(Signac)
library(Seurat)
library(reshape)
library(reshape2)
library(EnsDb.Hsapiens.v86)
library(tidyverse)

path_to_ATAC_merged_obj=''
Piece_ID_subtype_annotation=''

atac=readRDS(path_to_ATAC_merged_obj)

#Add assignment of normal duct cells for these cluters:
#13 - Lum mature
#15 - Basal prog
#30 - Lum prog

#Use cell_type_manual for cell type annotation
atac$cell_type_manual=ifelse(atac$seurat_clusters==13,"Lum mature",atac$cell_type_manual)
atac$cell_type_manual=ifelse(atac$seurat_clusters==15,"Basal progenitors",atac$cell_type_manual)
atac$cell_type_manual=ifelse(atac$seurat_clusters==30,"Lum progenitors",atac$cell_type_manual)


#Keep tumor and normal duct cells, and also remove cells annotated as tumor from NATs:
atac_all=atac
atac=subset(atac_all, cell_type_manual %in% c('Lum mature','Basal progenitors','Lum progenitors') |
(cell_type_manual=='Tumor' & Piece_ID!='HT384B1_M1'))


annot=read.table(Piece_ID_subtype_annotation, sep='\t', header=T)

luma=unique(annot$Sample_ID[annot$PAM50_sn %in% c('LumA')])
lumb=unique(annot$Sample_ID[annot$PAM50_sn %in% c('LumB')])
basal=unique(annot$Sample_ID[annot$PAM50_sn %in% c('Basal')])
her2=unique(annot$Sample_ID[annot$PAM50_sn %in% c('Her2')])


atac$test=atac$cell_type_manual
atac$test=ifelse(atac$Piece_ID %in% luma & atac$cell_type_manual=='Tumor','LumA_tumor',atac$test) 
atac$test=ifelse(atac$Piece_ID %in% lumb & atac$cell_type_manual=='Tumor','LumB_tumor',atac$test) 
atac$test=ifelse(atac$Piece_ID %in% basal & atac$cell_type_manual=='Tumor','Basal_tumor',atac$test) 
atac$test=ifelse(atac$Piece_ID %in% her2 & atac$cell_type_manual=='Tumor','Her2_tumor',atac$test) 


Idents(atac)=atac$test
Idents(atac)=factor(Idents(atac),levels=c("Basal progenitors", "Lum mature", "Lum progenitors", "LumA_tumor", 
"LumB_tumor", "Her2_tumor", "Basal_tumor"))

#Add annotation to the object:
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
genome(annotations) <- "NA"

#Change annotation to UCSC style
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

#Add the gene information to the object
Annotation(atac) <- annotations
annot=Annotation(atac)

#Plot coverage for promoters of all selected PAM50-genes (genes highlighted in Fig. 2C).

her2_genes=c('TMEM45B','CDC6','GRB7','ERBB2')
basal_genes=c('SFRP1','SFRP1','TYMS','MELK','MKI67','FOXC1','MYBL2','KIF2C','ANLN','CCNB1','UBE2T','EXO1',
'BIRC5','CENPF','PHGDH','UBE2C','MIA','CCNE1','CEP55','PTTG1','CDC20','RRM2','ACTR3B','EGFR','CDH3','KRT14',
'KRT5','KRT17','BAG1','MYC')
lum_genes=c('PGR','NAT1','ESR1','GPR160','FOXA1','FGFR4','MLPH','BCL2','MDM2','CXXC5','SLC39A6','MAPT','MMP11')

#These were selected for Fig. 2E:
#her2_genes=c('ERBB2','GRB7')
#basal_genes=c('SFRP1','KRT17')
#lum_genes=c('ESR1')

genes_selected=c(her2_genes,basal_genes,lum_genes)

for (gene in genes_selected){

if (gene %in% her2_genes){
   subt='Her2'
}else if (gene %in% basal_genes){
   subt='Basal'
}else{
   subt='Lum'
}
#Extract gene coordinates
region=as.data.frame(LookupGeneCoords(atac, gene, assay = NULL))

if(nrow(region)>0){
strand=as.character(unique(strand(annot[annot$gene_name==gene,])))[1]
chr=region$seqnames
if (strand=='+'){
st=region$start
}else if (strand=='-'){
st=region$end
}
en=region$end

#Add 2K bp before and after TSS
new_st=as.numeric(st)-2000
new_en=as.numeric(st)+2000
new_peak=paste(chr,new_st,new_en,sep='-')
print(paste(gene,' ',strand,sep=''))

cov_plot=CoveragePlot(
  object = atac,
  region = new_peak,
  annotation = FALSE,
  peaks = FALSE,
  links=FALSE
)
gene_plot <- AnnotationPlot(
  object = atac,
  region = new_peak
)
p=CombineTracks(
  plotlist = list(cov_plot,gene_plot)
)
pdf(paste("plots/",subt,"/Coverage_",new_peak,"_",gene,"_",subt,".pdf",sep=""),width=4.5,height=5.5,
useDingbats=FALSE)
print(p)
dev.off()
}
}
