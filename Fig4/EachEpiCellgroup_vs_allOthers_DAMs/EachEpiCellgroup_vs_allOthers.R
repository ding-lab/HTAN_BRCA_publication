library(Signac)
library(Seurat)
library(plyr)
library(dplyr)


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

#Subset object to tumor and normal duct cells:
ATAC=subset(atac, cell_type_manual %in% c('Tumor','Lum mature','Basal progenitors','Lum progenitors'))

ATAC$test=ifelse(ATAC$cell_type_manual=='Tumor',paste(ATAC$Piece_ID, ATAC$cell_type_manual,sep='_'),
ATAC$cell_type_manual)

annot=read.table(Piece_ID_subtype_annotation,sep='\t',header=T)
lum=unique(annot$Sample_ID[annot$PAM50_sn %in% c('LumA','LumB')])
basal=unique(annot$Sample_ID[annot$PAM50_sn %in% c('Basal')])
her2=unique(annot$Sample_ID[annot$PAM50_sn %in% c('Her2')])

lum=paste(lum, 'Tumor',sep='_')
basal=paste(basal, 'Tumor',sep='_')

ATAC$test=ifelse(ATAC$test %in% lum, "Lum", ATAC$test)
ATAC$test=ifelse(ATAC$test %in% basal, "Basal", ATAC$test)


DefaultAssay(ATAC) <- 'chromvar'

chromv= GetAssayData(object = ATAC)
jaspar=read.table('/diskmnt/Projects/HTAN_analysis/snATAC/Signac/CCRCC/JASPAR2020_motifs.txt',sep='\t',
header=TRUE)

mtx0=chromv
res=merge(mtx0,jaspar,by=0,all.x=TRUE)
rownames(res)=res$motif.name
res=res[,-1]
res=res[,1:(ncol(res)-2)]


ann_col0 = data.frame('cell_types' = ATAC$test)
ann_col0$cel_2=ann_col0$cell_types
ann_col0=ann_col0[order(ann_col0$cell_types),]

res=res[,rownames(ann_col0)]

cell_types=c("Basal", "Lum mature", "Lum progenitors", "Lum", "Basal progenitors")


#Run DAMs for one cell group vs all others
final_wilcoxon_stat=NULL

for (cell_t1 in cell_types){
    print (cell_t1)

    res_1=res[,colnames(res) %in% rownames(ann_col0)[ann_col0$cell_types==cell_t1]]
    res_2=res[,colnames(res) %in% rownames(ann_col0)[ann_col0$cell_types!=cell_t1]]
    all_wilcoxon_stat=NULL
        for (motif in 1:nrow(res)){
           mean_score1=mean(as.numeric(as.character(unlist(res_1[motif,]))),na.rm=TRUE)
           mean_score2=mean(as.numeric(as.character(unlist(res_2[motif,]))),na.rm=TRUE)
           w_test=wilcox.test(as.numeric(as.character(unlist(res_1[motif,]))),
as.numeric(as.character(unlist(res_2[motif,]))))
           stat=cbind(cell_t1,rownames(res)[motif],mean_score1,mean_score2,w_test$p.value)
           all_wilcoxon_stat=rbind(all_wilcoxon_stat,stat)
        }
        all_wilcoxon_stat=as.data.frame(all_wilcoxon_stat)
        all_wilcoxon_stat$V5=as.numeric(as.character(unlist(all_wilcoxon_stat$V5)))
        all_wilcoxon_stat=all_wilcoxon_stat[order(all_wilcoxon_stat$V5),]
        all_wilcoxon_stat$mean_score1=as.numeric(as.character(unlist(all_wilcoxon_stat$mean_score1)))
        all_wilcoxon_stat$mean_score2=as.numeric(as.character(unlist(all_wilcoxon_stat$mean_score2)))
        all_wilcoxon_stat$diff=all_wilcoxon_stat$mean_score1-all_wilcoxon_stat$mean_score2
        all_wilcoxon_stat=all_wilcoxon_stat[order(-all_wilcoxon_stat$diff),]

colnames(all_wilcoxon_stat)[1:2]=c('cell_t1','TF')
colnames(all_wilcoxon_stat)[5]=c('pvalue')
all_wilcoxon_stat$FDR=p.adjust(all_wilcoxon_stat$pvalue, method="fdr")

final_wilcoxon_stat=rbind(final_wilcoxon_stat,all_wilcoxon_stat)
}


final_wilcoxon_stat=as.data.frame(final_wilcoxon_stat)
colnames(final_wilcoxon_stat)[1:2]=c('cell_t1','TF_Name')

write.table(final_wilcoxon_stat,paste("out/Score_difference_EachCellGroup_vs_Others.tsv",
sep=""),quote=FALSE,sep="\t",row.names=FALSE)
