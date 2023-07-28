#post scnic

# filtering the regulon

library(SCopeLoomR)
library(Matrix)
library(AUCell)
library(SCENIC)
library(ComplexHeatmap)

# wdir = '/diskmnt/Projects/Users/ritalu/20220331_HTAN_BRCA/02_pySCENIC/object/'
# outdir = '/diskmnt/Projects/Users/ritalu/20220331_HTAN_BRCA/02_pySCENIC/20220418/object/'
setwd("")
vsnDir <- "."

loom <- open_loom('snRNA_pyscenic_output.loom')
# Read information from loom file:
regulons_incidMat <- get_regulons(loom,column.attr.name = "Regulons")

regulons <- regulonsToGeneLists(regulons_incidMat)
regulonsAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
regulonsAucThresholds <- get_regulon_thresholds(loom)
regulonsAucThresholds_2=names(regulonsAucThresholds)
names(regulonsAucThresholds_2)=regulonsAucThresholds
#embeddings <- get_embeddings(loom)

saveRDS(regulons, 'snRNA_Obj_regulons.20220504.rds')

###Follow the tutorial here http://htmlpreview.github.io/?https://github.com/aertslab/SCENIC/blob/master/inst/doc/SCENIC_Running.html

###Select the regulons:
res=getAUC(regulonsAUC)
#read cell annotation file
info=read.table('.CellInfo_subtype_tumor.0425.tsv',sep='\t',header=T)
rownames(info)=info$Barcode
info=info[colnames(res),]
# info$replicon_group=gsub('(.*)_(.*)_(.*)','\\1',info$Barcode)
# info$replicon_group=gsub('(.*)_(.*)','\\1',info$replicon_group)
# info$replicon_group=ifelse(info$Piece_ID %in% c("HT268B1-Th1H3", "HT029B1-S1PC", "HT035B1-S1PA","HT1408-06",
# "HT141B1-S1H1","HT206B1-S1H4", "HT271B1-S1H3"),'BRCA_Basal', info$replicon_group)

####First do this only for tumor cells:
#info_s=info[info$Cell_type=='Tumor',]
info_s = info
res_s=res[,info$Barcode]

cell_types=unique(info$replicon_group)

final_wilcoxon_stat=NULL
final_wilcoxon_stat_all=NULL
for (cell_t1 in cell_types){
    print(cell_t1)
    res_1=res_s[,colnames(res_s) %in% rownames(info_s)[info_s$replicon_group==cell_t1]]
    res_2=res_s[,colnames(res_s) %in% rownames(info_s)[info_s$replicon_group!=cell_t1]]
    all_wilcoxon_stat=NULL
    for (motif in 1:nrow(res_s)){
           mean_score1=mean(as.numeric(as.character(unlist(res_1[motif,]))),na.rm=TRUE)
           mean_score2=mean(as.numeric(as.character(unlist(res_2[motif,]))),na.rm=TRUE)
           w_test=wilcox.test(as.numeric(as.character(unlist(res_1[motif,]))),
as.numeric(as.character(unlist(res_2[motif,]))))
           stat=cbind(cell_t1,rownames(res_s)[motif],mean_score1,mean_score2,w_test$p.value)
           all_wilcoxon_stat=rbind(all_wilcoxon_stat,stat)
        }
        all_wilcoxon_stat=as.data.frame(all_wilcoxon_stat)
        all_wilcoxon_stat$V5=as.numeric(as.character(unlist(all_wilcoxon_stat$V5)))
        all_wilcoxon_stat=all_wilcoxon_stat[order(all_wilcoxon_stat$V5),]
        all_wilcoxon_stat$mean_score1=as.numeric(as.character(unlist(all_wilcoxon_stat$mean_score1)))
        all_wilcoxon_stat$mean_score2=as.numeric(as.character(unlist(all_wilcoxon_stat$mean_score2)))
        all_wilcoxon_stat$diff=all_wilcoxon_stat$mean_score1-all_wilcoxon_stat$mean_score2
        all_wilcoxon_stat=all_wilcoxon_stat[order(-all_wilcoxon_stat$diff),]

colnames(all_wilcoxon_stat)[1:2]=c('cell_t1','Regulon')
colnames(all_wilcoxon_stat)[5]=c('pvalue')
all_wilcoxon_stat$FDR=p.adjust(all_wilcoxon_stat$pvalue, method="fdr")
final_wilcoxon_stat=rbind(final_wilcoxon_stat,all_wilcoxon_stat)
}

final_wilcoxon_stat=as.data.frame(final_wilcoxon_stat)
#colnames(final_wilcoxon_stat)[1:2]=c('cell_t1','TF_Name')
final_wilcoxon_stat$cell_t2='Other'
write.table(final_wilcoxon_stat,'Regulons_AUC_difference.20220504.tsv',sep='\t',quote=F,row.names=F)

###Annotate regulons:
all_st=NULL
for (i in 1:length(regulons)){
    reg=names(regulons)[[i]]
    genes_n=length(regulons[[i]])
    st=cbind(reg,genes_n)
    all_st=rbind(all_st,st)
}
all_st=as.data.frame(all_st)
colnames(all_st)=c('Regulon','Genes_N')
write.table(all_st,'Regulon_annot.20220504.tsv',sep='\t',quote=F,row.names=F)




###try make bin act heatmaps using tutorial here: http://htmlpreview.github.io/?https://github.com/aertslab/SCENIC/blob/master/Tutorials_JupyterNotebooks/SCENIC_tutorial_2-ExploringOutput.html
set.seed(123)


#sample each cancer type by 300 cells
celltypes=unique(info_s$replicon_group)
nCells=300

all_cells_s=NULL
for (cell in celltypes){
    info_1=info_s[info_s$replicon_group==cell,]
    cells_s=sample(rownames(info_1),n_cells)
    all_cells_s=c(all_cells_s,cells_s)
}

#cellsSelected <- sample(colnames(regulonsAUC), nCells)
binarizeAUC <- function(auc, thresholds)
{
  thresholds <- thresholds[intersect(names(thresholds), rownames(auc))]
  regulonsCells <- setNames(lapply(names(thresholds),
                                   function(x) {
                                     trh <- thresholds[x]
                                     names(which(getAUC(auc)[x,]>trh))
                                   }),names(thresholds))

  regulonActivity <- reshape2::melt(regulonsCells)
  binaryRegulonActivity <- t(table(regulonActivity[,1], regulonActivity[,2]))
  class(binaryRegulonActivity) <- "matrix"

  return(binaryRegulonActivity)
}

binaryRegulonActivity <- binarizeAUC(regulonsAUC, regulonsAucThresholds_2)
binAct_subset <- binaryRegulonActivity[, which(colnames(binaryRegulonActivity) %in% all_cells_s)]
dim(binAct_subset)

write.table(binAct_subset,paste0('snRNA_',nCells,'_RandomCellsPercelltype.20220504.tsv'),sep='\t',quote=F)
write.table(binaryRegulonActivity,paste0('snRNA_',nCells,'_Percelltype.20220504.tsv'),sep='\t',quote=F)


############Try to identify Diff. active regulons based on binarized activity:
# binaryRegulonActivity

#info_s=info[info$Cell_type=='Tumor',]
res_s=binaryRegulonActivity[,info_s$Barcode]

cell_types=unique(info_s$replicon_group)

final_wilcoxon_stat=NULL
final_wilcoxon_stat_all=NULL
for (cell_t1 in cell_types){
    print(cell_t1)
    res_1=res_s[,colnames(res_s) %in% rownames(info_s)[info_s$replicon_group==cell_t1]]
    res_2=res_s[,colnames(res_s) %in% rownames(info_s)[info_s$replicon_group!=cell_t1]]
    all_wilcoxon_stat=NULL
    for (motif in 1:nrow(res_s)){
           mean_score1=mean(as.numeric(as.character(unlist(res_1[motif,]))),na.rm=TRUE)
           mean_score2=mean(as.numeric(as.character(unlist(res_2[motif,]))),na.rm=TRUE)
           w_test=wilcox.test(as.numeric(as.character(unlist(res_1[motif,]))),
as.numeric(as.character(unlist(res_2[motif,]))))
           stat=cbind(cell_t1,rownames(res_s)[motif],mean_score1,mean_score2,w_test$p.value)
           all_wilcoxon_stat=rbind(all_wilcoxon_stat,stat)
        }
        all_wilcoxon_stat=as.data.frame(all_wilcoxon_stat)
        all_wilcoxon_stat$V5=as.numeric(as.character(unlist(all_wilcoxon_stat$V5)))
        all_wilcoxon_stat=all_wilcoxon_stat[order(all_wilcoxon_stat$V5),]
        all_wilcoxon_stat$mean_score1=as.numeric(as.character(unlist(all_wilcoxon_stat$mean_score1)))
        all_wilcoxon_stat$mean_score2=as.numeric(as.character(unlist(all_wilcoxon_stat$mean_score2)))
        all_wilcoxon_stat$diff=all_wilcoxon_stat$mean_score1-all_wilcoxon_stat$mean_score2
        all_wilcoxon_stat=all_wilcoxon_stat[order(-all_wilcoxon_stat$diff),]

colnames(all_wilcoxon_stat)[1:2]=c('cell_t1','Regulon')
colnames(all_wilcoxon_stat)[5]=c('pvalue')
all_wilcoxon_stat$FDR=p.adjust(all_wilcoxon_stat$pvalue, method="fdr")
final_wilcoxon_stat=rbind(final_wilcoxon_stat,all_wilcoxon_stat)
}

final_wilcoxon_stat=as.data.frame(final_wilcoxon_stat)
#colnames(final_wilcoxon_stat)[1:2]=c('cell_t1','TF_Name')
final_wilcoxon_stat$cell_t2='Other'
write.table(final_wilcoxon_stat, paste0('Regulons_BinAct_difference_',nCells,'.20220504.tsv'),sep='\t',quote=F,row.names=F)


## re-filtering N_genes per regulon
## N_genes > 50
## replot 0510

library(ComplexHeatmap)
library(tidyverse)
col_cell_type = c("Basal_tumor" = "#E41A1C",
                  "Luminal_tumor" = "#984EA3",
                  "Basal_progenitor" = "#377EB8",
                  "Luminal_progenitor" = "#4DAF4A",
                  "Luminal_mature" = "#FF7F00")
col_sample = c('HT027B1' = '#6B5DDD' ,'HT029B1' = '#B2E977' ,'HT035B1' = '#824D8F' ,'HT036B1' = '#DCE6AE' ,'HT065B1' = '#6D9DDD' ,
               'HT088B1' = '#E77AA6' ,'HT110B1' = '#A434E8' ,'HT128B1' = '#67EBEC' ,'HT137B1' = '#C1E9D2' ,'HT141B1' = '#99E3AD',
               'HT163B1' = '#5A7A87' ,'HT206B1' = '#E0E644' ,'HT217B1' = '#9AA091' ,'HT235B1' = '#DBB67C' ,
               'HT243B1' = '#D7DD81' ,'HT262B1' = '#DA79E0' ,'HT263B1' = '#ECE1D3' ,'HT271B1' = '#BADCE6' ,
               'HT297B1' = '#CAC1DD' ,'HT305B1' = '#DD3FC3','HT308B1' = '#E5A4DC' ,'HT323B1' = '#6DEC95' ,
               'HT339B1' = '#4BD4B0' , 'HT384B1' = '#6AB8AE' )


RNA_binAct_subset3 = read.table('/Users/ritalu/Box/LD_Lab/20220331_HTAN_BRCA/02_pyscenic/20220504/Fig/celltypes_binAct_filter_300_0504.tsv',sep='\t')
colnames(RNA_binAct_subset3)<- sub(".", "-", colnames(RNA_binAct_subset3), fixed = TRUE)

regulon_N = read.table('/Users/ritalu/Box/LD_Lab/20220331_HTAN_BRCA/02_pyscenic/20220504/object/Regulon_annot.20220504.tsv',header = T,sep='\t')

#select regulon has 50 genes
nGenes = 50
RNA_select = RNA_binAct_subset3[rownames(RNA_binAct_subset3) %in% regulon_N[regulon_N$Genes_N >= nGenes,]$Regulon, ]

cellID = read.table('/Users/ritalu/Box/LD_Lab/20220331_HTAN_BRCA/02_pyscenic/20220418/object/CellInfo_subtype_tumor.0425.tsv',sep='\t',header=T)

cellID2 = cellID[cellID$Barcode %in% colnames(RNA_select),]

top = HeatmapAnnotation(cell_type= cellID2$replicon_group,
                        #Cluster=Annot5$Cluster,
                        Sample = cellID2$Sample,
                        
                        col=list(cell_type = col_cell_type,
                                 Sample = col_sample), 
                        
                        na_col = "white",
                        annotation_name_gp= gpar(fontsize = 12),
                        annotation_height =0.8) 
nCells = 300
rnaplot_super = draw(ComplexHeatmap::Heatmap(as.matrix(RNA_select[, cellID2$Barcode]) , name="Binarized activity", 
                                          col = c("white", "black"),top_annotation = top,
                                          cluster_rows = TRUE, cluster_columns=TRUE,
                                          clustering_method_rows = "ward.D2",
                                          clustering_method_columns = "ward.D2",
                                          show_column_names = FALSE, column_split =  cellID2$replicon_group,
                                          column_gap =  unit(1, "mm"),column_title = NULL,
                                          row_names_gp=grid::gpar(fontsize= 12)))

outfig = '/Users/ritalu/Box/LD_Lab/20220331_HTAN_BRCA/02_pyscenic/20220510_Final/'
pdf(paste0(outfig, 'Binarized_Activity_filtered_ngenes', nGenes,'_cells_',nCells,'_supervised_combo_0510.pdf'),width=12,height=8,useDingbats=FALSE)
draw(rnaplot_super ,adjust_annotation_extension = TRUE,annotation_legend_side = "left",heatmap_legend_side = "left")
dev.off()

