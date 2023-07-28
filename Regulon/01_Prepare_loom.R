#remake loom file

library(SCENIC)
library(SCopeLoomR)


wdir = '/diskmnt/Projects/Users/ritalu/20220331_HTAN_BRCA/02_pySCENIC/object/'
outdir = '/diskmnt/Projects/Users/ritalu/20220331_HTAN_BRCA/02_pySCENIC/20220418/object/'


exprMatrix <- read.table(paste0(wdir,'exprMatrix_snRNA_0418.tsv'),header=T)
colnames(exprMatrix)<- sub(".", "-", colnames(exprMatrix), fixed = TRUE)
cellInfo <- read.table(paste0(wdir,'cellInfo_snRNA_0418.tsv'), header=T)
row.names(cellInfo) <-cellInfo$Barcode
cellInfo_1 <- cellInfo[,c['Piece_ID', 'cell_type_specific','Barcode' )

loom <- build_loom(paste0(outdir , "snRNA_obj.loom"), dgem=exprMatrix)
loom <- add_cell_annotation(loom, cellInfo_1)
close_loom(loom)

#check
#loomPath <- system.file(package="SCENIC", "examples/mouseBrain_toy.loom"
loom2 <- open_loom('/diskmnt/Projects/Users/ritalu/20220331_HTAN_BRCA/02_pySCENIC/20220418/object/snRNA_obj_v2.loom')
exprMat <- get_dgem(loom2)
cellInfo2 <- get_cell_annotation(loom2)


#go to filter then pyscenic
