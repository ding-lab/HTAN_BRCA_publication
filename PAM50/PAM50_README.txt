PAM50 subtyping README
Author: Michael Iglesia

#-------------------------------------#
#     Background: PAM50 subtyping     #
#-------------------------------------#

The PAM50 subtyping algorithm assigns breast cancer samples to one of five subtypes, using 50 genes to make assignments employing a nearest centroid approach.
The five subtypes are Luminal A, Luminal B, HER2-enriched, Basal-like, and Normal-like.
Colors for plotting: #00008B, #ADD8E6, #FF69B4, FF0000, #008000.
This initial list of subtypes was originally published in Perou, et al. Nature 2000
Other useful papers for PAM50 subtyping background:
    Prat, et al. Molecular Oncology 2011
    Prat, et al. Breast Cancer Research & Treatment 2012
    The Cancer Genome Atlas Network, Nature 2012
    Prat, et al. The Breast 2015

In addition to the PAM50 subtyping method detailed below, the R package genefu can be used for subtyping from RNA-seq.


#-------------------------------------#
#        Formatting your data         #
#-------------------------------------#

Input data should be in a tab-delimited file with:

    first row = sample ID
    first column = gene symbol

For snRNA-seq data, each column should be an average gene value for a sample or Seurat cluster. To get these average values from a Seurat object, consider this approach:

Idents(obj) = obj$seurat_clusters
dataAll = AverageExpression(obj,assays="SCT",features=c("MIA","ACTR3B","ANLN","BAG1","BCL2","BIRC5","BLVRA","CCNB1","CCNE1","CDC20","CDC6","NUF2","CDH3","CENPF","CEP55","CXXC5","EGFR","ERBB2","ESR1","EXO1","FGFR4","FOXA1","FOXC1","GPR160","GRB7","KIF2C","NDC80","KRT14","KRT17","KRT5","MAPT","MDM2","MELK","MKI67","MLPH","MMP11","MYBL2","MYC","NAT1","ORC6","PGR","PHGDH","PTTG1","RRM2","SFRP1","SLC39A6","TMEM45B","TYMS","UBE2C","UBE2T"))

An example input file is provided: example_input.tsv

#-------------------------------------#
#          Running the code           #
#-------------------------------------#

The subtyping algorithm is called by running subtype_assignment.R after first editing parameters within this file.

Parameters to edit before running are:

    working dir:    This should point to the "PAM50_code" directory
    inputDir:       The directory where your input matrix of gene expression values is located
    short:          A short name which will be appended to all output files. Including a date here is wise.
    datatype:       Data type. Must be "sn" or "bulk"

Note that to run this code, you will need the ctc and heatmap.plus R packages installed.

#-------------------------------------#
#      Output and interpretation      #
#-------------------------------------#
Output

    [samplename]_pam50scores.txt Scores for each subtype, confidence value of assignment, and final assignment for all samples
    [samplename]__PAM50_normalized_heatmap.pdf Heatmap of PAM50 genes with color bars for subtype assignments and confidence of calls
    predictionScores_pam50RankCorrelation_1_[samplename].pdf For each subtype, shows box plots of correlation to each subtype centroid
    predictionScores_pam50RankCorrelation_2_[samplename].pdf PCA plot showing control and experimental assignments

Output files used for viewing interactive PAM50 heatmap in Java Treeview or similar application:

    [samplename]_PAM50_normalized_heatmap.atr
    [samplename]__PAM50_normalized_heatmap.cdt
    [samplename]__PAM50_normalized_heatmap.gtr

An example of script output can be found here:

/diskmnt/Projects/Users/miglesia/PAM50/example_files/

Caveats
    Clinical subtypes (ER, PR, and HER2-based) are a good sanity check and should be available in the clinical data. Most luminal A/B tumors should be ER+, most HER2-enriched tumors should be HER2+ (IHC 3+ or IHC 2+ and amplified), and most basal-like tumors should be triple-negative. A good example of expected concordance is here: Prat, et al. The Breast 2015