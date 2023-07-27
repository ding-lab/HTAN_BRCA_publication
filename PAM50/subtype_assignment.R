###
# input variables for the subtype prediction script
###
library(ctc)
library(heatmap.plus)

setwd("examplepath/PAM50_code") # Directory containing the "PAM50_code" scripts

paramDir<- "." # the location of unchanging files such as the function library and main program
inputDir<- "/examplepath/example_files"  # the location of the data matrix, and where output will be located

inputFile<- "example_input.tsv" # the input data matrix as a tab delimited text file
short<-"ExampleData" # short name that will be used for output files

datatype <- "sn" # must be sn or bulk
if(!(datatype %in% c("sn","bulk"))) {print("datatype must be sn or bulk")}
if(datatype == "sn") {calibrationParameters <- 24}
if(datatype == "bulk") {calibrationParameters <- 18}
#ER median sn = 24, bulk = 18 sn, -1 none	
#the column of the "mediansPerDataset.txt" file to use for calibration; 

# parameters we don't change
hasClinical<-FALSE
collapseMethod<-"mean"
####
# run the assignment algorithm
####

source(paste(paramDir,"subtypePrediction_functions_RWB.R",sep="/"))
source(paste(paramDir,"subtypePrediction_distributed.R",sep="/"))
