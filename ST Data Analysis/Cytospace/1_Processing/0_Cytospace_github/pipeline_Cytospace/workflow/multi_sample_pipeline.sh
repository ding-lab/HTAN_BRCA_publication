eval "$(conda shell.bash hook)"

# Variables: ------------------------------------------
# Change and assign following values per instructions!
# A. Project folder and sample info table----------------------------------
PROJECT_FOLDER="/PATH/TO/PROJECT/FOLDER/"
# Sample table should have 4 columns. Column name not needed in the file. See example file for details
# columns: snRNA_Name ST_Name ST_SeuratObject file_exists
SAMPLE_TABLE_PATH="$PROJECT_FOLDER/table/sample_table.csv"

## Part B: Preprocess snRNA-seq ----------------------------------
## B1. Parameter for split the MERGED single cell snRNA object
# NOTE this is optional! If already have the individual matching snRNA object. follow section B2
# RUN WITH SPLTTING MERGED SNRNA OBJECT - Comment this section out if following B2 section!
SPLIT_SN="True" # True or False
# Path to the MERGED single cell snRNA object 
sn_MERGED_OBJ_PATH="PATH/TO/integrated/MERGED_snRNA_OBJECT.qs"
# Colunm name for sample name in the MERGED single cell snRNA objec
SAMPLE_COLUMN="" # e.g. "Sample"
# Make root folder to output split snRNA objects
SN_ROOT="${PROJECT_FOLDER}/output/1_preprocess/1_Split_snRNA/" 
mkdir -p $SN_ROOT
## ----------------------------------------------------------------------

## B2. snRNA object pre-split already. Use this directly ----------------------------------
# If Already split the file, set SPLIT_SN=False, and make sure the SN_ROOT path point to the right folder
SPLIT_SN="False" # True or False
SN_ROOT="PATH/TO/SPLIT/SNRNA_OBJTECT/"
# The file should have structure: "${SN_ROOT}/individual/${SN_NAME}.qs"
## ----------------------------------------------------------------------

# C. Parameters for cell type transfer ----------------------------------
CELLTYPE_COLUMN="" # e.g. "cell_type_specific"


# Workflow starts here----------------------------------
# Part1 Read sample table
# Read the content of each column and store them as arrays
# Read as tsv or csv
#read -r header # read and discard the first line
while IFS= read -r line; do
    if [[ "$line" == *$'\t'* ]]; then
        IFS=$'\t' read -r snRNA_Name ST_Name ST_SeuratObject file_exists <<< "$line"
    else
        IFS=',' read -r snRNA_Name ST_Name ST_SeuratObject file_exists <<< "$line"
    fi
    sn_SAMPLE_NAMES+=("$snRNA_Name")
    ST_SAMPLE_NAMES+=("$ST_Name")
    ST_FILE_PATHS+=("$ST_SeuratObject")
    file_exists_array+=("$file_exists")
done < "$SAMPLE_TABLE_PATH"


# Print the arrays
echo "sn_SAMPLE_NAMES: ${sn_SAMPLE_NAMES[@]}"
echo "ST_SAMPLE_NAMES: ${ST_SAMPLE_NAMES[@]}"
echo "ST_FILE_PATHS: ${ST_FILE_PATHS[@]}"
echo "file_exists_array: ${file_exists_array[@]}"



# Constants ----------------------------------
# pipeline path. No need to change 
PREPROCESS_SCRIPT_ROOT="/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/Pipeline/RCTD/script/"


# Workflow Scripts ----------------------------------
### ----------- 1 Split snRNA object ----------- ###
# check if SPLIT_SN is true or false
if [ "$SPLIT_SN" = "True" ]; then
    echo "Splitting snRNA object"
    conda activate Signac
    echo "Running 1. Split snRNA object"
    RUN_SCRIPT="$PREPROCESS_SCRIPT_ROOT/1_Split_snRNA_object.r"
    #Rscript $RUN_SCRIPT --input $sn_MERGED_OBJ_PATH --output $SN_ROOT --sample_column $SAMPLE_COLUMN --num_cpu 50
else
    echo "SPLIT_SN is False. Not splitting snRNA object"
    continue
fi


### ----------- 2 Call Cytospace single sample pipeline ----------- ###
conda activate cytospace
echo "Running 2. Call Cytospace single sample pipeline"


# global variables for Cytospace
script_root="/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/Pipeline/Cytospace/script/"
workflow_root="/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/Pipeline/Cytospace/workflow/"
SINGLE_SAMPLE_SCRIPT="$workflow_root/single_sample_pipeline.sh"

# Loop throught snRNA + ST paires
for idx in "${!ST_SAMPLE_NAMES[@]}"
do
    # Assign variables
    ST_LIBRARYNAME=${ST_SAMPLE_NAMES[$idx]}
    ST_FILE=${ST_FILE_PATHS[$idx]}
    SN_NAME=${sn_SAMPLE_NAMES[$idx]}

    # Assign variables =================================
    ## ===== sample_variable ===== ##
    # General
    project_folder="${PROJECT_FOLDER}/Cytospace/${ST_LIBRARYNAME}"
    cell_type_column_name="${CELLTYPE_COLUMN}"

    # snRNA
    sn_sample_name="${SN_NAME}"
    sn_file_path="${SN_ROOT}/individual/${SN_NAME}.qs"

    # ST
    ST_sample_name="${ST_LIBRARYNAME}"
    ST_file_path="${ST_FILE}"
    # ========================================================
    
    echo "Processing ${SN_NAME} and ${ST_LIBRARYNAME}"
    sn_FILE=$SN_ROOT/individual/$SN_NAME.qs
    # first check if both sn_File and ST_File exists
    if [ ! -f "$sn_FILE" ]; then
        echo "$sn_FILE does not exist, skipping..."
        continue
    fi
    if [ ! -f "$ST_FILE" ]; then
        echo "$ST_FILE does not exist, skipping..."
        continue
    fi
    
    # run single sample script
    . "${SINGLE_SAMPLE_SCRIPT}"
done

