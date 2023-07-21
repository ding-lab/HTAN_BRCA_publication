eval "$(conda shell.bash hook)"

# Variables: ------------------------------------------
# Change and assign following values per instructions!
# A. Project folder and sample info table----------------------------------
PROJECT_FOLDER="/{PROJECT_PATH}/3_scRNA_by_subtype/4_runCytospace/"
# Sample table should have 4 columns. Column name not needed in the file. See example file for details
# columns: snRNA_Name ST_Name ST_SeuratObject file_exists
SAMPLE_TABLE_PATH="/{PROJECT_PATH}/1_RunCytospace/2_select_ST_with_subtype/out/ST_with_subtype.tsv"


# ## B3. Map use integrated snRNA object ----------------------------------
SPLIT_SN="False" # True or False
#sn_MERGED_OBJ_PATH="/diskmnt/Projects/HTAN_analysis_2/BRCA/scRNA/all_merged/scRNA_merged.rds"
SN_ROOT="/{PROJECT_PATH}/3_scRNA_by_subtype/1_split_sc_by_subtype/out/20230705/"
# Add new variable to check if using shared processed snRNA matrix
USE_SHARED_SN="True" # True or False

# C. Parameters for cell type transfer ----------------------------------
CELLTYPE_COLUMN="old_cell_type_specific"


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
PREPROCESS_SCRIPT_ROOT="./HTAN_BRCA_publication/Individual_scripts/Simon/1_Processing/0_preprocess"


# Workflow Scripts ----------------------------------
### ----------- 1 Split snRNA object ----------- ###
# check if SPLIT_SN is true or false
if [ "$SPLIT_SN" = "True" ]; then
    echo "Splitting snRNA object"
    conda activate Signac
    echo "Running 1. Split snRNA object"
    RUN_SCRIPT="$PREPROCESS_SCRIPT_ROOT/1_Split_snRNA_object.r"
    RUN_SCT="True"
else
    echo "SPLIT_SN is False. Not splitting snRNA object"
    RUN_SCT="False"
fi


### ----------- 2 Call Cytospace single sample pipeline ----------- ###
conda activate cytospace
echo "Running 2. Call Cytospace single sample pipeline"


# global variables for Cytospace
# Use modified Cytospace pipeline to allow using shared snRNA matrix
SCRIPT_ROOT="/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/BRCA_revision/2_github_organized/HTAN_BRCA_publication/Individual_scripts/Simon/1_Analysis/0_Cytospace_github/pipeline_Cytospace/"
script_root="${SCRIPT_ROOT}/script/"
workflow_root="${SCRIPT_ROOT}/workflow/"
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
    # Individual object mode
    sn_file_path="${SN_ROOT}/${SN_NAME}"
    # Merged object mode
    #sn_file_path=$sn_MERGED_OBJ_PATH

    # ST
    ST_sample_name="${ST_LIBRARYNAME}"
    ST_file_path="${ST_FILE}"
    # ========================================================
    
    echo "Processing ${SN_NAME} and ${ST_LIBRARYNAME}"
    sn_FILE=$sn_file_path
    
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

