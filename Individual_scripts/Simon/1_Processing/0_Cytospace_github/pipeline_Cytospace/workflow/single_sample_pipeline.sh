eval "$(conda shell.bash hook)"
script_root="/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/Pipeline/Cytospace/script/"

## ===== sample_variable ===== ##
# For single samples, we can hard code the variables here
# # General
# project_folder="PROJECT/FOLDER"
# cell_type_column_name="CELL_TYPE_COLUMN_NAME"

# # snRNA
# sn_sample_name="SN_SAMPLE_NAME"
# sn_file_path="SN_FILE_PATH"

# # ST
# ST_sample_name="ST_SAMPLE_NAME"
# ST_file_path="ST_FILE_PATH"


## ===== folder structures ===== ##
out_path="$project_folder/output/"
step_names=(
    "1_snRNA_input"
    "2_ST_input"
    "3_cytospace_output"
    "4_snRNA_with_SCT"
    "5_create_cytospace_seurat"
    "6_plot_cytospace_distribution"
)

for step_name in "${step_names[@]}"
do
    step_path="${out_path}/${step_name}"
    mkdir -p "$step_path"
done

## ===== workflow ===== ##
# STEP1: Prep snRNA
# Options:
#   -n SN_SAMPLE_NAME, --sn_sample_name=SN_SAMPLE_NAME
# 		sn sample name

# 	-f SN_FILE_PATH, --sn_file_path=SN_FILE_PATH
# 		sn file path

# 	-c CELL_TYPE_COLUMN_NAME, --cell_type_column_name=CELL_TYPE_COLUMN_NAME
# 		cell type column name

# 	-o OUT_PATH, --out_path=OUT_PATH
# 		out path

# 	-s SCRIPT_ROOT, --script_root=SCRIPT_ROOT
# 		script root

# Step1: Prep snRNA
# Check if output file scRNA_data.txt exists
if test -f "${out_path}/${step_names[0]}/scRNA_data.txt"; then
    echo "Step1 output file exists. Skip Step1"
else
    echo "Step1: Prep snRNA"
    conda activate seurat4.3

    Rscript "${script_root}/1_prep_snRNA.r" \
        -n "${sn_sample_name}" \
        -f "${sn_file_path}" \
        -c "${cell_type_column_name}" \
        -o "${out_path}/${step_names[0]}" \
        -s "${script_root}"
fi

# Step2: Prep ST ---------------------
# Options:
# 	-n SAMPLE_NAME, --sample_name=SAMPLE_NAME
# 		sample name

# 	-f ST_FILE_PATH, --ST_file_path=ST_FILE_PATH
# 		ST file path

# 	-o ST_OUT_DIR, --st_out_dir=ST_OUT_DIR
# 		ST out dir

# 	-s SCRIPT_ROOT, --script_root=SCRIPT_ROOT
# 		script root

# STEP2: Prep ST 
if test -f "${out_path}/${step_names[1]}/ST_data.txt"; then
    echo "Step2 output file exists. Skip Step2"
else
    echo "Step2: Prep ST"
    conda activate seurat4.3
    Rscript "${script_root}/2_prep_ST.r" \
        -n "${ST_sample_name}" \
        -f "${ST_file_path}" \
        -o "${out_path}/${step_names[1]}" \
        -s "${script_root}"
fi


# Step3: run cytospace ---------------------
echo "Step3: run cytospace"
conda activate cytospace
cytospace --scRNA-path "${out_path}/${step_names[0]}/scRNA_data.txt" \
    --cell-type-path "${out_path}/${step_names[0]}/cell_type_labels.txt" \
    --st-path "${out_path}/${step_names[1]}/ST_data.txt" \
    --coordinates-path "${out_path}/${step_names[1]}/Coordinates.txt" \
    -o "${out_path}/${step_names[2]}" \
    --solver-method lap_CSPR

# Step4: run SCT on snRNA/scRNA object ---------------------
# Note this step is to make sure original snRNA/scRNA object has the SCT assay calculated
# If existed, will just save in the new folder. If not will SCT and UMAP.

# Options:
# 	-s SN_FILE_PATH, --sn_file_path=SN_FILE_PATH
# 		snRNA object

# 	-n SN_SAMPLE_NAME, --sn_sample_name=SN_SAMPLE_NAME
# 		snRNA name

# 	-o OUT_PATH, --out_path=OUT_PATH
# 		out dir

echo "Step4: run SCT on snRNA/scRNA object"
conda activate seurat4.3
Rscript "${script_root}/4_run_SCT_snRNA.r" \
    -s "${sn_file_path}" \
    -n "${sn_sample_name}" \
    -o "${out_path}/${step_names[3]}"


# Step5: Create Cytospace Seurat Object ---------------------
# Options:
# 	-c CYTOSPACE_PATH, --cytospace_path=CYTOSPACE_PATH
# 		cytospace result path

# 	-s SNRNA_OBJ, --snRNA_obj=SNRNA_OBJ
# 		snRNA object

# 	-S ST_OBJ, --ST_obj=ST_OBJ
# 		ST object

# 	-o OUT_PATH, --out_path=OUT_PATH
# 		out dir

# 	-n ST_SAMPLE_NAME, --ST_sample_name=ST_SAMPLE_NAME
# 		ST sample name

# 	-r SCRIPT_ROOT, --script_root=SCRIPT_ROOT
# 		script root
echo "Step5: Create Cytospace Seurat Object"
conda activate seurat4.3
Rscript "${script_root}/5_create_cytospace_seurat_obj.r" \
    -c "${out_path}/${step_names[2]}/" \
    -s "${out_path}/${step_names[3]}/${sn_sample_name}.qs" \
    -S "${ST_file_path}" \
    -o "${out_path}/${step_names[4]}" \
    -n "${ST_sample_name}" \
    -r "${script_root}"

# Step6: Plot Cytospace distribution ---------------------
# Options:
# 	-s SAMPLE_NAME, --sample_name=SAMPLE_NAME
# 		sample name

# 	-c CYTOSPACE_SEURAT_PATH, --cytospace_seurat_path=CYTOSPACE_SEURAT_PATH
# 		cytospace seurat path

# 	-o OUT_PATH, --out_path=OUT_PATH
# 		out dir

# 	-l CELL_TYPE_COLUMN_NAME, --cell_type_column_name=CELL_TYPE_COLUMN_NAME
# 		cell type column name

# 	-p PT_SIZE_FACTOR, --pt_size_factor=PT_SIZE_FACTOR
# 		point size factor

# 	-a IMAGE_BG_ALPHA, --image_bg_alpha=IMAGE_BG_ALPHA
# 		image background alpha

# 	-t LABEL_ST, --label_st=LABEL_ST
# 		label ST

#   -r SCRIPT_ROOT, --script_root=SCRIPT_ROOT
# 		script root

echo "Step6: Plot Cytospace distribution"
conda activate seurat4.3
Rscript "${script_root}/6_plot_cytospace_distribution.r" \
    -s "${ST_sample_name}" \
    -c "${out_path}/${step_names[4]}/Cytospace_seurat_object.qs" \
    -o "${out_path}/${step_names[5]}" \
    -l "${cell_type_column_name}" \
    -p 1 \
    -a 0.1 \
    -t FALSE \
    -r "${script_root}"