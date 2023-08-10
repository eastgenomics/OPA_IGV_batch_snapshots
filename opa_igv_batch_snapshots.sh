#!/bin/bash

# Required filepaths and user input

read -rep $'Please enter file path to the OPA run folder from the "OPA Results" folder\ne.g "OPA 2023\\23-OPAD12_A_NP_02062023"\nFile path: ' run_input
opa_folder_filepath="\\\cuh_nas120\Haematology\Molecular Haemato-Oncology\A03_Test_Solid\OPA on Genexus\OPA Results\\${run_input}"
opa_id_full=$(basename "${opa_folder_filepath}")

read -rep $'Please enter file path for the Sample Sheet CSV from the "OPA" folder in the ClinGen drive\ne.g "23-OPAD12_02062023\\23-OPAD12_SampleSheet_02062023.csv"\nFile path: ' sample_sheet_input
sample_sheet="\\\clingen\CG\Regional Genetics Laboratories\Molecular Genetics\Data archive\Sequencing HT\Genexus\OPA\\${sample_sheet_input}"

opa_genelist="\\\cuh_nas120\Haematology\Molecular Haemato-Oncology\A03_Test_Solid\OPA on Genexus\OPA Results\OPA_IGV_batch_snapshots\OPA_GeneList_v1.bed"

# Error checks for user input -- required

sample_sheet_header=$(head -n 1 "$sample_sheet")
# Get first (opa_id) and last (opa_date) parts of OPA run name
opa_id=$(cut -d_ -f1 <<< $opa_id_full)
opa_date=${opa_id_full##*_}
# File will be named after OPA name with '_IGV.batch' extension
opa_batch_file="${opa_folder_filepath}\\${opa_id_full}_IGV.batch"

if [[ ! -d "$opa_folder_filepath" ]] || [[ ! -f "$sample_sheet" ]]; then
    echo -e "The filepath to the OPA samples and/or the sample sheet does not exist. Please enter a valid filepath."
    exit 1
elif [[ ! "$sample_sheet" == *.csv ]]; then
    echo "Sample sheet is not a csv."
    exit 1
elif ! echo $sample_sheet_header | grep -qwi "sample name" || ! echo $sample_sheet_header | grep -qwi "cancer type"; then
    echo "Sample sheet is missing required headings: 'Sample Sheet' and 'Cancer Type'."
    exit 1
elif [[ ! "$sample_sheet" == *"$opa_id"* ]] || [[ ! "$sample_sheet" == *"$opa_date"* ]]; then
    echo "Name of sample sheet file does not match OPA run name. Please check again."
    exit 1
elif [[ -f "$opa_batch_file" ]]; then
    echo "An IGV batch script already exists for this OPA run in" $run_input ". Please remove before creating a new batch script."
    exit 1
fi

# Find all OPA sample directories
sp_folders=$(find "$opa_folder_filepath" -mindepth 1 -maxdepth 1 -type d -printf '%f\n')

# Warning messages and user confirmation
samples_not_found=$(awk -F',' -v find_output="$sp_folders" 'NR > 1 { if ($1 != "" && index(find_output, $1) == 0) { print "Sample not found:", $1 } }' "$sample_sheet")
if [[ -n "$samples_not_found" ]]; then
    echo "$samples_not_found"
    while true; do
        read -r -p "Would you like to continue (y/n)? " proceed
        case "$proceed" in 
            y|Y) echo "Generating batch script..."; break;;
            n|N ) exit 0;;
            *) echo "Please enter either 'y' or 'n'";;
        esac
    done
fi

# Get column number for sample name and cancer type
sample_name_field=$(awk -v RS=',' -v IGNORECASE=1 '/Sample Name/{print NR}' "$sample_sheet")
cancer_type_field=$(awk -v RS=',' -v IGNORECASE=1 '/Cancer Type/{print NR}' "$sample_sheet")

# Set flanking sequences value
flanking=15

####################
# Creates the goto and snapshot commands of the IGV batch file
# Arguments:
#   Array of gene exons from OPA_GeneList.csv
# Outputs:
#   IGV batch file
####################
create_snapshots () {
    sample_id=$1
    shift
    exon_arr=("$@")
    for entry in "${exon_arr[@]}";
    do
        IFS='_' read -ra gene_exon <<< "$entry"
        gene="${gene_exon[0]}"
        exon="${gene_exon[1]}"
        chr="$(awk -F '\t' -v gene=$gene -v exon=$exon '$4 == gene && $6 == exon {print $1}' "$opa_genelist")"
        start="$(awk -F '\t' -v gene=$gene -v exon=$exon '$4 == gene && $6 == exon {print $2}' "$opa_genelist")"
        stop="$(awk -F '\t' -v gene=$gene -v exon=$exon '$4 == gene && $6 == exon {print $3}' "$opa_genelist")"
        start_flanked=$((start+1-flanking))
        stop_flanked=$((stop+flanking))
        echo goto chr${chr}:${start_flanked}-${stop_flanked}
        echo snapshot ${sample_id}_${gene}_Ex${exon}.png
    done >> ${opa_id_full}_IGV.batch
}

for sp_id in $sp_folders
do
    bam_file="${opa_folder_filepath}\\${sp_id}\\merged.bam.ptrim.bam"

    # Create header for IGV batch script and append
    batch_beginning="snapshotDirectory '\\\\${opa_folder_filepath}\\${sp_id}\\snapshots'\nnew\nload '\\\\$bam_file'\nsort position\ncollapse\ncolorBy READ_STRAND\npreference SAM.SHOW_SOFT_CLIPPED true"
    echo "$(echo -e $batch_beginning)" >> ${opa_id_full}_IGV.batch

    # Ignore any extensions in folder name after "_" if present
    # This is done for sample folders named after SP_ID_rpt, which are repeated runs
    # Allows for tumour type to be specified for these samples as sample names in sample sheet don't have the "_rpt" extension
    sp_id_no_rpt=${sp_id%_*}

    # Get tumour type for each sample
    tumour_type="$(awk -F, -v sp_id=$sp_id_no_rpt -v sample_field=$sample_name_field -v type_field=$cancer_type_field= '$sample_field == sp_id {print $type_field}' "$sample_sheet")"

    case ${tumour_type,,} in
        *"colorectal cancer"*)
            gene_exon_pair=("BRAF_11" "BRAF_15" "KRAS_2" "KRAS_3")
            ;;
        *"breast cancer"*)
            gene_exon_pair=("PIK3CA_2" "PIK3CA_8")
            ;;
        *"non-small cell lung cancer"*)
            gene_exon_pair=("EGFR_18" "EGFR_19" "EGFR_20" "MET_14")
            ;;
        *melanoma*)
            gene_exon_pair=("BRAF_11" "BRAF_15" "KIT_8" "KIT_9" "KIT_10" "KIT_11")
            ;;
        *gist*)
            gene_exon_pair=("BRAF_11" "BRAF_15" "KIT_8" "KIT_9" "KIT_10" "KIT_11" "PDGFRA_11" "PDGFRA_12" "PDGFRA_18")
            ;;
        *)
            echo 'Tumour type' $tumour_type 'assigned to' $sp_id 'does not require IGV images. Skipping.'
            continue
            ;;
    esac
    create_snapshots "$sp_id" "${gene_exon_pair[@]}"
done

# Move batch script to OPA run folder
mv ${opa_id_full}_IGV.batch "$opa_folder_filepath"

echo "Batch script generated."