#!/bin/bash

# Required filepaths and user input

read -r -p "Please enter filepath from the 'OPA Results' folder. TODO: provide example " run_input
drive_filepath="\\\cuh_nas120\Haematology\Molecular Haemato-Oncology\A03_Test_Solid\OPA on Genexus\OPA Results\\${run_input}"
opa_id=$(basename "${drive_filepath}")

read -r -p "Please enter filepath for the Sample Sheet CSV from the 'OPA' folder in the ClinGen drive. TODO: provide example " sample_sheet_input
sample_sheet="\\\clingen\CG\Regional Genetics Laboratories\Molecular Genetics\Data archive\Sequencing HT\Genexus\OPA\\${sample_sheet_input}"

# TODO: store somewhere and use file path
opa_genelist="OPA_GeneList.csv"

# Error checks for user input -- required

header=$(head -n 1 "$sample_sheet")
# Get first and last parts of OPA run name
first=$(cut -d_ -f1 <<< $opa_id)
last=${opa_id##*_}

if [[ ! -d "$drive_filepath" ]] || [[ ! -f "$sample_sheet" ]]; then
    echo -e "The filepath to the OPA samples and/or the sample sheet does not exist. Please enter a valid filepath."
    exit 1
elif [[ ! "$sample_sheet" == *.csv ]]; then
    echo "Sample sheet is not a csv."
    exit 1
#should this be case insensitive?
elif ! echo $header | grep -qwi "sample name" || ! echo $header | grep -qwi "cancer type"; then
    echo "Sample sheet is missing required headings: 'Sample Sheet' and 'Cancer Type'."
    exit 1
elif [[ ! "$sample_sheet" == *"$first"* ]] || [[ ! "$sample_sheet" == *"$last"* ]]; then
    echo "Name of sample sheet file does not match OPA run name. Please check again."
    exit 1
fi

# Find all OPA sample directories
sp_folders=$(find "$drive_filepath" -mindepth 1 -maxdepth 1 -type d -printf '%f\n')

# Warning messages and user confirmation
awk -F',' -v find_output="$sp_folders" 'NR > 1 { if ($1 != "" && index(find_output, $1) == 0) { print "Sample not found:", $1 } else { print "Sample found:", $1 } }' "$sample_sheet"
read -r -p "Batch script will only include found samples. Would you like to continue (y/n)? " proceed
case "$proceed" in 
  y|Y) echo "Generating batch script...";;
  n|N ) exit 0;;
  *) echo "Please enter either 'y' or 'n'";;
esac

# Get column number for sample name and cancer type
# depends on the answer to case insensitive above
sample_name_field=$(awk -v RS=',' '/Sample Name/{print NR}' "$sample_sheet")
cancer_type_field=$(awk -v RS=',' '/Cancer Type/{print NR}' "$sample_sheet")

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
    exon_arr=("$@")
    for entry in "${exon_arr[@]}";
    do
        IFS='_' read -ra gene_exon <<< "$entry"
        gene="${gene_exon[0]}"
        exon="${gene_exon[1]}"
        chr="$(awk -F, -v val1=$gene -v val2=$exon '$4 == val1 && $6 == val2 {print $1}' "$opa_genelist")"
        start="$(awk -F, -v val1=$gene -v val2=$exon '$4 == val1 && $6 == val2 {print $2}' "$opa_genelist")"
        stop="$(awk -F, -v val1=$gene -v val2=$exon '$4 == val1 && $6 == val2 {print $3}' "$opa_genelist")"
        start_flanked=$((start+1-flanking))
        stop_flanked=$((stop+flanking))
        echo goto chr${chr}:${start_flanked}-${stop_flanked}
        echo snapshot ${gene}_Ex${exon}.png
    done >> ${opa_id}_IGV.batch
}

for sp_id in $sp_folders
do
    bam_file="${drive_filepath}\\${sp_id}\\merged.bam.ptrim.bam"

    # Create header for IGV batch script and append
    batch_beginning="snapshotDirectory '\\\\${drive_filepath}\\${sp_id}\\snapshots'\nnew\nload '\\\\$bam_file'\nsort position\ncollapse\ncolorBy READ_STRAND\npreference SAM.SHOW_SOFT_CLIPPED true"
    echo "$(echo -e $batch_beginning)" >> ${opa_id}_IGV.batch

    # Ignore any extensions after _ if present (for SP-ID_rpt samples)
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
            echo 'Tumour type not found. Outputting no coordinates for' $sp_id
            continue
            ;;
    esac
    create_snapshots "${gene_exon_pair[@]}"
done

echo "Batch script generated."