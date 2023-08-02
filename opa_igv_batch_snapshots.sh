#!/bin/bash

# Get user input

read -r -p "Please enter filepath from the 'OPA Results' folder. TODO: provide example " run_input
# for testing, input would be "others\OPA_IGV_batch_snapshots_testing\23-OPAD19_B_NP_20062023_2"
drive_filepath="\\\cuh_nas120\Haematology\Molecular Haemato-Oncology\A03_Test_Solid\OPA on Genexus\OPA Results\\${run_input}"
# TODO: error checks: filepath is valid/exists, folder structure (sample folders with bam files)
opa_id=$(basename "${drive_filepath}")

read -r -p "Please enter filepath for the Sample Sheet CSV from the 'OPA' folder in the ClinGen drive. TODO: provide example " sample_sheet_input
# for testing, input would be "23-OPAD19_23-OPAR4_20062023\23-OPAD19_23-OPAR4_SampleSheet_20062023.csv"
sample_sheet="\\\clingen\CG\Regional Genetics Laboratories\Molecular Genetics\Data archive\Sequencing HT\Genexus\OPA\\${sample_sheet_input}"
# TODO: error checks: filepath/file is valid/exists, file is a csv, check for sample sheet headers

# TODO: not yet made dynamic, CSV needs to be finalised
opa_genelist="OPA_GeneList.csv"


####################
# Creates the goto and snapshot commands of the IGV batch file
# Arguments:
#   Array of gene exons from OPA_GeneList.csv
# Outputs:
#   IGV batch file
####################
create_snapshots () {
    gene_arr=("$@")
    for exon in "${gene_arr[@]}";
    do
        chr="$(awk -F, -v val=$exon '$4 == val {print $1}' "$opa_genelist")"
        start="$(awk -F, -v val=$exon '$4 == val {print $2}' "$opa_genelist")"
        stop="$(awk -F, -v val=$exon '$4 == val {print $3}' "$opa_genelist")"
        echo goto ${chr}:${start}-${stop}
        echo snapshot ${exon}.png
    done >> ${opa_id}_IGV.batch
}


# find all OPA sample directories from user input
sp_folders=$(find "$drive_filepath" -mindepth 1 -maxdepth 1 -type d -printf '%f\n')

for sp_id in $sp_folders
do
    bam_file="${drive_filepath}\\${sp_id}\\merged.bam.ptrim.bam"

    # create header for IGV batch script and append
    batch_beginning="snapshotDirectory '\\\\${drive_filepath}\\${sp_id}\\snapshots'\nnew\nload '\\\\$bam_file'\nsort position\ncollapse\ncolorBy READ_STRAND"
    echo "$(echo -e $batch_beginning)" >> ${opa_id}_IGV.batch

    # get tumour type for each sample
    tumour_type="$(awk -F, -v val=$sp_id '$1 == val {print $8}' "$sample_sheet")"

    if [[ ${tumour_type,,} == "colorectal cancer" ]]; then
        genes=("BRAF_Ex11" "BRAF_Ex15" "KRAS_Ex2" "KRAS_Ex3" "KRAS_Ex4" "NRAS_Ex2" "NRAS_Ex3" "NRAS_Ex4")
        create_snapshots "${genes[@]}"
    # TODO: incomplete regions for breast cancer, it's missing PIK3CA ex 21. Also fix bug for test sample that doesn't read breast cancer
    elif [[ ${tumour_type,,} == "breast cancer" ]]; then 
        genes=("PIK3CA_Ex2" "PIK3CA_Ex3" "PIK3CA_Ex5" "PIK3CA_Ex8" "PIK3CA_Ex10" "PIK3CA_Ex11" "PIK3CA_Ex20")
        create_snapshots "${genes[@]}"
    elif [[ ${tumour_type,,} == "lung cancer" || ${tumour_type,,} == "non-small cell lung cancer" ]]; then
        genes=("EGFR_Ex18" "EGFR_Ex19" "EGFR_Ex20" "EGFR_Ex21" "MET_Ex14")
        create_snapshots "${genes[@]}"
    elif [[ ${tumour_type,,} == "melanoma" ]]; then
        genes=("BRAF_Ex11" "BRAF_Ex15" "NRAS_Ex2" "NRAS_Ex3" "NRAS_Ex4" "KIT_Ex8" "KIT_Ex9" "KIT_Ex10" "KIT_Ex11" "KIT_Ex13" "KIT_Ex14" "KIT_Ex17" "KIT_Ex18")
        create_snapshots "${genes[@]}"
    elif [[ ${tumour_type,,} == "gist" ]]; then
        genes=("BRAF_Ex11" "BRAF_Ex15" "KIT_Ex8" "KIT_Ex9" "KIT_Ex10" "KIT_Ex11" "KIT_Ex13" "KIT_Ex14" "KIT_Ex17" "KIT_Ex18" "PDGFRA_Ex11" "PDGFRA_Ex12" "PDGFRA_Ex14" "PDGFRA_Ex18")
        create_snapshots "${genes[@]}"
    elif [[ ${tumour_type,,} == "brain cancer" ]]; then
        genes=("IDH1_Ex4" "IDH1_Ex6" "IDH1_Ex7" "IDH2_Ex4" "IDH2_Ex7")
        create_snapshots "${genes[@]}"
    else
        echo 'Tumour type not found. Outputting all coordinates.'
        all_genes_arr=( $( cut -d ',' -f4 OPA_GeneList.csv ) )
        create_snapshots "${all_genes_arr[@]}"
    fi
done