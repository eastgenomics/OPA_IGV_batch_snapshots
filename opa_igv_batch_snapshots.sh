#!/bin/bash

# get run_id as input
run_id="23-OPAD19_B_NP_20062023_2"
drive_filepath="\\\cuh_nas120\Haematology\Molecular Haemato-Oncology\A03_Test_Solid\OPA on Genexus\OPA Results\others\OPA_IGV_batch_snapshots_testing\\${run_id}"

for sp_folder in ./${run_id}/*/
do
    sp_id=$(basename $sp_folder)
    bam_file="${drive_filepath}\\${sp_id}\\merged.bam.ptrim.bam"
# create IGV batch script with chromosomal regions from OPA genelist
    batch_beginning="snapshotDirectory '\\\\${drive_filepath}\\${sp_id}\\snapshots'\nnew\nload '\\\\$bam_file'"
    echo "$(echo -e $batch_beginning | cat - OPA_GeneList.batch)" >> ${run_id}_IGV.batch
done