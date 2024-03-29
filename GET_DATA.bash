#!/bin/bash

for i in {1..1} ; do
    snakemake \
        --reason \
        --show-failed-logs \
        --printshellcmds \
        --cores 2 \
        --keep-going \
        -s workflow.smk \
        --config tsv_file=target_files/TMT_mod_files.tsv \
        --rerun-incomplete ${@} get_data
        # --verbose \
done

# target_files/HLA_groups.tsv \
# test_target_files/small_sample_info.tsv
# test_target_files/small_proteometools_ptm.tsv
# target_files/proteometools_ptm.tsv
# target_files/sample_info.tsv
# target_files/eval_peptidome_sample_info.tsv
