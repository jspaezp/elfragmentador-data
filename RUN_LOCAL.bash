#!/bin/bash

for i in {1..1} ; do
    /scratch/brown/jpaezpae/conda_envs/ef_data3/bin/poetry \
        run snakemake \
        --reason \
        --show-failed-logs \
        --printshellcmds \
        --cores 1 \
        -s workflow.smk \
        --config tsv_file=target_files/HLA_groups.tsv \
        --rerun-incomplete ${@}
        # --keep-going \
        # --verbose \
done

# target_files/HLA_groups.tsv \
# test_target_files/small_sample_info.tsv
# test_target_files/small_proteometools_ptm.tsv
# target_files/proteometools_ptm.tsv
# target_files/sample_info.tsv
# target_files/eval_peptidome_sample_info.tsv
