#!/bin/bash

for i in {1..1} ; do
    poetry run snakemake \
        --reason \
        --show-failed-logs \
        --printshellcmds \
        --cores 4 \
        -s workflow.smk \
        --config tsv_file=target_files/eval_peptidome_sample_info.tsv checkpoint="/home/jspaezp/Downloads/0.50.0b14_onecycle_10e_64_120_val_l=0.141270_epoch=009.ckpt" \
        --rerun-incomplete eval_all ${@}
        # --keep-going \
        # --verbose \
done

#        --config tsv_file=target_files/TMT_mod_files.tsv \
# test_target_files/small_sample_info.tsv
# test_target_files/small_proteometools_ptm.tsv
# target_files/proteometools_ptm.tsv
# target_files/sample_info.tsv
# target_files/eval_peptidome_sample_info.tsv
