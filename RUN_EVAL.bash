#!/bin/bash

for i in {1..20} ; do
    /scratch/brown/jpaezpae/conda_envs/ef_data3/bin/poetry \
        run snakemake \
        --reason \
        --show-failed-logs \
        --printshellcmds --skip-script-cleanup --show-failed-logs --verbose \
        --local-cores 2 \
        -j 800 --cluster-config cluster.json \
        --cluster "sbatch -A {cluster.account} --mem={cluster.mem} -t {cluster.time} --cpus-per-task={cluster.cpus-per-task} -N {cluster.N}" \
        -s workflow.smk \
        --config \
          tsv_file=target_files/eval_peptidome_sample_info.tsv \
          checkpoint="https://github.com/jspaezp/elfragmentador-modelzoo/raw/main/0.50.0b14/0.50.0b14_onecycle_10e_64_120_val_l%3D0.141270_epoch%3D009.ckpt" \
        --rerun-incomplete eval_all ${@}
        # --keep-going \
        # --verbose \
    sleep 60
done

#        --config tsv_file=target_files/TMT_mod_files.tsv \
# test_target_files/small_sample_info.tsv
# test_target_files/small_proteometools_ptm.tsv
# target_files/proteometools_ptm.tsv
# target_files/sample_info.tsv
# target_files/eval_peptidome_sample_info.tsv
