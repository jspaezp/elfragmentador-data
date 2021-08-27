
for i in {1..1} ; do
    /scratch/brown/jpaezpae/conda_envs/ef_data3/bin/poetry \
        run snakemake \
        -j 200 --cluster-config cluster.json \
        --cluster "sbatch -A {cluster.account} --mem={cluster.mem} -t {cluster.time} --cpus-per-task={cluster.cpus-per-task} -N {cluster.N}" \
        --local-cores 1 \
        -s workflow.smk \
        --config tsv_file=test_target_files/small_proteometools_ptm.tsv \
        --rerun-incomplete
done

# test_target_files/small_sample_info.tsv
# 
