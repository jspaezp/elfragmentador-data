
Run as ...

To display the graph

```
$ snakemake --cores 1 -s workflow.smk --dry-run --config tsv_file=./small_sample_info.tsv --dag | dot -Tsvg > dag.svg
```

When using locally

```
$ snakemake --cores 1 --use-conda --conda-prefix --config tsv_file=./small_sample_info.tsv 
```

When using in slurm:

```
$ `which snakemake` \
    -j 500 \
    --cluster-config cluster.json \
    --cluster "sbatch -A {cluster.account} \
    --mem={cluster.mem} \
    -t {cluster.time} \
    --cpus-per-task={cluster.cpus-per-task} \
    -N {cluster.N}" \
    --use-conda \
    --conda-prefix /scratch/brown/jpaezpae/conda_envs/irt \
    --local-cores 1 \
    -s workflow.smk \
    --printshellcmds \
    --config tsv_file=./small_sample_info.tsv   
```