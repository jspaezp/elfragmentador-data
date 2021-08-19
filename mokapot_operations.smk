def get_mokapot_ins(wildcards):
    outs = expand("comet/{sample}.pin", sample=get_samples(wildcards.experiment))
    return outs


def get_enzyme_regex(wildcards):
    out = exp_to_enzyme[wildcards.experiment]
    assert len(out) == 1
    return out[0]


rule mokapot:
    input:
        pin_files=get_mokapot_ins,
        fasta="experiment_fasta/{experiment}.fasta",
    output:
        "mokapot/{experiment}.mokapot.psms.txt",
        "mokapot/{experiment}.mokapot.peptides.txt",
    params:
        enzyme_regex=get_enzyme_regex,
    benchmark:
        "benchmarks/{experiment}.mokapot.benchmark.txt"
    shell:
        """
        set -e
        set -x

        mkdir -p mokapot
        # --proteins {input.fasta} \
        mokapot --verbosity 2 \
            --seed 2020 \
            --aggregate \
            --enzyme {params.enzyme_regex} \
            --decoy_prefix DECOY_ \
            --missed_cleavages 2 \
            --min_length 5 \
            --max_length 50 \
            -d mokapot \
            -r {wildcards.experiment} {input.pin_files} 
        """


# TODO check if havinf proteins is required ...
