from mokapot_utils import parse_mokapot_log


def get_mokapot_ins(base_dir="comet/", extension=".pin"):
    base_str = base_dir + "{sample}" + extension

    def _get_mokapot_ins(wildcards):
        outs = expand(base_str, sample=get_samples(wildcards.experiment))
        return outs

    return _get_mokapot_ins


def get_enzyme_regex(wildcards):
    out = exp_to_enzyme[wildcards.experiment]
    assert len(out) == 1
    return out[0]


def get_mokapot_command(target_dir="mokapot", addition=""):
    base_str = """
        set -e
        set -x

        trap 's=$?; echo >&2 "$0: Error on line "$LINENO": $BASH_COMMAND"; exit $s' ERR

        mkdir -p {target_dir}

        tot_line=$(cat {{input.pin_files}} | wc -l)
        extras=$(python -c "if 5000000 < $tot_line: print(f'--subset_max_train 5000000')")

        mokapot --verbosity 2 \
            --seed 2020 \
            --aggregate \
            --enzyme {{params.enzyme_regex}} \
            --decoy_prefix DECOY_ \
            --missed_cleavages 2 \
            --keep_decoys \
            --min_length 5 \
            --max_length 50 $extras \
            -d {target_dir} \
            -r {{wildcards.experiment}}{addition} \
            {{input.pin_files}} |& tee \
            {target_dir}/{{wildcards.experiment}}{addition}.mokapot.log

            # --subset_max_train 5000000 
            # This does not work right now on files with less than 5M psms....

        """
    out_str = base_str.format(target_dir=target_dir, addition=addition)

    return out_str


rule mokapot:
    input:
        pin_files=get_mokapot_ins("comet/", ".pin"),
    output:
        "mokapot/{experiment}.mokapot.psms.txt",
        "mokapot/{experiment}.mokapot.peptides.txt",
        "mokapot/{experiment}.mokapot.decoy.psms.txt",
        "mokapot/{experiment}.mokapot.decoy.peptides.txt",
        "mokapot/{experiment}.mokapot.log",
    params:
        enzyme_regex=get_enzyme_regex,
    benchmark:
        "benchmarks/{experiment}.mokapot.benchmark.txt"
    shell:
        get_mokapot_command(target_dir="mokapot")


rule mokapot_feature_weights:
    input:
        "mokapot/{experiment}.mokapot.log",
    output:
        "mokapot/{experiment}.mokapot.weights.csv",
    run:
        weights = parse_mokapot_log(str(input))
        weights.to_csv(str(output), index=False)
