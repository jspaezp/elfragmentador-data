from mokapot_utils import parse_mokapot_log
from mokapot_helpers import get_mokapot_command


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
