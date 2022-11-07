
from collections import defaultdict
import pandas as pd
import numpy as np
import random
from datetime import datetime
from mokapot_utils import split_mokapot_psms_file

random.seed(42)

# This prevents the command failing without an x server

now = datetime.now()  # current date and time
timestamp = now.strftime("%Y%m%d_%H%M")

IN_TSV = config["tsv_file"]
CHECKPOINT = config.get("checkpoint", None)

print(CHECKPOINT)

samples = pd.read_table(IN_TSV, comment="#")
samples["raw_sample"] = [x for x in samples["sample"]]
samples["sample"] = [x.replace(" ", "") for x in samples["sample"]]
samples=samples.set_index("sample", drop=False)
print(samples)

sample_to_raw = {k:v for k,v in zip(samples["sample"], samples["raw_sample"])}

def get_rawsample(sample):
    return sample_to_raw[sample]


def get_samples(experiment):
    return list(samples[samples["experiment"] == experiment]["sample"])


experiments = {k: v for k, v in zip(samples["experiment"], samples["comet_params"])}

exp_to_enzyme = defaultdict(lambda: list())
exp_to_fasta = defaultdict(lambda: list())
exp_to_sample = defaultdict(lambda: list())

for k, v1, v2, v3 in zip(
    samples["experiment"], samples["enzyme_regex"], samples["fasta"], samples["sample"]
):
    exp_to_enzyme[k].append(v1)
    exp_to_fasta[k].append(v2)
    exp_to_sample[k].append(v3)

exp_to_enzyme = {k: list(set(v)) for k, v in exp_to_enzyme.items()}
exp_to_fasta = {k: list(set(v)) for k, v in exp_to_fasta.items()}

samp_to_fasta = {}
samp_to_params = {}
samp_to_ftp = {}

for samp, fas, params, ftp in zip(
    samples["sample"], samples["fasta"], samples["comet_params"], samples["server"]
):
    samp_to_fasta[samp] = fas
    samp_to_params[samp] = params
    samp_to_ftp[samp] = ftp

def get_exp_spec_metadata(wildcards):
    samples = exp_to_sample[wildcards.experiment]
    out = ["raw_scan_metadata/" + sample + ".csv" for sample in samples]
    return out

# Mokapot related

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

# Search ops related

def get_fasta(wildcards):
    return samp_to_fasta[wildcards.sample]


def get_comet_params(wildcards):
    return samp_to_params[wildcards.sample]



UNIQ_EXP = np.unique(samples["experiment"])

localrules:
    # get list with $ grep -P "^rule" workflow.smk | sed -e "s/rule //g" | sed -e "s/:/,/g" 
    # and comment out the ones that should run in parallel in the cluster
    all,
    download_file,
    crap_fasta,
    decoy_db,
    biognosys_irt_fasta,
    human_fasta,
    yeast_fasta,
    arabidopsis_fasta,
    human_yeast_fasta,
    contam_fasta,
    added_irt_fasta,
    proteometools_fasta,
    proteometools2_fasta,
    experiment_fasta,
    comet_phospho_params,
    comet_gg_params,
    comet_proalanase_params,
    comet_peptidome_params,
    comet_nocleavage_params,
    comet_semitriptic_params,
    clip_comet_pin,
    mokapot_spectrast_in,
    aggregate_mokapot_sptxts,



module fasta_preparations:
    snakefile:
        "./snakemodules/fasta_preparations.smk"

use rule * from fasta_preparations

# module raw_file_operations:
#     snakefile:
#         "./snakemodules/raw_file_operations.smk"
# 
# use rule * from raw_file_operations
include: "./snakemodules/raw_file_operations.smk"

# module search_operations:
#     snakefile:
#         "./snakemodules/search_operations.smk"
# 
# use rule * from search_operations

include: "./snakemodules/search_operations.smk"

module spectrast_operations:
    snakefile:
        "./snakemodules/spectrast_operations.smk"

# module mokapot_operations:
#     snakefile:
#         "./snakemodules/mokapot_operations.smk"
# 
# use rule * from mokapot_operations

include: "./snakemodules/mokapot_operations.smk"

module irt_operations:
    snakefile:
        "./snakemodules/irt_operations.smk"

module train_data_operations:
    snakefile:
        "./snakemodules/train_data_operations.smk"

# module elfragmentador_operations:
#     snakefile:
#         "./snakemodules/elfragmentador_operations.smk"
# use rule * from elfragmentador_operations

include: "./snakemodules/elfragmentador_operations.smk"


# module ptm_operations:
#     snakefile:
#         "./snakemodules/ptm_operations.smk"

include: "./snakemodules/ptm_operations.smk"

module bibliospec_opts:
    snakefile:
        "./snakemodules/bibliospec_operations.smk"

# Note that this is a checkpoint, not a rule
checkpoint split_mokapot_psms:
    input:
        psms = "mokapot/{experiment}.mokapot.psms.txt",
        spec_metadata = get_exp_spec_metadata,
    output:
        split_files_dir=directory("split_mokapot/{experiment}/"),
    run:
        shell(f"mkdir -p {output.split_files_dir}")
        print(input.spec_metadata)
        split_mokapot_psms_file(
            input.psms, input.spec_metadata, output.split_files_dir
        )

use rule bibliospec from bibliospec_opts as bbspec_run with:
    input:
        psms = "split_mokapot/{experiment}/{experiment}.mokapot.psms.NCE{NCE}.txt",
        mzML = get_mokapot_ins("raw/", ".mzML"),
    output:
        ssl_file = "bibliospec/{experiment}.NCE{NCE}.ssl",
        library_name = "bibliospec/{experiment}.NCE{NCE}.blib",


def aggregate_input(wildcards):
    split_psms_output = checkpoints.split_mokapot_psms.get(**wildcards).output[0]
    globbing_path = f"{wildcards.experiment}.mokapot.psms"
    globbing_path = os.path.join(split_psms_output, globbing_path + ".NCE{NCE}.txt")
    glob_wildcards = glob_wildcards(globbing_path)
    glob_wildcards = glob_wildcards.NCE
    out = expand("bibliospec/{experiment}.NCE{NCE}.blib", experiment=wildcards.experiment, NCE=nce_glob)
    return out


rule all_split_blibspecs:
    input:
        aggregate_input,
        psms = "mokapot/{experiment}.mokapot.psms.txt",


common_inputs = [
    [
        f"aggregated/{experiment}/aggregated_concensus_{experiment}.mokapot.sptxt"
        for experiment in UNIQ_EXP
    ],
    [
        f"exp_aggregated_sptxt_csv/{experiment}.mokapot.sptxt.csv"
        for experiment in UNIQ_EXP
    ],
    [
        f"exp_aggregated_rt_sptxt_csv/{experiment}.mokapot.irt.sptxt.csv"
        for experiment in UNIQ_EXP
    ],
    [f"raw_scan_metadata/{sample}.csv" for sample in samples["sample"]],
    [f"comet/{sample}.pin" for sample in samples["sample"]],
    [f"comet/{sample}.xcorr.png" for sample in samples["sample"]],
    [f"comet/{sample}.lnexpect.png" for sample in samples["sample"]],
    [f"rt_csv/{experiment}.irt.csv" for experiment in UNIQ_EXP],
    [
        f"mokapot/{experiment}.mokapot.weights.csv"
        for experiment in UNIQ_EXP
    ],
]

    

all_inputs = [
    [
        f"aggregated_rt_sptxt_csv/{df_set}.mokapot.irt.sptxt.csv"
        for df_set in ("train", "test", "val")
    ],
    [
        f"aggregated_rt_sptxt_feather/{df_set}.mokapot.irt.sptxt.feather"
        for df_set in ("train", "test", "val")
    ],
    [
        f"zip_aggregated_rt_sptxt_csv/{timestamp}.{df_set}.mokapot.irt.sptxt.csv.gz"
        for df_set in ("train", "test", "val")
    ],
]

eval_inputs = [
    [f"ef_comet_pin/{sample}.elfragmentador.pin" for sample in samples["sample"]],
    [
        f"ef_evaluation/{experiment}.csv"
        for experiment in UNIQ_EXP
    ],
    [
        f"ef_evaluation/{experiment}.log"
        for experiment in UNIQ_EXP
    ],
    [
        f"ef_mokapot/{experiment}.elfragmentador.mokapot.psms.txt"
        for experiment in UNIQ_EXP
    ],
    [
        f"ef_mokapot/{experiment}.elfragmentador.mokapot.weights.csv"
        for experiment in UNIQ_EXP
    ],
    [
        f"ef_reports/{experiment}.report.html"
        for experiment in UNIQ_EXP
    ],
    [
        f"ef_reports/{experiment}.roc_curves.html"
        for experiment in UNIQ_EXP
    ],
    [
        f"ef_reports/{experiment}.plot_error_rates.html"
        for experiment in UNIQ_EXP
    ],
    [
        f"ef_reports/{experiment}.plot_error_rates_top1.html"
        for experiment in UNIQ_EXP
    ],
    [
        f"ef_reports/{experiment}.swapped.top.csv"
        for experiment in UNIQ_EXP
    ],
]


print(all_inputs)


rule all:
    input:
        *common_inputs,
        *all_inputs,


rule eval_all:
    input:
        *common_inputs,
        *eval_inputs,


rule get_raw_data:
    input:
        [f"raw/{sample}.raw" for sample in samples["sample"]],

rule get_data:
    input:
        [f"raw/{sample}.raw" for sample in samples["sample"]],
        [f"raw/{sample}.mzML" for sample in samples["sample"]],

rule mokapot_stuff:
    input:
        [
            f"mokapot/{experiment}.mokapot.weights.csv"
            for experiment in UNIQ_EXP
        ],
        [f"comet/{sample}.pin" for sample in samples["sample"]],

rule bibliospec_stuff:
    input:
        [
            f"bibliospec/{experiment}.blib"
            for experiment in UNIQ_EXP
        ],

rule scan_metadata:
    input:
        ["raw_scan_metadata/" + sample + ".csv" for sample in samples["sample"]]
