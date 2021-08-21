from collections import defaultdict
import pandas as pd
import numpy as np
import random
from datetime import datetime

random.seed(42)

# This prevents the command failing without an x server

now = datetime.now()  # current date and time
timestamp = now.strftime("%Y%m%d_%H%M")

in_tsv = config["tsv_file"]
samples = pd.read_table(in_tsv).set_index("sample", drop=False)


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


localrules:
    # get list with $ grep -P "^rule" workflow.smk | sed -e "s/rule //g" | sed -e "s/:/,/g" 
    # and comment out the ones that should run in parallel in the cluster
    all,
    crap_fasta,
    decoy_db,
    biognosys_irt_fasta,
    human_fasta,
    yeast_fasta,
    human_yeast_fasta,
    contam_fasta,
    added_irt_fasta,
    proteometools_fasta,
    comet_phospho_params,
    comet_gg_params,
    comet_proalanase_params,
    experiment_fasta,
    mokapot_spectrast_in,
    prosit_input,
    aggregate_mokapot_sptxts,


UNIQ_EXP = np.unique(samples["experiment"])


include: "./snakemodules/fasta_preparations.smk"
include: "./snakemodules/raw_file_operations.smk"
include: "./snakemodules/search_operations.smk"
include: "./snakemodules/spectrast_operations.smk"
include: "./snakemodules/mokapot_operations.smk"
include: "./snakemodules/irt_operations.smk"
include: "./snakemodules/train_data_operations.smk"
include: "./snakemodules/elfragmentador_operations.smk"


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
    [f"comet/{sample}.xcorr.png" for sample in samples["sample"]],
    [f"comet/{sample}.lnexpect.png" for sample in samples["sample"]],
    [f"rt_csv/{experiment}.irt.csv" for experiment in np.unique(samples["experiment"])],
    [
        f"mokapot/{experiment}.mokapot.weights.csv"
        for experiment in np.unique(samples["experiment"])
    ],
]

all_inputs = [
    [
        f"aggregated_rt_sptxt_csv/{df_set}.mokapot.irt.sptxt.csv"
        for df_set in ("train", "test", "val")
    ],
    [
        f"zip_aggregated_rt_sptxt_csv/{timestamp}.{df_set}.mokapot.irt.sptxt.csv.tar.gz"
        for df_set in ("train", "test", "val")
    ],
]

eval_inputs = [
    [f"ef_comet_pin/{sample}.elfragmentador.pin" for sample in samples["sample"]],
    [
        f"ef_mokapot/{experiment}.elfragmentador.mokapot.psms.txt"
        for experiment in UNIQ_EXP
    ],
    [
        f"ef_mokapot/{experiment}.elfragmentador.mokapot.weights.csv"
        for experiment in np.unique(samples["experiment"])
    ],
    [
        f"ef_reports/{experiment}.report.html"
        for experiment in np.unique(samples["experiment"])
    ],
    [
        f"ef_reports/{experiment}.roc_curves.html"
        for experiment in np.unique(samples["experiment"])
    ],
    [
        f"ef_reports/{experiment}.swapped.top.csv"
        for experiment in np.unique(samples["experiment"])
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


rule prosit_input:
    input:
        "spectrast/{experiment}.iproph.pp.sptxt",
    output:
        "prosit_in/{experiment}.iproph.pp.sptxt",
    shell:
        """
        set -x
        set -e

        mkdir -p prosit_in

        echo 'modified_sequence,collision_energy,precursor_charge' > {output}
        CE=$(grep -oP "CollisionEne.*? " {input} | uniq |  sed -e "s/CollisionEnergy=//g" | sed -e "s/\..*//g")
        grep -P "^Name" {input} | \
            sed -e "s/Name: //g" | \
            sed -e "s+/+,${{CE}},+g" >> {output}

        head {output}
        """
