
from collections import defaultdict
from pathlib import Path
import pandas as pd
import numpy as np
import subprocess
import pathlib
import shutil
from elfragmentador.spectra import sptxt_to_csv
from mokapot_utils import filter_mokapot_psm, psm_df_to_tsv, add_spectrast_ce_info, split_mokapot_spectrast_in
from spec_metadata import get_spec_metadata

# This prevents the command failing without an x server
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

from pyteomics import fasta
from tqdm.auto import tqdm
import mokapot

in_tsv = config['tsv_file']
samples = pd.read_table(in_tsv).set_index("sample", drop=False)

def get_samples(experiment):
    return list(samples[samples['experiment'] == experiment]['sample'])

experiments = {k: v for k, v in zip(samples['experiment'], samples['comet_params'])}

exp_to_enzyme = defaultdict(lambda: list())
exp_to_fasta = defaultdict(lambda: list())
exp_to_sample = defaultdict(lambda: list())

for k, v1, v2, v3 in zip(samples['experiment'], samples['enzyme_regex'], samples['fasta'], samples['sample']):
    exp_to_enzyme[k].append(v1)
    exp_to_fasta[k].append(v2)
    exp_to_sample[k].append(v3)

exp_to_enzyme = {k:list(set(v)) for k, v in exp_to_enzyme.items()}
exp_to_fasta = {k:list(set(v)) for k, v in exp_to_fasta.items()}

samp_to_fasta = {}
samp_to_params = {}
samp_to_ftp = {}

for samp, fas, params, ftp in zip(samples['sample'], samples['fasta'], samples['comet_params'], samples['server']):
    samp_to_fasta[samp] = fas
    samp_to_params[samp] = params
    samp_to_ftp[samp] = ftp

curr_dir = str(pathlib.Path(".").absolute())

if shutil.which("docker"):
    TPP_DOCKER=f"docker run -v {curr_dir}/:/data/ spctools/tpp "
else:
    TPP_DOCKER=f"singularity exec /scratch/brown/jpaezpae/opt/tpp.img "

localrules:
    # get list with $ grep -P "^rule" workflow.smk | sed -e "s/rule //g" | sed -e "s/:/,/g" 
    # and comment out the ones that should run in parallel in the cluster
    all,
    experiment_fasta,
    # mokapot_spectrast,
    # generate_sptxt_csv,
    prosit_input,

split_spectrast_ins = list(Path("split_spectrast_in/").glob("*.*.spectrast.mokapot.psms.tsv"))
outdir = Path("concensus_mokapot_spectrast/")
concensus_spectrast_outs = [
    "concensus_" + str(x.name).replace(".spectrast.mokapot.psms.tsv", ".mokapot.psms.spidx")
    for x in split_spectrast_ins
]
concensus_spectrast_outs = [str(outdir / x) for x in concensus_spectrast_outs]

rule all:
    input:
        concensus_spectrast_outs,


rule mokapot_spectrast:
    input:
        spectrast_usermods="spectrast_params/spectrast.usermods",
        mokapot_in="split_spectrast_in/{experiment}.{n}.spectrast.mokapot.psms.tsv"
    output:
        "mokapot_spectrast/{experiment}.{n}.mokapot.psms.sptxt",
        "mokapot_spectrast/{experiment}.{n}.mokapot.psms.splib",
        "mokapot_spectrast/{experiment}.{n}.mokapot.psms.pepidx",
        "mokapot_spectrast/{experiment}.{n}.mokapot.psms.spidx",
    benchmark:
        "benchmarks/{experiment}.{n}.mokapot_spectrast.benchmark.txt"
    shell:
        "set -x ; set -e ; mkdir -p mokapot_spectrast ; "
        f"{TPP_DOCKER}"
        " spectrast -V -cIHCD -M{input.spectrast_usermods}"
        " -Lmokapot_spectrast/{wildcards.experiment}.{wildcards.n}.mokapot.psms.log"
        " -cNmokapot_spectrast/{wildcards.experiment}.{wildcards.n}.mokapot.psms"
        " {input.mokapot_in} ;"


rule mokapot_spectrast_concensus:
    input:
        spectrast_usermods="spectrast_params/spectrast.usermods",
        splib="mokapot_spectrast/{experiment}.{n}.mokapot.psms.splib",
    output:
        "concensus_mokapot_spectrast/concensus_{experiment}.{n}.mokapot.psms.sptxt",
        "concensus_mokapot_spectrast/concensus_{experiment}.{n}.mokapot.psms.splib",
        "concensus_mokapot_spectrast/concensus_{experiment}.{n}.mokapot.psms.pepidx",
        "concensus_mokapot_spectrast/concensus_{experiment}.{n}.mokapot.psms.spidx"
    benchmark:
        "benchmarks/{experiment}.{n}.mokapot_spectrast_concensus.benchmark.txt"
    shell:
        "set -x ; set -e ; mkdir -p concensus_mokapot_spectrast ; "
        f"{TPP_DOCKER}"
        " spectrast -V -cr3 -cIHCD -cAC -c_DIS -M{input.spectrast_usermods}" 
        " -Lconcensus_mokapot_spectrast/concensus_{wildcards.experiment}.{wildcards.n}.mokapot.psms.log"
        " -cNconcensus_mokapot_spectrast/concensus_{wildcards.experiment}.{wildcards.n}.mokapot.psms"
        " {input.splib}"


rule generate_sptxt_csv:
    input:
        "spectrast/{experiment}.iproph.pp.sptxt",
    output:
        "sptxt_csv/{experiment}.iproph.pp.sptxt.csv"
    run:
        Path(str(output)).parent.mkdir(exist_ok=True)
        sptxt_to_csv(
            filepath=str(input),
            output_path=str(output),
            min_peaks=3,
            min_delta_ascore=20)


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
