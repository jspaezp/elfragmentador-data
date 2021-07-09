
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
    crap_fasta,
    decoy_db,
    biognosys_irt_fasta,
    human_fasta,
    yeast_fasta,
    human_yeast_fasta,
    contam_fasta,
    added_irt_fasta,
    proteometools_fasta,
    # download_file,
    # convert_file,
    # mzml_scan_metadata,
    comet_phospho_params,
    comet_gg_params,
    comet_proalanase_params,
    # comet_search,
    experiment_fasta,
    # mokapot,
    mokapot_spectrast_in,
    # split_mokapot_spectrast_in,
    # mokapot_spectrast,
    # interact,
    # peptideprophet,
    # indiv_spectrast,
    # iprophet,
    # spectrast,
    # generate_sptxt_csv,
    prosit_input,
    aggregate_mokapot_sptxts,

rule all:
    input:
        # [f"ind_spectrast/{x}.pp.sptxt" for x in samples["sample"]],
        # [f"spectrast/concensus_{x}.iproph.pp.sptxt" for x in samples["experiment"]],
        # [f"prosit_in/{x}.iproph.pp.sptxt" for x in samples["experiment"]],
        # [f"sptxt_csv/{x}.iproph.pp.sptxt.csv" for x in samples["experiment"]],
        # dynamic([f"concensus_mokapot_spectrast/concensus_{x}.{{ce}}.mokapot.psms.spidx" for x in np.unique(samples["experiment"])]),
        [f"aggregated/{experiment}/aggregated_concensus_{experiment}.mokapot.sptxt" for experiment in np.unique(samples["experiment"])],
        [f"raw_scan_metadata/{sample}.csv" for sample in samples["sample"]],
        [f"comet/{sample}.png" for sample in samples["sample"]],


rule crap_fasta:
    output:
        fasta="fasta/crap.fasta",
    shell:
        """
        mkdir -p fasta
        wget \
            --timeout=15 \
            --limit-rate=50m \
            --wait=5 \
            ftp://ftp.thegpm.org/fasta/cRAP/crap.fasta \
            -O ./fasta/crap.fasta
        """

rule decoy_db:
    input:
        "fasta/{file}.fasta"
    output:
        "fasta/{file}.decoy.fasta"
    run:
        # Note that it writes a concatenated decoy+target database, not only decoys
        fasta.write_decoy_db(str(input), str(output))


rule biognosys_irt_fasta:
    output:
        "fasta/irtfusion.fasta"
    shell:
        """
        mkdir -p fasta
        wget \
            https://biognosys.com/media.ashx/irtfusion.fasta \
            -O fasta/irtfusion.fasta
        """


rule human_fasta:
    output:
        "fasta/human.fasta"
    shell:
        """
        mkdir -p fasta
        wget \
            https://www.uniprot.org/uniprot/\?query\=proteome:UP000005640%20reviewed:yes\&format\=fasta \
            -O fasta/human.fasta
        """


rule yeast_fasta:
    output:
        "fasta/yeast.fasta"
    shell:
        """
        mkdir -p fasta
        wget \
            https://www.uniprot.org/uniprot/\?query\=proteome:UP000002311%20reviewed:yes\&format\=fasta \
            -O fasta/yeast.fasta
        """


rule human_yeast_fasta:
    input:
        "fasta/human.fasta",
        "fasta/yeast.fasta"
    output:
        "fasta/human_yeast.fasta"
    shell:
        """
        cat fasta/human.fasta fasta/yeast.fasta > fasta/human_yeast.fasta
        """

rule contam_fasta:
    input:
        "fasta/crap.fasta",
        "fasta/{basename}.fasta"
    output:
        "fasta/{basename}_contam.fasta"
    shell:
        """
        cat fasta/{wildcards.basename}.fasta fasta/crap.fasta > fasta/{wildcards.basename}_contam.fasta
        """

rule added_irt_fasta:
    input:
        "fasta/irtfusion.fasta",
        "fasta/{basename}.fasta"
    output:
        "fasta/{basename}_irt.fasta"
    shell:
        """
        cat fasta/{wildcards.basename}.fasta fasta/irtfusion.fasta > fasta/{wildcards.basename}_irt.fasta
        """

rule proteometools_fasta:
    input:
        "fasta/RT_QC.fasta",
        "fasta/crap.fasta",
        "fasta/proteometools/Packet_{basename}.fasta"
    output:
        "fasta/{basename}_contam.fasta"
    shell:
        """
        cat fasta/proteometools/Packet_{wildcards.basename}.fasta \
            fasta/RT_QC.fasta \
            fasta/crap.fasta > fasta/{wildcards.basename}_contam.fasta
        """

rule download_file:
    output:
        "raw/{sample}.raw"
    run:
        server = samp_to_ftp[wildcards.sample]
        shell("mkdir -p raw")
        shell((
            f"wget {server}"+"/{wildcards.sample}.raw -O "
            "{output} "
            "--timeout=15 "
            "--limit-rate=50m "
            "--wait=5 "
        ))
    

rule convert_file:
    input:
        "raw/{sample}.raw"
    output:
        "raw/{sample}.mzML"
    benchmark:
        "benchmarks/{sample}.conversion.benchmark.txt"
    run:
        # For some reaso, the dockerized version fails when running it directly
        # in this script, so you have to hack it this way ...
        subprocess.run(['zsh', 'msconvert.bash', str(input)])


rule mzml_scan_metadata:
    input:
        "raw/{sample}.mzML"
    output:
        "raw_scan_metadata/{sample}.csv"
    run:
        shell("mkdir -p raw_scan_metadata")
        df = get_spec_metadata(str(input))
        df.to_csv(str(output), index = False)


def get_fasta(wildcards):
    return samp_to_fasta[wildcards.sample]

def get_comet_params(wildcards):
    return samp_to_params[wildcards.sample]

rule comet_phospho_params:
    input:
        "comet_params/comet.params.high_high"
    output:
        "comet_params/comet_phospho.params.high_high"
    shell:
        """
        cat {input} | \
            sed -e "s/variable_mod02 = 0.0 X 0 3 -1 0 0/variable_mod02 = 79.966331 STY 0 3 -1 0 0/g" \
            | tee {output}
        """

rule comet_gg_params:
    input:
        "comet_params/comet.params.high_high"
    output:
        "comet_params/comet_gg.params.high_high"
    shell:
        """
        cat {input} | \
            sed -e "s/variable_mod02 = 0.0 X 0 3 -1 0 0/variable_mod02 = 114.042927 K 0 3 -1 0 0/g" \
            | tee {output}
        """

rule comet_proalanase_params:
    input:
        "comet_params/comet.params.high_high"
    output:
        "comet_params/comet.params.proalanase.high_high"
    shell:
        """
        cat {input} | \
            sed -e "s/^search_enzyme_number.*/search_enzyme_number = 10/g" | \
            sed -e "s/^10. Chymotrypsin.*/10. ProAlanase 1 PA -/g" \
            | tee {output}
        """

rule comet_search:
    input:
        raw="raw/{sample}.mzML",
        fasta=get_fasta,
        comet_params=get_comet_params,
    output:
        # Enable if using 2 in the decoy search parameter
        # decoy_pepxml = "comet/{sample}.decoy.pep.xml", 
        forward_pepxml = "comet/{sample}.pep.xml",
        forward_pin = "comet/{sample}.pin",
    benchmark:
        "benchmarks/{sample}.comet_search.benchmark.txt"
    run:
        shell("mkdir -p comet")
        cmd = (
            f"{TPP_DOCKER} "
            f"comet -P{str(input.comet_params)} "
            f"-D{str(input.fasta)} "
            f"{str(input.raw)} " )

        print(cmd)
        shell(cmd)
        shell(f"cp raw/{wildcards.sample}.pep.xml ./comet/.")
        shell(f"cp raw/{wildcards.sample}.pin ./comet/.")


rule comet_decoy_plot:
    input:
        forward_pepxml = "comet/{sample}.pep.xml",
        forward_pin = "comet/{sample}.pin",
    output:
        decoy_plot = "comet/{sample}.png"
    run:
        # Plotting the xcorr distribution of decoys and targets
        df = pd.read_csv(
            f"comet/{wildcards.sample}.pin",
            index_col=False,
            error_bad_lines=False,
            sep = "\t",
            usecols=['Xcorr', 'Label'],
            lineterminator='\n')

        print(df)

        tps = []
        fps = []
	tps.extend([s for s,l in zip(df['Xcorr'], df['Label']) if l > 0])
	fps.extend([s for s,l in zip(df['Xcorr'], df['Label']) if l < 0])

        plt.hist(tps, bins=20, alpha=0.8, color = 'cyan')
        plt.hist(fps, bins=20, alpha=0.8, color = 'magenta')
        plt.savefig(f"comet/{wildcards.sample}.png",)


def get_mokapot_ins(wildcards):
    outs = expand("comet/{sample}.pin", sample = get_samples(wildcards.experiment))
    return outs

def get_exp_fasta(wildcards):
    return exp_to_fasta[wildcards.experiment]

def get_enzyme_regex(wildcards):
    out = exp_to_enzyme[wildcards.experiment]
    assert len(out) == 1
    return out[0]

rule experiment_fasta:
    input:
        get_exp_fasta
    output:
        "experiment_fasta/{experiment}.fasta"
    shell:
        """
        set -x
        set -e

        mkdir -p experiment_fasta
        cat {input} > experiment_fasta/{wildcards.experiment}.fasta
        """

rule mokapot:
    input:
        pin_files=get_mokapot_ins,
        fasta="experiment_fasta/{experiment}.fasta"
    output:
        "mokapot/{experiment}.mokapot.psms.txt",
        "mokapot/{experiment}.mokapot.peptides.txt"
    params:
        enzyme_regex=get_enzyme_regex,
    benchmark:
        "benchmarks/{experiment}.mokapot.benchmark.txt"
    shell:
        """
        set -e
        set -x

        mkdir -p mokapot
        mokapot --verbosity 2 \
            --seed 2020 \
            --aggregate \
            --enzyme {params.enzyme_regex} \
            --decoy_prefix DECOY_ \
            --proteins {input.fasta} \
            --missed_cleavages 2 \
            --min_length 5 \
            --max_length 50 \
            -d mokapot \
            -r {wildcards.experiment} {input.pin_files} 
        """

rule mokapot_spectrast_in:
    input:
        psm_input="mokapot/{experiment}.mokapot.psms.txt",
        peptide_input="mokapot/{experiment}.mokapot.peptides.txt",
    output:
        "spectrast_in/{experiment}.spectrast.mokapot.psms.tsv",
    run:
        shell("mkdir -p spectrast_in")
        df = filter_mokapot_psm(input.psm_input, input.peptide_input)
        psm_df_to_tsv(df, str(output))

def get_exp_spec_metadata(wildcards):
    samples = exp_to_sample[wildcards.experiment]
    out = ["raw_scan_metadata/" + sample + ".csv" for sample in samples]
    return out 

# Note that this is a checkpoint, not a rule
checkpoint split_mokapot_spectrast_in:
    input:
        spectrast_in = "spectrast_in/{experiment}.spectrast.mokapot.psms.tsv",
        spec_metadata = get_exp_spec_metadata
    output:
        split_files=directory("split_spectrast_in/{experiment}")
    run:
        shell(f"mkdir -p {output.split_files}")
        print(input.spec_metadata)
        split_mokapot_spectrast_in(
            input.spectrast_in,
            input.spec_metadata,
            wildcards.experiment)
        

rule mokapot_spectrast:
    input:
        spectrast_usermods="spectrast_params/spectrast.usermods",
        mokapot_in="split_spectrast_in/{experiment}/{experiment}.{n}.spectrast.mokapot.psms.tsv"
    output:
        "mokapot_spectrast/{experiment}/{experiment}.{n}.mokapot.psms.sptxt",
        "mokapot_spectrast/{experiment}/{experiment}.{n}.mokapot.psms.splib",
        "mokapot_spectrast/{experiment}/{experiment}.{n}.mokapot.psms.pepidx",
        "mokapot_spectrast/{experiment}/{experiment}.{n}.mokapot.psms.spidx",
    benchmark:
        "benchmarks/{experiment}.{n}.mokapot_spectrast.benchmark.txt"
    shell:
        "set -x ; set -e ; mkdir -p mokapot_spectrast ; "
        f"{TPP_DOCKER}"
        " spectrast -V -cIHCD -M{input.spectrast_usermods}"
        " -Lmokapot_spectrast/{wildcards.experiment}/{wildcards.experiment}.{wildcards.n}.mokapot.psms.log"
        " -cNmokapot_spectrast/{wildcards.experiment}/{wildcards.experiment}.{wildcards.n}.mokapot.psms"
        " {input.mokapot_in} ;"


rule mokapot_spectrast_concensus:
    input:
        spectrast_usermods="spectrast_params/spectrast.usermods",
        splib="mokapot_spectrast/{experiment}/{experiment}.{n}.mokapot.psms.splib",
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
        " spectrast -V -cr1 -cIHCD -cAC -c_DIS -M{input.spectrast_usermods}" 
        " -Lconcensus_mokapot_spectrast/concensus_{wildcards.experiment}.{wildcards.n}.mokapot.psms.log"
        " -cNconcensus_mokapot_spectrast/concensus_{wildcards.experiment}.{wildcards.n}.mokapot.psms"
        " {input.splib}"


import os

def aggregate_input(wildcards):
    checkpoint_output = checkpoints.split_mokapot_spectrast_in.get(**wildcards).output[0]
    globbing_path = os.path.join(
            checkpoint_output,
            f"{wildcards.experiment}" + ".{n}.spectrast.mokapot.psms.tsv")
    print(globbing_path)
    globbed_wildcards = glob_wildcards(globbing_path)

    print(globbed_wildcards)
    globbed_wildcards = globbed_wildcards.n

    print(globbed_wildcards)

    out = expand(
        "concensus_mokapot_spectrast/concensus_{experiment}.{n}.mokapot.psms.sptxt",
         experiment=wildcards.experiment,
         n=globbed_wildcards)


    return out


rule aggregate_mokapot_sptxts:
    input:
        aggregate_input
    output:
        "aggregated/{experiment}/aggregated_concensus_{experiment}.mokapot.sptxt",
    shell:
        "cat {input} > {output}"


rule interact:
    input:
        # rv_file="comet/{sample}.decoy.pep.xml",
        fw_file="comet/{sample}.pep.xml",
        fasta_file=get_fasta
    output:
        "pp/{sample}.pep.xml"
    shell:
        "set -x ; set -e ; mkdir -p pp ; "
        f"{TPP_DOCKER} "
        "xinteract -G -N{output} -nP "
        "{input.fw_file}"
        # "{input.rv_file} {input.fw_file}"

rule peptideprophet:
    input:
        "pp/{sample}.pep.xml"
    output:
        "pp/{sample}.pp.pep.xml"
    shell:
        """
        set -e
        set -x
        cp {input} {output}
        
        """ + 
        f"{TPP_DOCKER} " +
        "PeptideProphetParser {output} ACCMASS DECOY=DECOY_ DECOYPROBS NONPARAM"

rule indiv_spectrast:
    input:
        "pp/{sample}.pp.pep.xml"
    output:
        "ind_spectrast/{sample}.pp.sptxt",
        "ind_spectrast/{sample}.pp.splib",
        "ind_spectrast/{sample}.pp.pepidx",
        "ind_spectrast/{sample}.pp.spidx"
    shell:
        "set -x ; set -e ; mkdir -p ind_spectrast ; "
        f"{TPP_DOCKER}"
        " spectrast -c_RDYDECOY_ -cP0.9 -cq0.01 -cIHCD" # " -Mspectrast.usermods"
        " -Lind_spectrast/{wildcards.sample}.pp.log"
        " -cNind_spectrast/{wildcards.sample}.pp {input}"

def get_iproph_ins(wildcards):
    outs = expand("pp/{sample}.pp.pep.xml", sample = get_samples(wildcards.experiment))
    return outs

rule iprophet:
    input: get_iproph_ins
    output:
        "ip/{experiment}.iproph.pp.pep.xml"

    shell:
        "set -x ; set -e ; "
        f"{TPP_DOCKER}"
        "InterProphetParser THREADS=4 DECOY=DECOY NONSS NONSE"
        " {input} {output}"

# TODO check how to fix the issue where using interprophet destrys the model ...
rule spectrast:
    input:
        "ip/{experiment}.iproph.pp.pep.xml"
    output:
        "spectrast/{experiment}.iproph.pp.sptxt",
        "spectrast/{experiment}.iproph.pp.splib",
        "spectrast/{experiment}.iproph.pp.pepidx",
        "spectrast/{experiment}.iproph.pp.spidx",
        "spectrast/concensus_{experiment}.iproph.pp.sptxt",
        "spectrast/concensus_{experiment}.iproph.pp.splib",
        "spectrast/concensus_{experiment}.iproph.pp.pepidx",
        "spectrast/concensus_{experiment}.iproph.pp.spidx"
    shell:
        "set -x ; set -e ; mkdir -p spectrast ; "
        f"{TPP_DOCKER}"
        " spectrast -c_RDYDECOY_ -cP0.9 -cq0.01 -cIHCD" # " -Mspectrast.usermods"
        " -Lspectrast/{wildcards.experiment}.iproph.pp.log"
        " -cNspectrast/{wildcards.experiment}.iproph.pp {input} ;"
        f"{TPP_DOCKER}"
        " spectrast -cr1 -cAC -c_DIS " 
        " -Lspectrast/concensus_{wildcards.experiment}.iproph.pp.log"
        " -cNspectrast/concensus_{wildcards.experiment}.iproph.pp "
        " spectrast/{wildcards.experiment}.iproph.pp.splib"

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
