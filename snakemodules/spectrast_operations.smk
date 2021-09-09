import os
import pathlib
import pandas as pd
from pandas.errors import EmptyDataError

from elfragmentador.spectra import SptxtReader
from elfragmentador.annotate import canonicalize_seq
from mokapot_utils import (
    filter_mokapot_psm,
    psm_df_to_tsv,
    add_spectrast_ce_info,
    split_mokapot_spectrast_in,
)


include: "./env_setup.smk"
include: "./raw_file_operations.smk"


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


# Note that this is a checkpoint, not a rule
checkpoint split_mokapot_spectrast_in:
    input:
        spectrast_in="spectrast_in/{experiment}.spectrast.mokapot.psms.tsv",
        spec_metadata=get_exp_spec_metadata,
    output:
        split_files=directory("split_spectrast_in/{experiment}"),
    run:
        shell(f"mkdir -p {output.split_files}")
        print(input.spec_metadata)
        split_mokapot_spectrast_in(
            input.spectrast_in, input.spec_metadata, wildcards.experiment
        )


rule mokapot_spectrast:
    input:
        spectrast_usermods="spectrast_params/spectrast.usermods",
        mokapot_in="split_spectrast_in/{experiment}/{experiment}.{n}.spectrast.mokapot.psms.tsv",
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
        "concensus_mokapot_spectrast/concensus_{experiment}.{n}.mokapot.psms.spidx",
    benchmark:
        "benchmarks/{experiment}.{n}.mokapot_spectrast_concensus.benchmark.txt"
    shell:
        "set -x ; set -e ; mkdir -p concensus_mokapot_spectrast ; "
        f"{TPP_DOCKER}"
        " spectrast -V -cr1 -cIHCD -cAC -c_DIS -M{input.spectrast_usermods}"
        " -Lconcensus_mokapot_spectrast/concensus_{wildcards.experiment}.{wildcards.n}.mokapot.psms.log"
        " -cNconcensus_mokapot_spectrast/concensus_{wildcards.experiment}.{wildcards.n}.mokapot.psms"
        " {input.splib} ; "


rule mokapot_sptxt_add_ce:
    input:
        "concensus_mokapot_spectrast/concensus_{experiment}.{n}.mokapot.psms.sptxt",
    output:
        "ce_concensus_mokapot_spectrast/{experiment}.{n}.concensus.ce.mokapot.psms.sptxt",
    run:
        add_string = f"{wildcards.n}".replace("_", ".").replace("UNKNOWN", "NaN")
        add_string = f"CollisionEnergy={add_string}"
        command_string = (
            f"sed -e 's/Comment: /Comment: {add_string} /g' {input} > {output}"
        )
        print(command_string)
        shell(command_string)


def aggregate_input(wildcards):
    checkpoint_output = checkpoints.split_mokapot_spectrast_in.get(**wildcards).output[
        0
    ]
    globbing_path = os.path.join(
        checkpoint_output, f"{wildcards.experiment}" + ".{n}.spectrast.mokapot.psms.tsv"
    )
    print(globbing_path)
    globbed_wildcards = glob_wildcards(globbing_path)

    print(globbed_wildcards)
    globbed_wildcards = globbed_wildcards.n

    print(globbed_wildcards)

    out = expand(
        "ce_concensus_mokapot_spectrast/{experiment}.{n}.concensus.ce.mokapot.psms.sptxt",
        experiment=wildcards.experiment,
        n=globbed_wildcards,
    )

    return out


rule aggregate_mokapot_sptxts:
    input:
        aggregate_input,
    output:
        "aggregated/{experiment}/aggregated_concensus_{experiment}.mokapot.sptxt",
    shell:
        "cat {input} > {output}"


# Defined as ...
#     if the string in the key matches the experiment name...
#     run the command on the value with the sptxt as an argument
ammendments = {
    "Kmod_Trimethyl": 'sed -e "s+K.Acetyl.+K[TRIMETHYL]+g" -i ',
}


# TODO this task is trivially parallelizable ... make it happen ...
rule generate_sptxt_csv:
    input:
        "aggregated/{experiment}/aggregated_concensus_{experiment}.mokapot.sptxt",
    output:
        csv="exp_aggregated_sptxt_csv/{experiment}.mokapot.sptxt.csv",
        feather="exp_aggregated_sptxt_feather/{experiment}.mokapot.sptxt.feather",
    run:
        for k, v in ammendments.items():
            if k in str(input):
                shell(f"{v} {str(input)}")

        pathlib.Path(str(output.csv)).parent.mkdir(exist_ok=True)
        pathlib.Path(str(output.feather)).parent.mkdir(exist_ok=True)
        reader=SptxtReader(filepath=str(input))
        df = reader.to_df(
            min_peaks=3,
            min_delta_ascore=20,
        )
        df.to_csv(str(output.csv), index=False)
        df.reset_index(drop=True).to_feather(str(output.feather))



rule add_rt_to_sptxt_csv:
    input:
        sptxt_df="exp_aggregated_sptxt_feather/{experiment}.mokapot.sptxt.feather",
        irt_df="rt_csv/{experiment}.irt.csv",
    output:
        csv="exp_aggregated_rt_sptxt_csv/{experiment}.mokapot.irt.sptxt.csv",
        feather="exp_aggregated_rt_sptxt_feather/{experiment}.mokapot.irt.sptxt.feather",
    run:
        sptxt_df = pd.read_feather(input.sptxt_df)

        try:
            rt_df = pd.read_csv(input.irt_df)
            del sptxt_df["iRT"]
            rt_df["ModSequences"] = [canonicalize_seq(x) for x in rt_df["StripPeptide"]]
            rt_df = rt_df[["ModSequences", "mean"]].rename(columns={"mean": "iRT"})
            sptxt_df = sptxt_df.merge(rt_df, how="left", on="ModSequences")
        except EmptyDataError as e:
            print("Not appending retention times.... because there are none")
            print(e)

        Path(str(output.csv)).parent.mkdir(exist_ok=True)
        Path(str(output.feather)).parent.mkdir(exist_ok=True)
        sptxt_df.to_csv(str(output.csv), index=False)
        sptxt_df.reset_index(drop=True).to_feather(str(output.feather))
