import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt


include: "./env_setup.smk"


rule comet_phospho_params:
    input:
        "comet_params/comet.params.high_high",
    output:
        "comet_params/comet_phospho.params.high_high",
    shell:
        """
        cat {input} | \
            sed -e "s/variable_mod02 = 0.0 X 0 3 -1 0 0/variable_mod02 = 79.966331 STY 0 3 -1 0 0/g" \
            | tee {output}
        """


rule comet_gg_params:
    input:
        "comet_params/comet.params.high_high",
    output:
        "comet_params/comet_gg.params.high_high",
    shell:
        """
        cat {input} | \
            sed -e "s/variable_mod02 = 0.0 X 0 3 -1 0 0/variable_mod02 = 114.042927 K 0 3 -1 0 0/g" \
            | tee {output}
        """


rule comet_proalanase_params:
    input:
        "comet_params/comet.params.high_high",
    output:
        "comet_params/comet.params.proalanase.high_high",
    shell:
        """
        cat {input} | \
            sed -e "s/^search_enzyme_number.*/search_enzyme_number = 10/g" | \
            sed -e "s/^10. Chymotrypsin.*/10. ProAlanase 1 PA -/g" \
            | tee {output}
        """


rule comet_peptidome_params:
    input:
        "comet_params/comet.params.high_high",
    output:
        "comet_params/comet.params.peptidome.high_high",
    shell:
        """
        cat {input} | \
            sed -e "s/^search_enzyme_number.*/search_enzyme_number = 0/g" | \
            sed -e "s/^digest_mass_range = 600.0 5000.0/digest_mass_range = 700.0 2700.0/g" \
            | tee {output}
        """


rule comet_nocleavage_params:
    input:
        "comet_params/comet.params.high_high",
    output:
        "comet_params/comet.nocleavage.params.high_high"
    shell:
        """
        cat {input} | \
            sed -e "s/^search_enzyme_number.*/search_enzyme_number = 10/g" | \
            sed -e "s/^10. Nothing.*/10. Nothing 1 Z -/g" \
            | tee {output}
        """


rule comet_semitriptic_params:
    input:
        "comet_params/comet.params.high_high",
    output:
        "comet_params/comet.semitriptic.params.high_high"
    shell:
        """
        cat {input} | \
            sed -r "s/num_enzyme_termini.*/num_enzyme_termini = 9 /g" \
            | tee {output}
        """

def get_fasta(wildcards):
    return samp_to_fasta[wildcards.sample]


def get_comet_params(wildcards):
    return samp_to_params[wildcards.sample]


rule comet_search:
    input:
        raw="raw/{sample}.mzML",
        fasta=get_fasta,
        comet_params=get_comet_params,
    output:
        # Enable if using 2 in the decoy search parameter
        # decoy_pepxml = "comet/{sample}.decoy.pep.xml", 
        forward_pepxml="comet/{sample}.pep.xml",
        forward_pin="comet/{sample}.pin",
    benchmark:
        "benchmarks/{sample}.comet_search.benchmark.txt"
    run:
        shell("mkdir -p comet")
        cmd = (
            f"{TPP_DOCKER} "
            f"comet -P{str(input.comet_params)} "
            f"-D{str(input.fasta)} "
            f"{str(input.raw)} "
        )

        print(cmd)
        shell(cmd)
        shell(f"cp raw/{wildcards.sample}.pep.xml ./comet/.")
        shell(f"cp raw/{wildcards.sample}.pin ./comet/.")


rule clip_comet_pin:
    input:
        forward_pin="comet/{sample}.pin",
    output:
        clipped_pin="comet/{sample}.clipped.pin",
    shell:
        """
	grep -oP "^(\S+(\s|$)+){{28}}" \
            {input} | \
            sed -e "s/\s+$//g" > {output}
        """

rule comet_decoy_plot:
    input:
        forward_pepxml="comet/{sample}.pep.xml",
        clipped_pin="comet/{sample}.clipped.pin",
    output:
        xcorr_plot="comet/{sample}.xcorr.png",
        expect_score="comet/{sample}.lnexpect.png",
    run:
        df = pd.read_csv(
            str(input.clipped_pin),
            index_col=False,
            error_bad_lines=False,
            sep="\t",
            usecols=["Xcorr", "lnExpect", "Label"],
            lineterminator="\n",
            low_memory=False,
        )

        print(df)

        decoys = df[df["Label"] < 0]
        targets = df[df["Label"] > 0]

        plt.hist(targets["Xcorr"], bins=20, alpha=0.8, color="cyan")
        plt.hist(decoys["Xcorr"], bins=20, alpha=0.8, color="magenta")
        plt.title("Xcorr Distribution (Targets are Cyan)")
        plt.savefig(
            f"comet/{wildcards.sample}.xcorr.png",
        )
        plt.close()

        plt.hist(targets["lnExpect"], bins=20, alpha=0.8, color="cyan")
        plt.hist(decoys["lnExpect"], bins=20, alpha=0.8, color="magenta")
        plt.title("lnExpect Distribution (Targets are Cyan)")
        plt.savefig(
            f"comet/{wildcards.sample}.lnexpect.png",
        )
        plt.close()
