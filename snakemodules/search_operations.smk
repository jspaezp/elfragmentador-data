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


rule comet_decoy_plot:
    input:
        forward_pepxml="comet/{sample}.pep.xml",
        forward_pin="comet/{sample}.pin",
    output:
        xcorr_plot="comet/{sample}.xcorr.png",
        expect_score="comet/{sample}.lnexpect.png",
    run:
        # Plotting the xcorr distribution of decoys and targets
        df = pd.read_csv(
            f"comet/{wildcards.sample}.pin",
            index_col=False,
            error_bad_lines=False,
            sep="\t",
            usecols=["Xcorr", "lnExpect", "Label"],
            lineterminator="\n",
        )

        print(df)

        tps = []
        fps = []

        tps.extend([s for s, l in zip(df["Xcorr"], df["Label"]) if l > 0])
        fps.extend([s for s, l in zip(df["Xcorr"], df["Label"]) if l < 0])

        plt.hist(tps, bins=20, alpha=0.8, color="cyan")
        plt.hist(fps, bins=20, alpha=0.8, color="magenta")
        plt.savefig(
            f"comet/{wildcards.sample}.xcorr.png",
        )

        tps = []
        fps = []

        tps.extend([s for s, l in zip(df["lnExpect"], df["Label"]) if l > 0])
        fps.extend([s for s, l in zip(df["lnExpect"], df["Label"]) if l < 0])

        plt.hist(tps, bins=20, alpha=0.8, color="cyan")
        plt.hist(fps, bins=20, alpha=0.8, color="magenta")
        plt.savefig(
            f"comet/{wildcards.sample}.lnexpect.png",
        )
