from pyteomics import fasta


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
        "fasta/{file}.fasta",
    output:
        "fasta/{file}.decoy.fasta",
    run:
        # Note that it writes a concatenated decoy+target database, not only decoys
        fasta.write_decoy_db(str(input), str(output))


rule biognosys_irt_fasta:
    output:
        "fasta/irtfusion.fasta",
    shell:
        """
        mkdir -p fasta
        wget \
            https://biognosys.com/media.ashx/irtfusion.fasta \
            -O fasta/irtfusion.fasta
        """


rule arabidopsis_fasta:
    output:
        "fasta/arabidopsis.fasta"
    shell:
        """
        mkdir -p fasta
        wget \
            https://www.uniprot.org/uniprot/?query=proteome:UP000006548%20reviewed:yes\&format\=fasta \
            -O fasta/arabidopsis.fasta
        """


rule human_fasta:
    output:
        "fasta/human.fasta",
    shell:
        """
        mkdir -p fasta
        wget \
            https://www.uniprot.org/uniprot/\?query\=proteome:UP000005640%20reviewed:yes\&format\=fasta \
            -O fasta/human.fasta
        """


rule yeast_fasta:
    output:
        "fasta/yeast.fasta",
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
        "fasta/yeast.fasta",
    output:
        "fasta/human_yeast.fasta",
    shell:
        """
        cat fasta/human.fasta fasta/yeast.fasta > fasta/human_yeast.fasta
        """


rule contam_fasta:
    input:
        "fasta/crap.fasta",
        "fasta/{basename}.fasta",
    output:
        "fasta/{basename}_contam.fasta",
    shell:
        """
        cat fasta/{wildcards.basename}.fasta fasta/crap.fasta > fasta/{wildcards.basename}_contam.fasta
        """


rule added_irt_fasta:
    input:
        "fasta/irtfusion.fasta",
        "fasta/{basename}.fasta",
    output:
        "fasta/{basename}_irt.fasta",
    shell:
        """
        cat fasta/{wildcards.basename}.fasta fasta/irtfusion.fasta > fasta/{wildcards.basename}_irt.fasta
        """


rule proteometools_fasta:
    input:
        "fasta/RT_QC.fasta",
        "fasta/crap.fasta",
        "fasta/proteometools/Packet_{basename}.fasta",
    output:
        "fasta/{basename}_contam.fasta",
    shell:
        """
        cat fasta/proteometools/Packet_{wildcards.basename}.fasta \
            fasta/RT_QC.fasta \
            fasta/crap.fasta > fasta/{wildcards.basename}_contam.fasta
        """


def get_exp_fasta(wildcards):
    print(wildcards)
    return exp_to_fasta[wildcards.experiment]


rule experiment_fasta:
    input:
        get_exp_fasta,
    output:
        "experiment_fasta/{experiment}.fasta",
    shell:
        """
        set -x
        set -e

        mkdir -p experiment_fasta
        cat {input} > experiment_fasta/{wildcards.experiment}.fasta
        """
