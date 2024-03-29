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


rule generic_uniprot_fasta:
    output:
        "fasta/Uniprot/{accession,UP[0-9]+}.fasta",
    shell:
        """
        mkdir -p fasta/Uniprot
        wget \
            https://www.uniprot.org/uniprot/?query=proteome:{wildcards.accession}%20reviewed:yes\&format\=fasta \
            -O {output}
        """


rule arabidopsis_fasta:
    output:
        "fasta/arabidopsis.fasta",
    shell:
        """
        mkdir -p fasta
        wget \
            'https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28proteome%3AUP000006548%29%20AND%20%28reviewed%3Atrue%29' \
            -O fasta/arabidopsis.fasta
        """


rule human_fasta:
    output:
        "fasta/human.fasta",
    shell:
        """
        mkdir -p fasta
        wget \
            'https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28proteome%3AUP000005640%29%20AND%20%28reviewed%3Atrue%29' \
            -O fasta/human.fasta
        """


rule yeast_fasta:
    output:
        "fasta/yeast.fasta",
    shell:
        """
        mkdir -p fasta
        wget \
            'https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28proteome%3AUP000002311%29%20AND%20%28reviewed%3Atrue%29' \
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


rule proteometools2_fasta:
    input:
        "fasta/proteometools2/INDIVIDUAL_PEPTIDE_PROCAL.fasta",
        "fasta/proteometools2/{basename}.fasta",
    output:
        "fasta/{basename}_indiv_irt.fasta",
    shell:
        # No crap is added because this will be searched without cleavages,
        # I would need to add every peptide as as fasta entry ...
        """
        cat {input} > {output}
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
