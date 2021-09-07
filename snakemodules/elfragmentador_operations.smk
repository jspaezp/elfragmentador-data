
include: "./mokapot_operations.smk"
include: "./search_operations.smk"


import pathlib


rule elfragmentador_pin:
    input:
        "comet/{sample}.pin",
    output:
        pin="ef_comet_pin/{sample}.elfragmentador.pin",
        log="ef_comet_pin/{sample}.elfragmentador.pin.log",
    params:
        checkpoint=f"{CHECKPOINT}",
    threads: 8
    shell:
        """
        set -e
        set -x
        trap 's=$?; echo >&2 "$0: Error on line "$LINENO": $BASH_COMMAND"; exit $s' ERR

        mkdir -p ef_comet_pin

        poetry run elfragmentador_append_pin \
            --pin {input} \
            --out {output.pin} \
            --threads {threads} \
            --model_checkpoint {params.checkpoint} |& tee {output.log} 

        """


rule elfragmentador_mokapot:
    input:
        pin_files=get_mokapot_ins("ef_comet_pin/", ".elfragmentador.pin"),
    output:
        "ef_mokapot/{experiment}.elfragmentador.mokapot.psms.txt",
        "ef_mokapot/{experiment}.elfragmentador.mokapot.peptides.txt",
        "ef_mokapot/{experiment}.elfragmentador.mokapot.decoy.psms.txt",
        "ef_mokapot/{experiment}.elfragmentador.mokapot.decoy.peptides.txt",
        "ef_mokapot/{experiment}.elfragmentador.mokapot.log",
    params:
        enzyme_regex=get_enzyme_regex,
    shell:
        get_mokapot_command(target_dir="ef_mokapot", addition=".elfragmentador")


rule ef_mokapot_feature_weights:
    input:
        "ef_mokapot/{experiment}.elfragmentador.mokapot.log",
    output:
        "ef_mokapot/{experiment}.elfragmentador.mokapot.weights.csv",
    run:
        weights = parse_mokapot_log(str(input))
        weights.to_csv(str(output), index=False)


rule generate_report:
    input:
        elfragmentador_pin=get_mokapot_ins("ef_comet_pin/", ".elfragmentador.pin"),
        model_importance_ef="ef_mokapot/{experiment}.elfragmentador.mokapot.weights.csv",
        decoy_psms_ef="ef_mokapot/{experiment}.elfragmentador.mokapot.decoy.psms.txt",
        psms_ef="ef_mokapot/{experiment}.elfragmentador.mokapot.psms.txt",
        model_importance="mokapot/{experiment}.mokapot.weights.csv",
        decoy_psms="mokapot/{experiment}.mokapot.decoy.psms.txt",
        psms="mokapot/{experiment}.mokapot.psms.txt",
    output:
        html="ef_reports/{experiment}.report.html",
        output_swapped_psms="ef_reports/{experiment}.swapped.csv",
        output_top_swapped_psms="ef_reports/{experiment}.swapped.top.csv",
    params:
        checkpoint=f"{CHECKPOINT}"
    run:
        cmd = [
                " set -x  ;                                                               ",
                " set -e  ;                                                               ",
                " mkdir -p ef_reports ;                                                   ",
                " R -e \"rmarkdown::render(                                               ",
                "         'templates/plot_spectrum_comp_ef.Rmd',                          ",
                "         params = list(                                                  ",
                "             elfragmentador_pin='{input.elfragmentador_pin}',            ",
                "             decoy_psms_ef='{input.decoy_psms_ef}',                      ",
                "             decoy_psms='{input.decoy_psms}',                            ",
                "             psms_ef='{input.psms_ef}',                                  ",
                "             psms='{input.psms}',                                        ",
                "             model_importance='{input.model_importance}',                ",
                "             model_importance_ef='{input.model_importance_ef}',          ",
                "             output_swapped_psms='{output.output_swapped_psms}',         ",
                "             output_top_swapped_psms='{output.output_top_swapped_psms}'  ",
                "         ),                                                              ",
                "         output_file = '{output.html}',                                  ",
                "         clean = FALSE,                                                  ",
                "         knit_root_dir = getwd(), intermediates_dir=tempdir(),           ",
                "         output_dir = 'ef_reports')\"                                    ",
            ]

        cmd = "".join([x.strip() for x in cmd])
        print(cmd)
        shell(cmd)

        cmd = (
            "poetry run python ./scripts/plot_top_n.py "
            "--checkpoint {params.checkpoint} "
            "--csv {output.output_top_swapped_psms} "
            "--search_path . "
            "--prefix ef_reports/{wildcards.experiment}"
        )
        shell(cmd)


rule generate_roc_curves:
    input:
        psms_ef="ef_mokapot/{experiment}.elfragmentador.mokapot.psms.txt",
        psms="mokapot/{experiment}.mokapot.psms.txt",
    output:
        html="ef_reports/{experiment}.roc_curves.html",
    run:
        cmd = [
            " set -x ;                                         ",
            " set -e ;                                         ",
            " mkdir -p ef_reports ;                            ",
            " R -e  \"rmarkdown::render(                       ",
            "         'templates/plot_roc_curves.Rmd',         ",
            "         params = list(                           ",
            "             psms_ef='{input.psms_ef}',           ",
            "             psms='{input.psms}'                  ",
            "         ),                                       ",
            "         output_file = '{output.html}',           ",
            "         clean = FALSE,                           ",
            "         knit_root_dir = getwd(),                 ",
            "         output_dir = 'ef_reports')\"             ",
        ]
        cmd = "".join([x.strip() for x in cmd])
        print(cmd)
        shell(cmd)


rule evaluation_on_psms:
    input:
        mokapot_psms="mokapot/{experiment}.mokapot.psms.txt",
    output:
        out_csv="ef_evaluation/{experiment}.csv",
        log="ef_evaluation/{experiment}.log"
    params:
        checkpoint=f"{CHECKPOINT}",
    threads: 8
    run:
        cmd = [
            "mkdir -p ef_evaluation ; ",
            " elfragmentador_evaluate --input {input.mokapot_psms} ",
            " --screen_nce='-10,-8,-6,-4,-2,0,2,4,6,8,10' --batch_size 50 --max_spec 100000 ",
            " --out_csv {output.out_csv} ",
            " --model_checkpoint {params.checkpoint} ",
            " --threads {threads} | tee {output.log}",
        ]

        cmd = "".join(cmd)
        print(cmd)
        shell(cmd)

    
