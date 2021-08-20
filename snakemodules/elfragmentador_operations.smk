
include: "./mokapot_operations.smk"
include: "./search_operations.smk"


import pathlib


rule elfragmentador_pin:
    input:
        "comet/{sample}.pin",
    output:
        pin="ef_comet_pin/{sample}.elfragmentador.pin",
        log="ef_comet_pin/{sample}.elfragmentador.pin.log",
    shell:
        """
        set -e
        set -x
        trap 's=$?; echo >&2 "$0: Error on line "$LINENO": $BASH_COMMAND"; exit $s' ERR

        mkdir -p ef_comet_pin

        elfragmentador_append_pin \
            --pin {input} \
            --out {output.pin} |& tee {output.log} 

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
        decoy_psms_ef="mokapot/{experiment}.mokapot.decoy.psms.txt",
        decoy_psms="mokapot/{experiment}.mokapot.psms.txt",
        psms_ef="mokapot/{experiment}.mokapot.psms.txt",
        psms="mokapot/{experiment}.mokapot.decoy.psms.txt",
        model_importance="mokapot/{experiment}.mokapot.weights.csv",
        model_importance_ef="ef_mokapot/{experiment}.elfragmentador.mokapot.weights.csv",
    output:
        html="ef_reports/{experiment}.report.html",
        output_swapped_psms="ef_reports/{experiment}.swapped.csv",
        output_top_swapped_psms="ef_reports/{experiment}.swapped.top.csv",
    run:
        cmd = """
                set -x
                set -e

                mkdir -p ef_reports
                R -e \
                    "rmarkdown::render(\
                        'templates/plot_spectrum_comp_ef.Rmd',\
                        params = list( \
                            elfragmentador_pin='{input.elfragmentador_pin}', \
                            decoy_psms_ef='{input.decoy_psms_ef}', \
                            decoy_psms='{input.decoy_psms}', \
                            psms_ef='{input.psms_ef}', \
                            psms='{input.psms}', \
                            model_importance='{input.model_importance}', \
                            model_importance_ef='{input.model_importance_ef}', \
                            output_swapped_psms='{output.output_swapped_psms}', \
                            output_top_swapped_psms='{output.output_top_swapped_psms}' \
                        ), \
                        output_file = '{output.html}', \
                        clean = FALSE, \
                        knit_root_dir = getwd(), \
                        output_dir = 'ef_reports')" 
                """
        print(cmd)
        shell(cmd)
