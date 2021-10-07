
def get_mokapot_command(target_dir="mokapot", addition=""):
    base_str = """
        set -e
        set -x

        trap 's=$?; echo >&2 "$0: Error on line "$LINENO": $BASH_COMMAND"; exit $s' ERR

        mkdir -p {target_dir}

        tot_line=$(cat {{input.pin_files}} | wc -l)
        extras=$(python -c "if 8000000 < $tot_line: print(f'--subset_max_train 5000000')")

        mokapot --verbosity 2 \
            --seed 2020 \
            --aggregate \
            --enzyme {{params.enzyme_regex}} \
            --decoy_prefix DECOY_ \
            --missed_cleavages 2 \
            --keep_decoys \
            --min_length 5 \
            --max_length 50 $extras \
            -d {target_dir} \
            -r {{wildcards.experiment}}{addition} \
            {{input.pin_files}} |& tee \
            {target_dir}/{{wildcards.experiment}}{addition}.mokapot.log

            # --subset_max_train 5000000 
            # This does not work right now on files with less than 5M psms....

        """
    out_str = base_str.format(target_dir=target_dir, addition=addition)

    return out_str
