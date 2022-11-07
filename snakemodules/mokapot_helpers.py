
def get_mokapot_command(target_dir="mokapot", addition=""):
    base_str = """
        set -e
        set -x

        trap 's=$?; echo >&2 "$0: Error on line "$LINENO": $BASH_COMMAND"; exit $s' ERR

        mkdir -p {target_dir}


        mokapot --verbosity 2 \
            --seed 2020 \
            --aggregate \
            --enzyme {{params.enzyme_regex}} \
            --decoy_prefix DECOY_ \
            --missed_cleavages 2 \
            --keep_decoys \
            --min_length 5 \
            --max_length 50 \
            --subset_max_train 5000000 \
            -d {target_dir} \
            -r {{wildcards.experiment}}{addition} \
            {{input.pin_files}} |& tee \
            {target_dir}/{{wildcards.experiment}}{addition}.mokapot.log


        """
    out_str = base_str.format(target_dir=target_dir, addition=addition)

    return out_str
