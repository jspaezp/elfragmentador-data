
rule combine_and_split_train:
    input:
        [
            f"exp_aggregated_rt_sptxt_csv/{experiment}.mokapot.irt.sptxt.csv"
            for experiment in UNIQ_EXP
        ],
    output:
        train_df="aggregated_rt_sptxt_csv/train.mokapot.irt.sptxt.csv",
        test_df="aggregated_rt_sptxt_csv/test.mokapot.irt.sptxt.csv",
        val_df="aggregated_rt_sptxt_csv/val.mokapot.irt.sptxt.csv",
    run:
        random.seed(42)
        df = pd.concat([pd.read_csv(str(x)) for x in input])

        seq_numbers = {k: random.random() for k in set(df["Sequences"])}

        df_train = df[[seq_numbers[x] <= 0.4 for x in df["Sequences"]]]
        df_test = df[
            [seq_numbers[x] > 0.4 and seq_numbers[x] <= 0.7 for x in df["Sequences"]]
        ]
        df_val = df[
            [seq_numbers[x] > 0.7 and seq_numbers[x] <= 1.0 for x in df["Sequences"]]
        ]
        Path(str(output[0])).parent.mkdir(exist_ok=True)

        df_train.to_csv(str(output.train_df), index=False)
        df_test.to_csv(str(output.test_df), index=False)
        df_val.to_csv(str(output.val_df), index=False)


rule zip_and_timestamp:
    input:
        train_df="aggregated_rt_sptxt_csv/train.mokapot.irt.sptxt.csv",
        test_df="aggregated_rt_sptxt_csv/test.mokapot.irt.sptxt.csv",
        val_df="aggregated_rt_sptxt_csv/val.mokapot.irt.sptxt.csv",
    output:
        train_df="zip_aggregated_rt_sptxt_csv/{timestamp}.train.mokapot.irt.sptxt.csv.gz",
        test_df="zip_aggregated_rt_sptxt_csv/{timestamp}.test.mokapot.irt.sptxt.csv.gz",
        val_df="zip_aggregated_rt_sptxt_csv/{timestamp}.val.mokapot.irt.sptxt.csv.gz",
    run:
        Path(str(output[0])).parent.mkdir(exist_ok=True)
        for x, y in zip(input, output):
            command = f"gzip -c {str(y)} > {str(x)}"
            print(command)
            shell(command)
