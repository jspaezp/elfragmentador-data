from elfragmentador.datamodules import _convert_tensor_columns_df


rule combine_and_split_train:
    input:
        [
            f"exp_aggregated_rt_sptxt_csv/{experiment}.mokapot.irt.sptxt.csv"
            for experiment in UNIQ_EXP
        ],
    output:
        train_csv="aggregated_rt_sptxt_csv/train.mokapot.irt.sptxt.csv",
        test_csv="aggregated_rt_sptxt_csv/test.mokapot.irt.sptxt.csv",
        val_csv="aggregated_rt_sptxt_csv/val.mokapot.irt.sptxt.csv",
        train_feather="aggregated_rt_sptxt_feather/train.mokapot.irt.sptxt.feather",
        test_feather="aggregated_rt_sptxt_feather/test.mokapot.irt.sptxt.feather",
        val_feather="aggregated_rt_sptxt_feather/val.mokapot.irt.sptxt.feather",
    run:
        random.seed(42)
        print("Reading Data")
        df = pd.concat(
            [
                _convert_tensor_columns_df(pd.read_csv(str(x)), verbose=False)
                for x in input
            ]
        )

        print("Generating random strings")
        seq_numbers = {k: random.random() for k in set(df["Sequences"])}

        df_train = df[[seq_numbers[x] <= 0.4 for x in df["Sequences"]]]
        df_train = df_train.reset_index(drop=True)
        print(f"Train dataset generated of length {len(df_train)}, {len(df)/len(df_train)}%")

        df_test = df[
            [seq_numbers[x] > 0.4 and seq_numbers[x] <= 0.7 for x in df["Sequences"]]
        ].reset_index(drop=True)
        print(f"Test dataset generated of length {len(df_test)}, {len(df)/len(df_test)}%")

        df_val = df[
            [seq_numbers[x] > 0.7 and seq_numbers[x] <= 1.0 for x in df["Sequences"]]
        ].reset_index(drop=True)
        print(f"Val dataset generated of length {len(df_val)}, {len(df)/len(df_val)}%")

        Path(str(output[0])).parent.mkdir(exist_ok=True)

        print("Writting Files")
        df_train.to_csv(str(output.train_csv), index=False)
        df_train.to_feather(str(output.train_feather))

        df_test.to_csv(str(output.test_csv), index=False)
        df_test.to_feather(str(output.test_feather))

        df_val.to_csv(str(output.val_csv), index=False)
        df_val.to_feather(str(output.val_feather))


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
        for inp, outp in zip(input, output):
            command = f"gzip -c {str(inp)} > {str(outp)}"
            print(command)
            shell(command)
