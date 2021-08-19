import json
import elfragmentador.constants as CONSTANTS
from elfragmentador.math_utils import polyfit, apply_polyfit

def get_experiment_metadata_files(wildcards):
    outs = expand(
        "raw_scan_metadata/{sample}.csv", sample=get_samples(wildcards.experiment)
    )

    return outs

def _read_spectrast_in(spectrast_in, only_irt=False):
    df = pd.read_csv(spectrast_in, sep="\t", usecols=["SpecId", "ScanNr", "Peptide"])

    # R.PCYCSSGCGSSCCQSSCCK.S/2 to 'PCYCSSGCGSSCCQSSCCK'
    df["StripPeptide"] = [x[2:-4] for x in df["Peptide"]]

    if only_irt:
        df = df[[x in CONSTANTS.IRT_PEPTIDES for x in df["StripPeptide"]]].reset_index(
            drop=True
        )
        df["iRT"] = [CONSTANTS.IRT_PEPTIDES[x]["irt"] for x in df["StripPeptide"]]

    df["File"] = [pathlib.Path(x).stem for x in df["SpecId"]]
    return df


def calculate_irt_coefficients(spectrast_in, metadata_files):
    """
    Calculates the coefficients of a linear regression that fits
    iRT ~ RT using irt peptides as reference
    """

    df = _read_spectrast_in(spectrast_in=spectrast_in, only_irt=True)

    coefficients_dict = {}
    for tmp_file_name, tmp_df in df.groupby("File"):
        if len(tmp_df) < 5:
            continue
        curr_metadata_file_name = [x for x in metadata_files if tmp_file_name in x][0]
        tmp_metadata_df = pd.read_csv(curr_metadata_file_name)
        scan_times = {
            k: v
            for k, v in zip(tmp_metadata_df["ScanNr"], tmp_metadata_df["RetentionTime"])
        }
        tmp_df["RT"] = [scan_times[x] for x in tmp_df["ScanNr"]]
        fit_coefficients = polyfit(x=tmp_df["RT"], y=tmp_df["iRT"])
        coefficients_dict[tmp_file_name] = fit_coefficients

    return coefficients_dict


def calculate_irt_spectrast_in(spectrast_in, coefficients_dict, metadata_files):
    """Calculates the iRT of the spectra in a spectrast in, product from mokapot"""

    df = _read_spectrast_in(spectrast_in=spectrast_in, only_irt=False)

    outs = []

    for tmp_file_name, tmp_df in df.groupby("File"):
        curr_metadata_file_name = [x for x in metadata_files if tmp_file_name in x][0]
        tmp_metadata_df = pd.read_csv(curr_metadata_file_name)
        scan_times = {
            k: v
            for k, v in zip(tmp_metadata_df["ScanNr"], tmp_metadata_df["RetentionTime"])
        }
        tmp_df["RT"] = [scan_times[x] for x in tmp_df["ScanNr"]]
        poly = coefficients_dict.get(tmp_file_name, None)
        if poly is None:
            continue

        poly = poly["polynomial"]
        tmp_df["iRT"] = apply_polyfit(tmp_df["RT"], polynomial=poly)
        outs.append(tmp_df)

    concat_out = pd.concat(outs).reset_index(drop=True)

    return concat_out


rule calculate_irt_coefficients:
    input:
        spectrast_in_tsv="spectrast_in/{experiment}.spectrast.mokapot.psms.tsv",
        metadata_files=get_experiment_metadata_files,
    output:
        "irt_coefficients/{experiment}.coefficients.json",
    run:
        coefficients_dict = calculate_irt_coefficients(
            spectrast_in=input.spectrast_in_tsv, metadata_files=input.metadata_files
        )

        with open(str(output), "w") as f:
            json.dump(coefficients_dict, fp=f, indent=2)


rule generate_irt_csv:
    input:
        spectrast_in_tsv="spectrast_in/{experiment}.spectrast.mokapot.psms.tsv",
        metadata_files=get_experiment_metadata_files,
        coefficients_json="irt_coefficients/{experiment}.coefficients.json",
    output:
        "rt_csv/{experiment}.irt.csv",
    run:
        with open(input.coefficients_json, "r") as f:
            coefficients_dict = json.load(f)

        try:
            out_df = calculate_irt_spectrast_in(
                spectrast_in=input.spectrast_in_tsv,
                coefficients_dict=coefficients_dict,
                metadata_files=input.metadata_files,
            )

            out_df = out_df.groupby(["StripPeptide"])["iRT"]
            out_df = out_df.aggregate(["mean", "min", "max", "count"]).reset_index()

            out_df.to_csv(str(output), index=False)
        except ValueError as e:
            print(f"Handling value error {e}, creating empty file")
            pathlib.Path(str(output)).touch()
