import io
from pathlib import Path

import numpy as np
import pandas as pd
from loguru import logger as lg_logger
from pyteomics.mzml import PreIndexedMzML
from tqdm.auto import tqdm


def filter_mokapot_psm(psm_input: str, peptide_input: str):
    """
    psm_input = "mokapot/MannLabPhospho.mokapot.psms.txt"
    peptide_input = "mokapot/MannLabPhospho.mokapot.peptides.txt"
    filter_mokapot_psm(psm_input, peptide_input)
    """

    df = pd.read_table(str(psm_input))
    pep_df = pd.read_table(str(peptide_input))
    print(f"Filtering peptides for {str(psm_input)} {str(peptide_input)}")
    print(f"Pep df length {len(pep_df)}")
    pep_df = pep_df[(pep_df["mokapot q-value"] < 0.01) & (pep_df["mokapot PEP"] < 0.01)]
    passed_peptides = set(pep_df["Peptide"])
    print(f"Unique Peptides Passed {len(passed_peptides)}")

    if len(passed_peptides) == 0:
        raise ValueError(
            f"Numbber of peptides that passed is 0 for {psm_input} {peptide_input}"
        )

    print(f"Number of PSMs {len(df)}")
    df = df[(df["mokapot q-value"] < 0.01) & (df["mokapot PEP"] < 0.01)]

    df = df[[x in passed_peptides for x in df["Peptide"]]]
    print(f"Number of PSMs {len(df)} after peptide filtering")
    return df


def psm_df_to_tsv(df: pd.DataFrame, output: str):
    col_neworder = [
        "SpecId",
        "ScanNr",
        "Peptide",
        "mokapot PEP",
        "mokapot score",
        "Proteins",
    ]
    df["Peptide"] = [x.replace("[", "[+") for x in df["Peptide"]]
    df["Charge"] = ["_".join(line.split("_")[-2]) for line in df["SpecId"]]
    df["Peptide"] = [f"{x}/{y}" for x, y in zip(df["Peptide"], df["Charge"])]
    df["SpecId"] = ["_".join(line.split("_")[:-3]) for line in df["SpecId"]]
    df["mokapot PEP"] = 1 - df["mokapot PEP"]

    df = df[col_neworder]
    df.to_csv(str(output), sep="\t", index=False)
    print(df)


def split_mokapot_spectrast_in(
    spectrast_in: str, spec_metadata_list: list, experiment: str
) -> None:
    """
    spectrast_in = "spectrast_in/ProteomeToolsSP.spectrast.mokapot.psms.tsv"
    spec_metadata_list = [
        "raw_scan_metadata/01709a_GA9-TUM_second_pool_1_01_01-2xIT_2xHCD-1h-R1.csv",
        "raw_scan_metadata/01748a_BF2-TUM_second_pool_48_01_01-3xHCD-1h-R1.csv",
        "raw_scan_metadata/01748a_BG2-TUM_second_pool_49_01_01-2xIT_2xHCD-1h-R1.csv",
        "raw_scan_metadata/01748a_BG2-TUM_second_pool_49_01_01-3xHCD-1h-R1.csv"
    ]
    split_mokapot_spectrast_in(spectrast_in, spec_metadata_list, experiment = "testing")
    """

    outdir = f"split_spectrast_in/{experiment}"
    required_cols = [
        "SpecId",
        "ScanNr",
        "Peptide",
        "mokapot PEP",
        "mokapot score",
        "Proteins",
    ]
    # TODO add generation of ssl file for bibliospec
    # https://skyline.ms/wiki/home/software/BiblioSpec/page.view?name=BiblioSpec%20input%20and%20output%20file%20formats
    ssl_columns = ["file", "scan", "charge", "sequence", "score", "retention-time"]
    column_aliases = {
        "SpecId": "file",
        "ScanNr": "scan",
        "Peptide": "sequence",
        "score": "mokapot PEP",
        "RetentionTime": "retention-time",
    }

    print("Reading Main DF")
    main_df = pd.read_csv(str(spectrast_in), sep="\t")
    main_df = main_df.set_index(["SpecId", "ScanNr"])
    print(f"Columns: {list(main_df)}")

    print("Reading Secondary DFs")
    for sec_spec in tqdm(spec_metadata_list):
        sec_df = pd.read_csv(str(sec_spec))
        sec_df = sec_df.set_index(["SpecId", "ScanNr"])
        main_df = main_df.join(sec_df, lsuffix="_extra")

        extra_cols = [x for x in list(main_df) if x.endswith("_extra")]
        for ec in extra_cols:
            curr_vals = main_df[ec.replace("_extra", "")]
            extra_vals = main_df[ec]
            missing_vals = np.isnan(curr_vals)
            print(f"{missing_vals.sum()} Missing Values")

            # This makes it so it tries to replace missing values with the ones in
            # extra_vals (suffixed column), otherwise keeps the old one.
            replace_val = np.where(missing_vals, extra_vals, curr_vals)
            main_df[ec.replace("_extra", "")] = replace_val
            del main_df[ec]

    print(f"Columns: {list(main_df)}")
    print("Getting Unique collision Energies")
    unique_ce = np.unique(main_df["CollisionEnergy"])
    unique_ce = unique_ce[~np.isnan(unique_ce)]
    print(f"uniques: {unique_ce}")

    for ce in tqdm(unique_ce):
        out_name = f"{outdir}/{str(experiment)}.{str(ce).replace('.', '_')}.spectrast.mokapot.psms.tsv"
        tmp_df = main_df[ce == main_df["CollisionEnergy"]].reset_index()[required_cols]
        print(f"Saving {out_name}")
        # TODO add here removal of non-replicated spectra
        # TODO also consider if you could split when the file is too large
        tmp_df.to_csv(out_name, sep="\t", index=False)

    out_name = f"{outdir}/{str(experiment)}.UNKNOWN.spectrast.mokapot.psms.tsv"
    print(f"Saving {out_name}")
    tmp_df = main_df[np.isnan(main_df["CollisionEnergy"])].reset_index()[required_cols]
    tmp_df.to_csv(out_name, sep="\t", index=False)


def add_sptxt_ce_info_meta():
    pass


def add_spectrast_ce_info(mzml_directory, in_sptxt, out_sptxt):
    """
    add_ce_info("raw/", "./mokapot_spectrast/MannLabPhospho.mokapot.psms.sptxt", "foo.txt")
    """
    id_template = "controllerType=0 controllerNumber=1 scan={scan_number}"
    add_template = "RetentionTime={rt} CollisionEnergy={ce}"
    mzml_dict = {}
    with open(str(in_sptxt), "r") as infile:
        with open(str(out_sptxt), "w") as outfile:
            for i, line in enumerate(tqdm(infile)):
                line = line.strip()
                if line.startswith("Comment"):
                    l_list = line.strip().split(" ")

                    extra_left, extra_right = l_list[:-3], l_list[-3:]
                    spectrum = (
                        [x for x in extra_left if x.startswith("RawSpectrum")][0]
                        .split("=")[1]
                        .split(".")
                    )

                    raw_file = ".".join(spectrum[:-2])
                    spec_id = spectrum[-1]

                    if raw_file not in mzml_dict:
                        mzml_file = f"{mzml_directory}/{raw_file}.mzML"
                        mzml_dict[raw_file] = PreIndexedMzML(mzml_file)

                    spec = mzml_dict[raw_file][id_template.format(scan_number=spec_id)]
                    assert spec["ms level"] != 1

                    ce = spec["precursorList"]["precursor"]
                    ce = ce[0]["activation"]["collision energy"]
                    rt = spec["scanList"]["scan"][0]["scan start time"]
                    addition = add_template.format(rt=rt, ce=ce)
                    line = " ".join(extra_left + [addition] + extra_right)

                outfile.write(line + "\n")


def parse_mokapot_log(path):
    FEATURE_START_CONTENT = "Normalized feature weights in the learned model"
    # FEATURE_START_CONTENT = "Feature   Weight"
    FEATURE_END_CONTENT = "Done training."

    in_feature_context = False
    fold = 0
    feature_list = []
    df_list = []

    with open(path, "r") as f:
        for line in f:
            if not in_feature_context:
                if FEATURE_START_CONTENT in line:
                    in_feature_context = True
                    fold += 1
            else:
                if FEATURE_END_CONTENT in line:
                    in_feature_context = False
                    iostring = io.StringIO("\n".join(feature_list))
                    tmp_df = pd.read_csv(iostring, delim_whitespace=True)
                    tmp_df["Fold"] = fold
                    df_list.append(tmp_df)
                    feature_list = []
                else:
                    feature_list.append(line.replace("[INFO]", "").strip())

    return pd.concat(df_list)


def split_mokapot_psms_file(psms_file, spec_metadata_list, out_dir):
    """

    Example section of metadata file
        ScanNr,CollisionEnergy,RetentionTime,SpecId
        1032,25.0,9.2224227088,raw/01625b_GA1-TUM_first_pool_1_01_01-3xHCD-1h-R1
        1033,30.0,9.223171026667,raw/01625b_GA1-TUM_first_pool_1_01_01-3xHCD-1h-R1
        1034,35.0,9.224075944267,raw/01625b_GA1-TUM_first_pool_1_01_01-3xHCD-1h-R1

    Example section of psms file
        SpecId  Label   ScanNr  ExpMass CalcMass        Peptide mokapot score   mokapot
        raw/03210a_BF4-TUM_aspn_16_01_01-3xHCD-1h-R1_3465_4_2   False   3465    1825.04
        raw/03210a_BE5-TUM_aspn_5_01_01-3xHCD-1h-R1_12584_2_1   False   12584   1327.60
        raw/03210a_BF4-TUM_aspn_16_01_01-3xHCD-1h-R1_3055_4_2   False   3055    1825.04
        raw/03210a_BE10-TUM_aspn_10_01_01-3xHCD-1h-R1_11237_2_1 False   11237   1327.
    """
    out_name_template = Path(psms_file).stem + ".NCE{}.txt"

    # read both files
    lg_logger.info(f"Reading metadata file {spec_metadata_list}")
    metadata_df = pd.concat([
        pd.read_csv(meta_file, sep="\t", dtype={"ScanNr": np.int32})
        for meta_file in spec_metadata_list
    ])
    lg_logger.info(f"Reading psms file {psms_file}")
    psms_df = pd.read_csv(
        psms_file,
        sep="\t",
        dtype={
            "SpecId": str,
            "Label": bool,
            "ScanNr": np.int32,
            "ExpMass": np.float64,
            "CalcMass": np.float64,
            "Peptide": str,
            "mokapot score": np.float32,
            "mokapot q-value": np.float32,
            "mokapot PEP": np.float32,
            "Proteins": str,
        },
    )

    # Process psms file to get the spec id equivalent from the spec id
    lg_logger.info("Processing psms file")
    psms_df["RawFile"] = [SPEC_ID_REGEX.match(SpecId).groups()[0] for SpecId in psms_df.SpecId]

    # Left join the two files
    metadata_df.rename(columns = {'SpecId':'RawFile'}, inplace = True)
    lg_logger.info("Joining metadata and psms files")
    joined_df = psms_df.merge(metadata_df, on="RawFile", how="left")
    
    # Write the output
    lg_logger.info("Writing output")
    for nce in joined_df.CollisionEnergy.unique():
        out_name = out_name_template.format(nce)
        out_name = Path(out_dir) / out_name
        lg_logger.info(f"Writing output to {out_name}")
        joined_df[joined_df.CollisionEnergy == nce].to_csv(
            out_name, sep="\t", index=False, header=True
        )

if __name__ == "__main__":
    # spectrast_in = "spectrast_in/ProteomeToolsSP.spectrast.mokapot.psms.tsv"
    # spec_metadata_list = [
    #     "raw_scan_metadata/01709a_GA9-TUM_second_pool_1_01_01-2xIT_2xHCD-1h-R1.csv",
    #     "raw_scan_metadata/01748a_BF2-TUM_second_pool_48_01_01-3xHCD-1h-R1.csv",
    #     "raw_scan_metadata/01748a_BG2-TUM_second_pool_49_01_01-2xIT_2xHCD-1h-R1.csv",
    #     "raw_scan_metadata/01748a_BG2-TUM_second_pool_49_01_01-3xHCD-1h-R1.csv",
    # ]
    # experiment = "testing"
    # split_mokapot_spectrast_in(spectrast_in, spec_metadata_list, experiment)

    print(parse_mokapot_log("./mokapot/PEPTIDOME_EF2.mokapot.log"))
