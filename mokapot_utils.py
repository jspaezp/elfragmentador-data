
from pyteomics.mzml import PreIndexedMzML
import pandas as pd
import numpy as np
from tqdm.auto import tqdm


def filter_mokapot_psm(psm_input: str, peptide_input: str):
    """
    psm_input = "mokapot/MannLabPhospho.mokapot.psms.txt"
    peptide_input = "mokapot/MannLabPhospho.mokapot.peptides.txt"
    filter_mokapot_psm(psm_input, peptide_input)
    """

    df = pd.read_table(str(psm_input))
    pep_df = pd.read_table(str(peptide_input))
    print(f"Pep df length {len(pep_df)}")
    pep_df = pep_df[
        (pep_df['mokapot q-value'] < 0.001) &
        (pep_df['mokapot PEP'] < 0.01)
    ]
    passed_peptides = set(pep_df['Peptide'])

    print(f"Unique Peptides {len(passed_peptides)}")
    print(f"Number of PSMs {len(df)}")
    df = df[
        (df['mokapot q-value'] < 0.001) &
        (df['mokapot PEP'] < 0.01)
    ]

    df = df[[x in passed_peptides for x in df['Peptide']]]
    print(f"Number of PSMs {len(df)} after peptide filtering")
    return(df)


def psm_df_to_tsv(df: pd.DataFrame, output: str):
    # TODO change this so it uses column names
    index_order = [0, 2, 5, 9, 7, 6]
    df['Peptide'] = [x.replace("[", "[+") for x in df['Peptide']]
    df['Charge'] = ["_".join(line.split("_")[-2]) for line in df['SpecId']]
    df['Peptide'] = [f"{x}/{y}" for x, y in zip(df['Peptide'], df['Charge'])]
    df['SpecId'] = ["_".join(line.split("_")[:-3]) for line in df['SpecId']]
    df['mokapot PEP'] = 1-df['mokapot PEP']

    col_neworder = [list(df)[x] for x in index_order]
    df = df[col_neworder]
    df.to_csv(str(output), sep='\t', index=False)
    print(df)

def split_mokapot_spectrast_in(spectrast_in: str, spec_metadata_list: list, experiment: str) -> None:
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

    outdir = "split_spectrast_in"
    required_cols = ["SpecId",  "ScanNr",  "Peptide", "mokapot PEP", "mokapot score", "Proteins"]
    # TODO add generation of ssl file for bibliospec
    # https://skyline.ms/wiki/home/software/BiblioSpec/page.view?name=BiblioSpec%20input%20and%20output%20file%20formats
    ssl_columns = ["file", "scan", "charge", "sequence", "score", "retention-time"]
    column_aliases = {
        'SpecId':'file',
        'ScanNr':'scan',
        'Peptide':'sequence',
        'score':'mokapot PEP',
        'RetentionTime':'retention-time'}

    print("Reading Main DF")
    main_df = pd.read_csv(str(spectrast_in), sep = '\t')
    main_df = main_df.set_index(['SpecId', 'ScanNr'])

    print("Reading Secondary DFs")
    for sec_spec in tqdm(spec_metadata_list):
        sec_df = pd.read_csv(str(sec_spec))
        sec_df = sec_df.set_index(['SpecId', 'ScanNr'])
        main_df = main_df.join(sec_df, lsuffix = '_extra')

        extra_cols = [x for x in list(main_df) if x.endswith('_extra')]
        for ec in extra_cols:
            curr_vals = main_df[ec.replace('_extra', '')]
            extra_vals = main_df[ec]
            missing_vals = np.isnan(curr_vals)
            print(f"{missing_vals.sum()} Missing Values")

            # This makes it so it tries to replace missing values with the ones in
            # extra_vals (suffixed column), otherwise keeps the old one.
            replace_val = np.where(missing_vals, extra_vals, curr_vals)
            main_df[ec.replace('_extra', '')] = replace_val
            del main_df[ec]

    print("Getting Unique collision Energies")
    unique_ce = np.unique(main_df['CollisionEnergy'])
    unique_ce = unique_ce[~np.isnan(unique_ce)]
    print(f"uniques: {unique_ce}")

    for ce in tqdm(unique_ce):
        out_name = f"{outdir}/{str(experiment)}.{str(ce).replace('.', '_')}.spectrast.mokapot.psms.tsv"
        print(f"Saving {out_name}")
        tmp_df = main_df[ce == main_df['CollisionEnergy']].reset_index()[required_cols]
        # TODO add here removal of non-replicated spectra
        # TODO also consider if you could split when the file is too large
        tmp_df.to_csv(out_name, sep = "\t", index = False)

    out_name = f"{outdir}/{str(experiment)}.UNKNOWN.spectrast.mokapot.psms.tsv"
    print(f"Saving {out_name}")
    tmp_df = main_df[np.isnan(main_df['CollisionEnergy'])].reset_index()[required_cols]
    tmp_df.to_csv(out_name, sep = "\t", index = False)


def add_sptxt_ce_info_meta():
    pass

def add_spectrast_ce_info(mzml_directory, in_sptxt, out_sptxt):
    """
    add_ce_info("raw/", "./mokapot_spectrast/MannLabPhospho.mokapot.psms.sptxt", "foo.txt")
    """
    id_template = 'controllerType=0 controllerNumber=1 scan={scan_number}'
    add_template = 'RetentionTime={rt} CollisionEnergy={ce}'
    mzml_dict = {}
    with open(str(in_sptxt), 'r') as infile:
        with open(str(out_sptxt), 'w') as outfile:
            for i, line in enumerate(tqdm(infile)):
                line = line.strip()
                if line.startswith('Comment'):
                    l_list = line.strip().split(" ")

                    extra_left, extra_right = l_list[:-3], l_list[-3:]
                    spectrum = [
                        x for x in extra_left
                        if x.startswith('RawSpectrum')
                    ][0].split('=')[1].split('.')

                    raw_file = '.'.join(spectrum[:-2])
                    spec_id = spectrum[-1]
                    
                    if raw_file not in mzml_dict:
                        mzml_file = f"{mzml_directory}/{raw_file}.mzML"
                        mzml_dict[raw_file] = PreIndexedMzML(mzml_file)
                 
                    spec = mzml_dict[raw_file][id_template.format(scan_number=spec_id)]
                    assert spec['ms level'] != 1

                    ce = spec['precursorList']['precursor']
                    ce = ce[0]['activation']['collision energy']
                    rt = spec['scanList']['scan'][0]['scan start time']
                    addition = add_template.format(rt=rt, ce=ce)
                    line = " ".join(extra_left+[addition]+extra_right)

                outfile.write(line + "\n")

if __name__ == "__main__":
    spectrast_in = "spectrast_in/ProteomeToolsSP.spectrast.mokapot.psms.tsv"
    spec_metadata_list = [
        "raw_scan_metadata/01709a_GA9-TUM_second_pool_1_01_01-2xIT_2xHCD-1h-R1.csv",
        "raw_scan_metadata/01748a_BF2-TUM_second_pool_48_01_01-3xHCD-1h-R1.csv",
        "raw_scan_metadata/01748a_BG2-TUM_second_pool_49_01_01-2xIT_2xHCD-1h-R1.csv",
        "raw_scan_metadata/01748a_BG2-TUM_second_pool_49_01_01-3xHCD-1h-R1.csv"
    ]
    experiment = "testing"
    split_mokapot_spectrast_in(spectrast_in, spec_metadata_list, experiment)
