
from pyteomics.mzml import PreIndexedMzML
import pandas as pd
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


def add_spectrast_ce_info(mzml_directory, in_sptxt, out_sptxt):
    """
    add_ce_info("raw/", "./mokapot_spectrast/MannLabPhospho.mokapot.psms.sptxt", "foo.txt")
    """
    id_template = 'controllerType=0 controllerNumber=1 scan={scan_number}'
    add_template = 'RetentionTime={rt} CollisionEnergy={ce}'
    mzml_dict = {}
    with open(in_sptxt, 'r') as infile:
        with open(out_sptxt, 'w') as outfile:
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

