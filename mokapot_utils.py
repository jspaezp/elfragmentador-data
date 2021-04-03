
import pandas as pd

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
    df = pd.read_table(str(input))
    df['Peptide'] = [x.replace("[", "[+") for x in df['Peptide']]
    df['Charge'] = ["_".join(line.split("_")[-2]) for line in df['SpecId']]
    df['Peptide'] = [f"{x}/{y}" for x, y in zip(df['Peptide'], df['Charge'])]
    df['SpecId'] = ["_".join(line.split("_")[:-3]) for line in df['SpecId']]
    df['mokapot PEP'] = 1-df['mokapot PEP']

    col_neworder = [list(df)[x] for x in index_order]
    df = df[col_neworder]
    df.to_csv(str(output), sep='\t', index=False)
    print(df)
