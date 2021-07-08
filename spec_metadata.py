
from pathlib import Path
import pandas as pd
from pyteomics.mzml import PreIndexedMzML
from tqdm.auto import tqdm


def get_spec_metadata(mzml_path: str) -> pd.DataFrame:
    """
    raw_file = "20200708_EXPL5_nLC3_DS_SA_HeLa_ProAlanase_IS_02"
    mzml_path = f"raw/{raw_file}.mzML"
    get_spec_metadata(mzml_path)
    """
    filename = Path(mzml_path).parent / Path(mzml_path).stem
    mzml = PreIndexedMzML(mzml_path)
    df_dict = {
        "ScanNr": [],
        "CollisionEnergy": [],
        "RetentionTime": [],
        "SpecId": str(filename)
    }

    for spec in tqdm(mzml):
        if spec['ms level'] == 1:
            continue

        id = spec['id']
        controller, cont_number, scan = [x.split('=')[1] for x in id.split(' ')]
        ce = spec['precursorList']['precursor']
        ce = ce[0]['activation']['collision energy']
        rt = spec['scanList']['scan'][0]['scan start time']

        df_dict['RetentionTime'].append(rt)
        df_dict['CollisionEnergy'].append(ce)
        df_dict['ScanNr'].append(scan)

    return(pd.DataFrame(df_dict))

if __name__ == "__main__":
    raw_file = "20200708_EXPL5_nLC3_DS_SA_HeLa_ProAlanase_IS_02"
    mzml_path = f"raw/{raw_file}.mzML"
    print(get_spec_metadata(mzml_path))
