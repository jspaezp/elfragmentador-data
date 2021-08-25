
from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np 
import spectrum_utils.plot as sup
import spectrum_utils.spectrum as sus
from pyteomics import mzml

import elfragmentador as ef
import elfragmentador.model as efm

from elfragmentador.spectra import Spectrum

model = efm.PepTransformerModel.load_from_checkpoint(ef.DEFAULT_CHECKPOINT)

def read_spectrum_number(scan_id, file, sequence):
    reader = mzml.PreIndexedMzML(file)
    mzml_index = f'controllerType=0 controllerNumber=1 scan={scan_id}'
    spectrum = reader.get_by_id(mzml_index)
    charge = int(spectrum["precursorList"]["precursor"][0]["selectedIonList"]["selectedIon"][0]['charge state'])
    nce = float(spectrum["precursorList"]["precursor"][0]["activation"]["collision energy"])
    precursor_mz = spectrum["precursorList"]["precursor"][0]["selectedIonList"]["selectedIon"][0]['selected ion m/z']

    spec = Spectrum(
        sequence=sequence,
        charge=charge,
        nce=nce,
        parent_mz=precursor_mz,
        mzs = spectrum['m/z array'],
        intensities = spectrum['intensity array'])
        
    return spec


def plot_scan_mirror(scan_id, peptide, file, scaling="root", nce_offset = 0, prefix: str=""):
    spectrum_top = read_spectrum_number(scan_id=scan_id, sequence=peptide, file=file)

    # note that I am adding 5 to the nce because the instrument is a QE and empirically
    # there is an offser in 5 with respect to the fusion lumos NCE.
    prediction = model.predict_from_seq(peptide, spectrum_top.charge, nce=spectrum_top.nce + nce_offset, as_spectrum=True)

    fig, ax = plt.subplots(figsize=(7, 5), constrained_layout=True)
    spectrum_top.plot(
        mirror=prediction,
        ax=ax)

    plt.title(f"{peptide}/{spectrum_top.charge}+, ScanID={scan_id}, mz={spectrum_top.parent_mz:.2f}, scaling={scaling}")

    save_name = f"{prefix}_{peptide}_{scan_id}.png"
    print(save_name)
    plt.savefig(save_name)
    plt.close()


def clean_file_name(fname):
    fname = fname[:-fname[::-1].index("_")-1]
    fname = fname[:-fname[::-1].index("_")-1]
    fname = fname[:-fname[::-1].index("_")-1]
    return fname


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--csv")
    parser.add_argument("--search_path")
    parser.add_argument("--prefix")

    args = parser.parse_args()
    df = pd.read_csv(args.csv)

    for scan_id, peptide, spec_id in zip(df["ScanNr"], df["Peptide"], df["SpecId"]):
        file = clean_file_name(spec_id)
        file = Path(args.search_path) / (file + ".mzML")
        assert file.is_file(), file

        peptide = peptide[peptide.index(".")+1:]
        peptide = peptide[:peptide.index(".")]

        plot_scan_mirror(scan_id, peptide, str(file), None, prefix=args.prefix)


