from pathlib import Path

import elfragmentador as ef
import elfragmentador.model as efm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import spectrum_utils.plot as sup
import spectrum_utils.spectrum as sus
from elfragmentador.annotate import get_theoretical_mass
from elfragmentador.spectra import Spectrum
from pyteomics import mzml


def read_spectrum_number(scan_id, file, sequence):
    reader = mzml.PreIndexedMzML(file)
    mzml_index = f"controllerType=0 controllerNumber=1 scan={scan_id}"
    spectrum = reader.get_by_id(mzml_index)

    nce = float(
        spectrum["precursorList"]["precursor"][0]["activation"]["collision energy"]
    )
    precursor_mz = spectrum["precursorList"]["precursor"][0]["selectedIonList"][
        "selectedIon"
    ][0]["selected ion m/z"]

    try:
        charge = int(
            spectrum["precursorList"]["precursor"][0]["selectedIonList"]["selectedIon"][
                0
            ]["charge state"]
        )
    except KeyError as e:
        charge = round(get_theoretical_mass(sequence) / precursor_mz)

    try:
        spec = Spectrum(
            sequence=sequence,
            charge=charge,
            nce=nce,
            parent_mz=precursor_mz,
            mzs=spectrum["m/z array"],
            intensities=spectrum["intensity array"],
        )
    except ValueError as e:
        print(sequence)
        raise ValueError(e)

    return spec


def plot_scan_mirror(
    scan_id, peptide, file, scaling="root", nce_offset=0, prefix: str = ""
):
    spectrum_top = read_spectrum_number(scan_id=scan_id, sequence=peptide, file=file)

    # note that I am adding 5 to the nce because the instrument is a QE and empirically
    # there is an offser in 5 with respect to the fusion lumos NCE.
    prediction = model.predict_from_seq(
        peptide,
        spectrum_top.charge,
        nce=spectrum_top.nce + nce_offset,
        as_spectrum=True,
        enforce_length=False,
    )

    fig, ax = plt.subplots(figsize=(7, 5), constrained_layout=True)
    spectrum_top.plot(mirror=prediction, ax=ax)

    plt.title(
        f"{peptide}/{spectrum_top.charge}+, ScanID={scan_id},"
        f" mz={spectrum_top.parent_mz:.2f}, scaling={scaling}"
    )

    save_name = f"{prefix}_{peptide}_{scan_id}.png"
    print(save_name)
    plt.savefig(save_name)
    plt.close()


def clean_file_name(fname):
    fname = fname[: -fname[::-1].index("_") - 1]
    fname = fname[: -fname[::-1].index("_") - 1]
    fname = fname[: -fname[::-1].index("_") - 1]
    return fname


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--csv")
    parser.add_argument("--search_path")
    parser.add_argument("--prefix")
    parser.add_argument("--checkpoint")

    args = parser.parse_args()

    model = efm.PepTransformerModel.load_from_checkpoint(args.checkpoint)
    df = pd.read_csv(args.csv)

    for i, (scan_id, peptide, spec_id) in enumerate(
        zip(df["ScanNr"], df["Peptide"], df["SpecId"])
    ):
        print(peptide)
        file = clean_file_name(spec_id)
        file = Path(args.search_path) / (file + ".mzML")
        assert file.is_file(), file

        peptide = peptide[peptide.index(".") + 1 :]
        peptide = peptide[: -peptide[::-1].index(".") - 1]

        plot_scan_mirror(
            scan_id, peptide, str(file), None, prefix=str(args.prefix) + f"_{i+1}"
        )
