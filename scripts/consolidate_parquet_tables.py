from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import uniplot
from loguru import logger as lg_logger
from ms2ml import AnnotatedPeptideSpectrum, Config, Peptide, Spectrum
from ms2ml.data.parsing.bibliospec import _decompress_peaks
from ms2ml.landmarks import IRT_PEPTIDES
from sklearn.linear_model import HuberRegressor
from tqdm.auto import tqdm

"""
bibliospec_tables/AspN/mods.parquet
[15.9949]
bibliospec_tables/FirstPoolSemi/mods.parquet
[15.9949]
bibliospec_tables/GGMannLab/mods.parquet
[ 15.9949 114.0429]
bibliospec_tables/HelaProalanase/mods.parquet
[15.9949]
bibliospec_tables/HLA_1/mods.parquet
[15.9949]
bibliospec_tables/HLA2_1/mods.parquet
[15.9949]
bibliospec_tables/HLA2_2/mods.parquet
[15.9949]
bibliospec_tables/HLA_2/mods.parquet
[15.9949]
bibliospec_tables/HLA_3/mods.parquet
[15.9949]
bibliospec_tables/Human_Hela/mods.parquet
[15.9949]
bibliospec_tables/HumanTestisPhospho/mods.parquet
[15.9949 79.9663]
bibliospec_tables/Kmod_Acetyl/mods.parquet
[15.9949 42.0106]
bibliospec_tables/Kmod_Biotinyl/mods.parquet
[ 15.9949 226.0776]
bibliospec_tables/Kmod_Dimethyl/mods.parquet
[15.9949 28.0313]
bibliospec_tables/Kmod_Formyl/mods.parquet
[15.9949 27.9949]
bibliospec_tables/Kmod_Methyl/mods.parquet
[14.0157 15.9949]
bibliospec_tables/Kmod_Trimethyl/mods.parquet
[15.9949 42.047 ]
bibliospec_tables/Kmod_Ubiquitinyl/mods.parquet
[ 15.9949 114.0429]
bibliospec_tables/Kmod_Unmod/mods.parquet
[15.9949]
bibliospec_tables/LysN/mods.parquet
[15.9949]
bibliospec_tables/MannLabPhospho2/mods.parquet
[15.9949 79.9663]
bibliospec_tables/MannLabPhospho/mods.parquet
[15.9949 79.9663]
bibliospec_tables/Pmod_Hydroxyproline/mods.parquet
[15.9949]
bibliospec_tables/ProteomeToolsFP2/mods.parquet
[15.9949]
bibliospec_tables/ProteomeToolsFP/mods.parquet
[15.9949]
bibliospec_tables/ProteomeToolsIF2/mods.parquet
[15.9949]
bibliospec_tables/ProteomeToolsIF/mods.parquet
[15.9949]
bibliospec_tables/ProteomeToolsMF/mods.parquet
[15.9949]
bibliospec_tables/ProteomeToolsPTMT/mods.parquet
[15.9949]
bibliospec_tables/ProteomeToolsSA/mods.parquet
[15.9949]
bibliospec_tables/ProteomeToolsSP2/mods.parquet
[15.9949]
bibliospec_tables/ProteomeToolsSP/mods.parquet
[15.9949]
bibliospec_tables/ProteomeToolsSRMP/mods.parquet
[15.9949]
bibliospec_tables/ProteomeToolsTP/mods.parquet
[15.9949]
bibliospec_tables/Rmod_Citrullin/mods.parquet
[ 0.984  15.9949]
bibliospec_tables/Rmod_Dimethyl_asymm/mods.parquet
[15.9949 28.0313]
bibliospec_tables/Rmod_Dimethyl_symm/mods.parquet
[15.9949 28.0313]
bibliospec_tables/Rmod_Methyl/mods.parquet
[14.0157 15.9949]
bibliospec_tables/Rmod_Unmod/mods.parquet
[15.9949]
bibliospec_tables/SecondPoolSemi/mods.parquet
[15.9949]
bibliospec_tables/TMT_ASPN/mods.parquet
[ 15.9949 229.1629 245.1578 458.3258]
bibliospec_tables/TMT_HLA_1/mods.parquet
[ 15.9949 229.1629 245.1578 458.3258]
bibliospec_tables/TMT_HLA2_1/mods.parquet
[ 15.9949 229.1629 245.1578 458.3258]
bibliospec_tables/TMT_HLA2_2/mods.parquet
[ 15.9949 229.1629 245.1578 458.3258]
bibliospec_tables/TMT_HLA_2/mods.parquet
[ 15.9949 229.1629 245.1578 458.3258]
bibliospec_tables/TMT_LYSN/mods.parquet
[ 15.9949 229.1629 245.1578 458.3258]
bibliospec_tables/Ymod_Nitrotyr/mods.parquet
[15.9949 44.9851]
bibliospec_tables/Ymod_Phospho/mods.parquet
[15.9949 79.9663]
bibliospec_tables/Ymod_Unmod/mods.parquet
[15.9949]
"""

# wtf ... ok ... the .ssl should move all terminal mods to the
# first aminoacid ...
#        'N[+245.2]M[+229.2]KPTGTDPRILSIAAEVA',
#         n[229.1629]MK[229.1629]PTGTDPRILSIAAEVA
#         n[229.1629]M[15.9949]K[229.1629]PTGTDPRILSIAAEV

#       Cool, found a bug in bibliospec, reported ...
#       'N[+458.3]KAAAAAAALQA'
#        n[229.1629]K[229.1629]AAAAAAALQA


# But translating to unimod is helpful ...
mod_alias_dict = {
    "[+0.984]": "[U:7]",
    "[+15.9949]": "[U:35]",
    "[+79.9663]": "[U:21]",
    "[+28.0313]": "[U:36]",  # dimethyl 36
    "[+44.9851]": "[U:354]",  # nitro 354
    "[+229.1629]": "[U:737]",  # tmt6 737
    "[+14.0157]": "[U:34]",  # methyl 34
    "[+114.0429]": "[U:121]",  # GG 121
    "[+42.0106]": "[U:1]",  # Acetyl
    "[+42.047]": "[U:37]",  # Acetyl
}

import re

TERM_MOD = re.compile("^(.)(\[\+[0-9.]+\])(.*)$")
TEMPLATES = {
    "[+229.1629]": "[+229.1629]-{aa}",
    "[+458.3258]": "[+229.1629]-{aa}[+229.1629]",
    "[+245.1578]": "[+229.1629]-{aa}[+15.9949]",
}
EXCLUSIONS = {
    ("K", "[+229.1629]"),
}


def move_terminal_mods(seq):
    """Move terminal mods to first aminoacid"""
    match = TERM_MOD.match(seq)
    if not match:
        return seq

    aa, mod, rest = match.groups()
    if mod in TEMPLATES:
        if (aa, mod) in EXCLUSIONS:
            return seq
        mod = TEMPLATES[mod].format(aa=aa)
        seq = mod + rest
    return seq


def test_move_terminal_mods():
    assert move_terminal_mods("K[+229.1629]AA") == "K[+229.1629]AA"
    assert move_terminal_mods("K[+458.3258]AA") == "[+229.1629]-K[+229.1629]AA"
    assert move_terminal_mods("K[+229.1629]AA") == "K[+229.1629]AA"
    assert move_terminal_mods("M[+245.1578]AA") == "[+229.1629]-M[+15.9949]AA"


def main(base_path):
    lg_logger.info("Starting the processing of spectra ... ")
    lg_logger.info("Using the following parameters: ")
    lg_logger.info("base_path: %s" % base_path)

    spec_meta_df = pd.read_parquet(base_path / "spec_meta.parquet")[
        ["RawFile", "CollisionEnergy", "ScanNr"]
    ]
    spec_sourcefiles = pd.read_parquet(base_path / "spec_sourcefiles.parquet")
    spec_data_df = pd.merge(spec_meta_df, spec_sourcefiles)

    spec_df = pd.read_parquet(base_path / "spec.parquet")
    spec_df["SpecIDinFile"] = spec_df["SpecIDinFile"].astype(int)

    df = pd.merge(
        spec_data_df,
        spec_df,
        left_on=["id", "ScanNr"],
        right_on=["fileID", "SpecIDinFile"],
    )
    # Set this for interactive use
    pd.set_option("display.max_columns", None)
    lg_logger.info(df)

    split_lgl = [x in IRT_PEPTIDES for x in df.peptideModSeq]
    irt_df = df[split_lgl].copy()
    df = df[[not x for x in split_lgl]].copy()

    irt_df["irt"] = [IRT_PEPTIDES[x]["irt"] for x in irt_df.peptideModSeq]
    if (found_peps := len(np.unique(irt_df.peptideModSeq))) > 3:
        lg_logger.info(
            f"Fitting huber regressor to irt peptides, found {found_peps} unique"
            " peptides"
        )
        my_svr = HuberRegressor()
        my_svr = my_svr.fit(X=irt_df.retentionTime.values.reshape(-1, 1), y=irt_df.irt)

        df["pred_irt"] = my_svr.predict(df.retentionTime.values.reshape(-1, 1))

        uniplot.plot(df["pred_irt"], df["retentionTime"])
        uniplot.plot(irt_df["irt"], irt_df["retentionTime"], color="red")
    else:
        lg_logger.info(
            f"Skipping irt calibration, found onlu {found_peps} uniqueirt peptides"
        )
        df["pred_irt"] = float("nan")

    config = Config.from_toml("ms2ml_config.toml")

    """
    I need these cols:
        Parameters:
            seq (Tensor): Encoded peptide sequences
            mods (Tensor): Modification encodings for the sequence
            charge (Tensor): Long tensor with the charges
            nce (Tensor): Normalized collision energy
            irt (Tensor): Tensor containing normalized irt predictions
            spectra (Tensor): Tensor containing encoded predicted spectra
            weight (Union[Tensor, None]): Weight of the element for the loss
    """

    outs = []
    skipped_sequnces = []

    MOD_REGEX = re.compile(r"(.?\[[A-Z]?[+0-9:]+\]-?)")
    mod_count = defaultdict(int)

    for (modseq, nce, charge, ims), x_df in tqdm(
        df.groupby(
            ["peptideModSeq", "CollisionEnergy", "precursorCharge", "ionMobility"]
        )
    ):
        # If needed the peptide here van be recycled to speedup the computations
        # basicaly doing a nested for loop where the peptide is shared between the
        # different collision energies
        out_dict = {}
        modseq = move_terminal_mods(modseq)
        for k, v in mod_alias_dict.items():
            modseq = modseq.replace(k, v)

        for x in MOD_REGEX.findall(modseq):
            mod_count[x] += 1

        try:
            if len(x_df) < 3:
                continue
            pep = Peptide.from_sequence(f"{modseq}/{charge}", config=config)
            out_dict["seq"] = pep.aa_to_vector()
            out_dict["mods"] = pep.mod_to_vector()
            out_dict["charge"] = np.array([charge], dtype=np.int32)
            decompressed = [
                (*_decompress_peaks(x, y, z), pmz)
                for x, y, z, pmz in zip(
                    x_df["peakMZ"],
                    x_df["peakIntensity"],
                    x_df["numPeaks"],
                    x_df["precursorMZ"],
                )
            ]
            specs = [
                Spectrum(
                    mz=mz,
                    intensity=intensity,
                    ms_level=2,
                    precursor_mz=pmz,
                    config=config,
                )
                for mz, intensity, pmz in decompressed
            ]
            specs = [spec.annotate(pep) for spec in specs]
            specs = [spec.encode_fragments() for spec in specs]
            specs = [spec / (spec.max() + 1e-8) for spec in specs]
            out_dict["nce"] = np.array([nce], dtype=np.float32)
            out_dict["spectra"] = np.stack(specs).mean(axis=0)
            out_dict["weight"] = np.log1p(np.array([len(specs)]))
            out_dict["irt"] = np.array([x_df["pred_irt"].mean()], dtype=np.float32)
            out_dict["ims"] = np.array([ims], dtype=np.float32)
            outs.append(out_dict)
        except ValueError:
            skipped_sequnces.append(modseq)
            lg_logger.warning(f"Skipped {modseq}")

    lg_logger.info(f"Found {len(outs)} spectra")
    lg_logger.info(f"Skipped {len(skipped_sequnces)} sequences")
    lg_logger.info(f"{skipped_sequnces}")
    lg_logger.info(f"Found {len(mod_count)} unique modifications")
    lg_logger.info(f"{mod_count}")

    out_df = pd.DataFrame(outs)
    lg_logger.info(f"{out_df.head()}")
    out_df.to_parquet(base_path / "processed.parquet")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("base_path", type=Path)
    args = parser.parse_args()
    main(args.base_path)
