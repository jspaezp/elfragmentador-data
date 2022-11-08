import pandas as pd
from ms2ml import Peptide
from ms2ml.data.parsing.bibliospec import _decompress_peaks
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

import re
# wtf ... ok ... the .ssl should move all terminal mods to the
# first aminoacid ...
#        'N[+245.2]M[+229.2]KPTGTDPRILSIAAEVA',
#         n[229.1629]MK[229.1629]PTGTDPRILSIAAEVA
#         n[229.1629]M[15.9949]K[229.1629]PTGTDPRILSIAAEV

#       Cool, found a bug in bibliospec
#       'N[+458.3]KAAAAAAALQA'
#        n[229.1629]K[229.1629]AAAAAAALQA

special_replacements = {
    "N[+458.3258]": "[+229.1629]-N[+229.1629]",
}

unique_mods = np.unique(pd.read_parquet('bibliospec_tables/AspN/mods.parquet')['mass'])
mod_alias_dict = {f"[+{x:.01f}]": f"[+{x:.04f}]" for x in unique_mods}

spec_meta_df = pd.read_parquet("bibliospec_tables/Kmod_Dimethyl/spec_meta.parquet")[["RawFile", "CollisionEnergy", "ScanNr"]]
spec_sourcefiles = pd.read_parquet("bibliospec_tables/Kmod_Dimethyl/spec_sourcefiles.parquet")
spec_data_df = pd.merge(spec_meta_df, spec_sourcefiles)

spec_df = pd.read_parquet("bibliospec_tables/Kmod_Dimethyl/spec.parquet")
spec_df["SpecIDinFile"] = spec_df["SpecIDinFile"].astype(int)

df = pd.merge(spec_data_df, spec_df, left_on=["id", "ScanNr"], right_on=["fileID", "SpecIDinFile"])
for i in df.groupby(["peptideModSeq", "CollisionEnergy"]):
    mod_seq = i.peptideModSeq
    for k, v in mod_alias_dict.items():
        mod_seq = mod_seq.replace(k, c)
    break

[_decompress_peaks(x,y,z) for x, y, z in tqdm(zip(df["peakMZ"], df["peakIntensity"], df["numPeaks"]))]
