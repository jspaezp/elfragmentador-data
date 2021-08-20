from pyteomics.mzml import PreIndexedMzML

raw_file = "20200708_EXPL5_nLC3_DS_SA_HeLa_ProAlanase_IS_02"
spec_index = 516

mzml = PreIndexedMzML(f"raw/{raw_file}.mzML")

# mzml['controllerType=0 controllerNumber=1 scan=11832']

"""
# print([x for x in mzml['controllerType=0 controllerNumber=1 scan=11832']])
{
    'index': 11831,
    'id': 'controllerType=0 controllerNumber=1 scan=11832',
    'defaultArrayLength': 551,
    'scanList': {
        'count': 1,
        'scan': [{
            'scanWindowList': {
                'count': 1,
                'scanWindow': [{
                    'scan window lower limit': 100.0,
                    'scan window upper limit': 2014.24755859375}]},
            'scan start time': 19.047732521248, ## x['scanList']['scan'][0]['scan start time']
            # ['scanList']['scan'][0]['filter string']
            'filter string': 'FTMS + p NSI d Full ms2 651.0649@hcd30.00 [100.0000-2014.2476]',
            'preset scan configuration': 2.0,
            'ion injection time': 54.000001400709,
            '[Thermo Trailer Extra]Monoisotopic M/Z:': 650.7300730828281}],
        'no combination': ''},
    'precursorList': {
        'count': 1, # TODO assert this
        'precursor': [{
            'spectrumRef': 'controllerType=0 controllerNumber=1 scan=11820',
            'isolationWindow': {
                'isolation window target m/z': 651.064880371094,
                'isolation window lower offset': 0.649999976158,
                'isolation window upper offset': 0.649999976158,
                'ms level': 1},
            'selectedIonList': {
                'count': 1,
                'selectedIon': [{
                    'selected ion m/z': 650.730073082828,
                    'charge state': 3.0,
                    'peak intensity': 3156310.75}]},
	    'activation': { # ['precursorList']['precursor'][0]['actvation']
	        'beam-type collision-induced dissociation': '',
	        'collision energy': 30.0}
        }]
    },
    'MSn spectrum': '',
    'ms level': 2,
    'positive scan': '',
    'centroid spectrum': '',
    'base peak m/z': 129.1024756,
    'base peak intensity': 375724.96875,
    'total ion current': 7282736.0,
    'lowest observed m/z': 101.071174621582,
    'highest observed m/z': 1926.151000976563,
    'count': 2
}
"""

print(
    mzml["controllerType=0 controllerNumber=1 scan=11832"]["precursorList"][
        "precursor"
    ][0]["activation"]
)
"""
{'beam-type collision-induced dissociation': '', 'collision energy': 30.0}
"""

print(
    mzml["controllerType=0 controllerNumber=1 scan=11832"]["precursorList"][
        "precursor"
    ][0]["activation"]["collision energy"]
)
"""
30.0
"""

print(
    mzml["controllerType=0 controllerNumber=1 scan=11832"]["scanList"]["scan"][0][
        "filter string"
    ]
)
"""
FTMS + p NSI d Full ms2 651.0649@hcd30.00 [100.0000-2014.2476]
"""

print(
    mzml["controllerType=0 controllerNumber=1 scan=11832"]["scanList"]["scan"][0][
        "scan start time"
    ]
)
"""
19.047732521248
"""

from tqdm.auto import tqdm
import shlex
import time


def add_ce_info(mzml_directory, in_sptxt, out_sptxt):
    id_template = "controllerType=0 controllerNumber=1 scan={scan_number}"
    add_template = "RetentionTime={rt} CollisionEnergy={ce}"
    mzml_dict = {}
    with open(in_sptxt, "r") as infile:
        with open(out_sptxt, "w") as outfile:
            for i, l in enumerate(tqdm(infile)):
                l = l.strip()
                if l.startswith("Comment"):
                    l_list = l.strip().split(" ")
                    """
                    try:
                        l_list = shlex.split(l.strip())
                    except ValueError:
                        print(l)
                        print((
                            "Line did not have a closing quote"
                            " will fall back to string splitting"
                        ))
                        # import pdb ; pdb.set_trace()
                    """

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
                    l = " ".join(extra_left + [addition] + extra_right)

                outfile.write(l + "\n")


add_ce_info("raw/", "./mokapot_spectrast/MannLabPhospho.mokapot.psms.sptxt", "foo.txt")
