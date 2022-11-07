
from snakemake.utils import min_version
min_version("6.0")
from contextlib import contextmanager

from loguru import logger as lg_logger

import tempfile

import numpy as np
import pandas as pd
from pathlib import Path

@contextmanager
def temporary_links(files, target_dir):
    target_dir = Path(target_dir)
    linked_files = []
    for mzml in files:
        mzml = Path(mzml).resolve(strict=True)
        # The actual path for the link
        link = target_dir / Path(mzml).name
        link.symlink_to(mzml)
        lg_logger.debug(f"Created link {link} to {mzml}")
        linked_files.append(link)
    yield linked_files
    for link in linked_files:
        if link.is_symlink():
            link.unlink()
            lg_logger.debug(f"Removed link {link}")
        else:
            lg_logger.error(f"Link {link} is not a link anymore")
            raise RuntimeError(f"Link {link} is not a link anymore")

SPEC_ID_REGEX = re.compile(r"^(.*)_(\d+)_(\d+)_(\d+)$")

def _parse_line(line):
    outs = SPEC_ID_REGEX.match(line.SpecId).groups()
    file_name, spec_number, charge, _ = outs
    file_name = Path(file_name).stem + ".mzML"

    first_period = 0
    last_period = 0
    for i, x in enumerate(line.Peptide):
        if x == ".":
            if first_period == 0:
                first_period = i
            else:
                last_period = i

    sequence = line.Peptide[first_period + 1 : last_period]

    line_out = {
        "file": file_name,
        "scan": spec_number,
        "charge": charge,
        "sequence": sequence,
        "score-type": "PERCOLATOR QVALUE",
        "score": line["mokapot q-value"],
    }
    return line_out


def convert_to_ssl(input_file, output_file):
    lg_logger.info(f"Reading File {input_file}")
    df = pd.read_csv(
        input_file,
        sep="\t",
        dtype={
            "SpecId": str,
            "Label": bool,
            "ScanNr": np.int32,
            "ExpMass": np.float16,
            "CalcMass": np.float16,
            "Peptide": str,
            "mokapot score": np.float16,
            "mokapot q-value": np.float32,
            "mokapot PEP": np.float32,
            "Proteins": str,
        },
    )

    lg_logger.info("Processing File")
    out_df = pd.DataFrame([_parse_line(x) for _, x in df.iterrows()])
    lg_logger.info(f"Writting output: {output_file}")
    out_df.to_csv(output_file, sep="\t", index=False, header=True)
    lg_logger.info("Done")


def split_mokapot_psms(psms_file, spec_metadata_list, out_dir):
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
    metadata_df = pd.read_csv(spec_metadata_list, sep="\t", dtype={"ScanNr": np.int32})
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


# Note that this is a checkpoint, not a rule
checkpoint split_mokapot_psms:
    input:
        psms = "mokapot.psms.txt",
        spec_metadata=["spec_metadata.csv","spec_metadata2.csv"],
    output:
        split_files_dir=directory("split_mokapot/{experiment}"),
    run:
        shell(f"mkdir -p {output.split_files_dir}")
        print(input.spec_metadata)
        split_mokapot_psms(
            input.psms, input.spec_metadata, output.split_files_dir
        )

rule bibliospec:
    """Runs bibliospec to build the library
    Currently not implemented... because I cannot coimplile it in mac :(
    """
    input:
        psms = "results/{experiment}/mokapot/mokapot.psms.txt",
        mzML = lambda x: exp_files_map[x.experiment]['mzML']
    output:
        ssl_file = "results/{experiment}/bibliospec/{experiment}.ssl",
        library_name = "results/{experiment}/bibliospec/{experiment}.blib",
    run:
        lg_logger.info("Converting psms to ssl")
        convert_to_ssl(input.psms, output.ssl_file)
        lg_logger.info("Done Converting psms to ssl")

        shell_cmd = [
            "LC_ALL='C' BlibBuild "
            "-C 2G", # minimum size to start caching
            "-c 0.99",
            "-m 500M", # sqlite cache size
            "-A", # warns ambiguous spectra
            str(output.ssl_file),
            str(output.library_name),
        ]
        shell_cmd = " ".join(shell_cmd)
        lg_logger.info("Running: {shell_cmd}")
        # This creates soft links in the so bibliospec can find the raw spectra

        out_parent_dir = Path(str(output.ssl_file)).parent
        with temporary_links(input, str(out_parent_dir)) as linked_files:
            lg_logger.info("Running bibliospec: %s", shell_cmd)
            shell(shell_cmd)

rule filter_bibliospec:
    input:
        library = "results/{experiment}/bibliospec/{experiment}.blib",
    output:
        filtered_library = "results/{experiment}/bibliospec/{experiment}.filtered.blib",
    run:
        """
        BlibFilter [options] <redundant-library> <filtered-library> 
        -m [ --memory-cache ] <size> SQLite memory cache size in Megs. Default 250M.
        -n [ --min-peaks ] <num> Only include spectra with at least this many peaks. Default 20.
        -s [ --min-score ] <score> Best spectrum must have at least this average score to be included. Default 0.
        -p [ --parameter-file ] <file> File containing search parameters. Command line values override file values.
        -v [ --verbosity ] <level> Control the level of output to stderr. (silent, error, status, warn, debug, detail, all) Default status.
        """

        shell(
            "LC_ALL='C' BlibFilter "
            "--memory-cache 500M "
            "--min-peaks 20 "
            "--min-score 0 " # Filtering should be done creating the first lib
            "-v status "
            "{input.library}"
            "{output.filtered_library} "
        )

