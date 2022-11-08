
from snakemake.utils import min_version
min_version("6.0")
from contextlib import contextmanager

from loguru import logger as lg_logger

import tempfile

import numpy as np
import pandas as pd
from pathlib import Path
from tqdm.auto import tqdm
tqdm.pandas()

@contextmanager
def temporary_links(files, target_dir):
    target_dir = Path(target_dir)
    linked_files = []
    try:
      for mzml in files:
          mzml = Path(mzml).resolve(strict=True)
          # The actual path for the link
          link = target_dir / Path(mzml).name
          link.symlink_to(mzml)
          lg_logger.debug(f"Created link {link} to {mzml}")
          linked_files.append(link)
      yield linked_files
    finally:
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



rule bibliospec:
    """Runs bibliospec to build the library
    Currently not implemented... because I cannot coimplile it in mac :(
    """
    input:
        psms = "results/{experiment}/mokapot/mokapot.psms.txt",
        peptides = "results/{experiment}/mokapot/mokapot.peptides.txt",
        mzML = lambda x: exp_files_map[x.experiment]['mzML']
    output:
        ssl_file = "results/{experiment}/bibliospec/{experiment}.ssl",
        library_name = "results/{experiment}/bibliospec/{experiment}.blib",
    run:
        lg_logger.info("Converting psms to ssl")
        convert_to_ssl(input.psms, input.peptides, output.ssl_file)
        lg_logger.info("Done Converting psms to ssl")

        shell_cmd = [
            "LC_ALL='C' BlibBuild ",
            "-H ", # More than one decimal to describe mod
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

def convert_to_ssl(input_file, input_peptides, output_file):
    lg_logger.info(f"Reading File {input_peptides}")
    pep_df = pd.read_csv(input_peptides, sep="\t")
    pep_df = pep_df[pep_df["mokapot q-value"] < 0.01]
    peps = set(pep_df["Peptide"])

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

    lg_logger.info(f"Processing File {df.shape}")
    df = df[df["mokapot q-value"] < 0.01]
    lg_logger.info(f"Processing File {df.shape} filtered pval")
    df = df[[pep in peps for pep in tqdm(df["Peptide"])]]
    lg_logger.info(f"Processing File {df.shape} filtered peptides pval")
    df2 = []
    for i, x in tqdm(df.groupby('Peptide')):
      tmp = x.sort_values(by = 'mokapot score', ascending = False).head(20)
      df2.append(tmp)

    df = pd.concat(df2)

    lg_logger.info(f"Processing File {df.shape} filtered top")
    out_df = pd.DataFrame([_parse_line(x) for _, x in tqdm(df.iterrows())])

    lg_logger.info(f"Writting output: {output_file}")
    out_df.to_csv(output_file, sep="\t", index=False, header=True)

    lg_logger.info("Done")

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
