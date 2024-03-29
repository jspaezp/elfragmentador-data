import subprocess
from spec_metadata import get_spec_metadata


rule download_file:
    output:
        "raw/{sample}.raw",
    run:
        server = samp_to_ftp[wildcards.sample]
        rawsample = sample_to_raw[wildcards.sample]
        shell("mkdir -p raw")
        shell(
            (
                f"wget '{server}" + "/{rawsample}.raw' -O "
                "'{output}' -nc "
                "--timeout=15 "
                "--limit-rate=50m "
                "--wait=5 "
            )
        )


rule convert_file:
    input:
        ancient("raw/{sample}.raw"),
    output:
        "raw/{sample}.mzML",
    benchmark:
        "benchmarks/{sample}.conversion.benchmark.txt"
    run:
        # For some reaso, the dockerized version fails when running it directly
        # in this script, so you have to hack it this way ...
        cmd = ["zsh", "msconvert.bash", str(input).replace(" ", "\ ")]
        print(cmd)
        subprocess.run(cmd, check=True, capture_output=True)


rule mzml_scan_metadata:
    input:
        "raw/{sample}.mzML",
    output:
        "raw_scan_metadata/{sample}.csv",
    run:
        shell("mkdir -p raw_scan_metadata")
        df = get_spec_metadata(str(input))
        df.to_csv(str(output), index=False)
