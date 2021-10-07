import shutil
import pathlib

curr_dir = str(pathlib.Path(".").absolute())

if shutil.which("docker"):
    TPP_DOCKER = f"docker run -v {curr_dir}/:/data/ spctools/tpp "
else:
    TPP_DOCKER = f"singularity exec /scratch/brown/jpaezpae/opt/tpp.img "
