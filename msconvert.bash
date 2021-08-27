#!/usr/bin/zsh

set -x
set -e

real_path=$(realpath $1)
real_parent=$(dirname $real_path)
mzml_path="$(echo $real_path | sed -e 's/.raw/.mzML/g')"

if [[ ! -z "$(command -v docker)" ]] ; then
    echo "Docker found, running native docker"
    docker run \
        -it --rm -e WINEDEBUG=-all \
        -v $PWD/:/data chambm/pwiz-skyline-i-agree-to-the-vendor-licenses wine msconvert.exe \
        --zlib \
         -o /data/ --verbose \
        --filter "peakPicking true 1-" \
        --filter "activation HCD" \
        --filter "analyzer FT" /data/$1
else
    echo "Docker not found, attempting singfularity run"

    # Needs running mkdir /scratch/brown/jpaezpae/pwiz_sandbox/pwiz_sandbox/mywineprefix
    singularity exec \
      -B ${real_parent}/:/data \
      -B `mktemp -d /dev/shm/wineXXX`:/mywineprefix \
      -w /scratch/brown/jpaezpae/pwiz_sandbox/pwiz_sandbox \
      mywine msconvert \
         --zlib \
         --filter "peakPicking true 1-" \
         --filter "activation HCD" \
         --filter "analyzer FT" \
         -o /data/ --verbose \
         /data/$(basename $real_path)


    # singularity exec --writable \
    #     -S /mywineprefix/ ${CLUSTER_SCRATCH}/pwiz_sandbox/pwiz_sandbox \
    #     mywine msconvert \
    #     --zlib \
    #     --filter "peakPicking true 1-" \
    #     --filter "activation HCD" \
    #     --filter "analyzer FT" $1
fi


ls -lcth ${real_path}
# cat ${file_base} > ./raw/${file_base}
# rm -rf ${file_base}
ls -lcth ${mzml_path}
# touch "$(echo $1 | sed -e 's/.raw/.mzML/g')"
