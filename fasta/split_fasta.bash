#!bin/bash

set -e
set -x

inpath=${PWD}/$1
outpath=$(echo $1 | sed -e "s/\..*//g")

mkdir -p ${outpath}
cd ${outpath}
csplit -s -z ${inpath} '/>/' '{*}'
for i in xx* ; do \
  n=$(sed 's/>// ; s/ .*// ; 1q' "$i") ; \
  mv "$i" "$n.fasta" ; \
done
