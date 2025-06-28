#!/bin/bash

set -e
set -u
set -o pipefail

bin_path=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

## load lib
source ${bin_path}/wfun_lib.sh

# get monocot and dicot species list
ls d*_mod.EDTA.intact.gff3.array_boundaries | awk -F/ '{print $NF}' | awk -F "." '{print $1}' > species.list
ls l*_mod.EDTA.intact.gff3.array_boundaries | awk -F/ '{print $NF}' | awk -F "." '{print $1}' >> species.list

  while read -r file_001; do
    species=$(echo ${file_001} | awk -F/ '{print $NF}' | awk -F "." '{print $1}')
    file_002=$(ls -htl ${species}*_mod.EDTA.intact.gff3.array_boundaries | awk '{print $NF}')
    Rscript sat_invasion_tally.R ${file_002} ${file_001}
  done < <(ls -htl *.mod.EDTA.intact.fa.rexdb-plant.cls.tsv | awk '{print $NF}' | grep -f species.list)
