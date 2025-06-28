#!/bin/bash

set -e
set -u
set -o pipefail

bin_path=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

## load lib
source ${bin_path}/wfun_lib.sh


while read -r file_000; do
  # cent array file
  species=$(echo ${file_000##*/} | awk -F "_" '{print $1}')
  echo "${species}"

  Rscript structural_homology_ratio.R "${file_000}" > ${species}_structural_homology_ratio.log

done < <(ls -htl *.fa_edta_filtered.csv.reassigned.array_boundaries | awk '{print $NF}')
