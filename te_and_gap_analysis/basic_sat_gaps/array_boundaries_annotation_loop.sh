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

   # comment out as needed
   # reassigned
   # file_001="${species}_edta_filtered.csv.reassigned"
   # intact
   file_001="${species}.mod.EDTA.intact.gff3"

   Rscript array_boundaries_annotation.R "${file_000}" "${file_001}" > ${file_001##*\/}_array_boundaries_annotation.log

  done < <(ls -htl *fa_all_centromeric_arrays.tsv | awk '{print $NF}')
