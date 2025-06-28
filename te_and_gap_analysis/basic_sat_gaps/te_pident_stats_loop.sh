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

   Rscript te_pident_stats_intact.R "${file_000}" > ${species}_te_pident_stats_intact.log

  done < <(ls -htl *fa_mod.EDTA.intact.gff3.array_boundaries | awk '{print $NF}')
