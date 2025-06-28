#!/bin/bash

set -e
set -u
set -o pipefail

  bin_path=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

  ## load lib
  source ${bin_path}/wfun_lib.sh

  while read -r file_001; do
    Rscript gap_infer_prop_all_gaps.R ${file_001}
  done < <(ls -htl *.TEanno.split.all.gaps | tr -s ' ' | awk '$5 > 0{print $NF}')
