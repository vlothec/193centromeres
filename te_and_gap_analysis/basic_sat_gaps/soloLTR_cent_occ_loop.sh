
#!/bin/bash

set -e
set -u
set -o pipefail

bin_path=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

## load lib
source ${bin_path}/wfun_lib.sh

  while read -r file_001; do
    species=$(echo ${file_001} | awk -F/ '{print $NF}' | awk -F "." '{print $1}')
    file_002=$(ls -htl ${species}*.fa_all_centromeric_arrays.tsv | awk '{print $NF}')

    Rscript soloLTR_cent_occ.R ${file_001} ${file_002}

  done < <(ls -htl *fa_wga_soloLTR.gff3 | awk '{print $NF}')
