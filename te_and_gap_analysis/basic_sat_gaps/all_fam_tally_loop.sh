#!/bin/bash

set -e
set -u
set -o pipefail

bin_path=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

## load lib
source ${bin_path}/wfun_lib.sh

 for species in drHedHeli1 daBalNigr1 daSheArve1 daMisOron1 daLinVulg1 ddEupPepu3; do
  Rscript all_fam_tally.R "${species}"
 done
