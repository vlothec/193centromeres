#!/bin/bash

set -e
set -u
set -o pipefail

bin_path=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

## load lib
source ${bin_path}/wfun_lib.sh

for species in drHedHeli1 daBalNigr1 daSheArve1 daMisOron1 daLinVulg1; do

  # soloLTRseeker with TEsorter classification
  # loop for running it sequentially
  # change it to run them in parallel
  nohup /path/to/soloLTRseeker -t ${species}.1.fa.mod.EDTA.intact.gff3 ${species}.1.fa.gz &

done
