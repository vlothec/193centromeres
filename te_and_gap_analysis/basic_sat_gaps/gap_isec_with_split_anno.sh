#!/bin/bash

set -e
set -u
set -o pipefail

  bin_path=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

  ## load lib
  source ${bin_path}/wfun_lib.sh

  # intersect gaps
  while read -r file_001; do
    # set var
    species=$(echo ${file_001/\.fa*/})
    file_002=$(ls -htl ${species}_edta_filtered.csv.reassigned | awk '{print $NF}')

    # sort input file 1
    awk 'NR>1' ${file_001} \
       | awk 'BEGIN{OFS="\t"}{if($NF !~ /_/){$NF=NR"_"$NF; print}else{print}}' \
       | sort -k1,1 -k2,2n > temp_file_001
    mv temp_file_001 ${file_001}
    # sort input file 2

    awk -F, '$10 !~ /repeat_region_/' ${file_002} \
      | awk -F, 'BEGIN{OFS="\t"}{print $2, $5, $6, $4, $7, $8}' \
      | awk 'BEGIN{OFS="\t"}{if($2 < 0){$2 = 0; print}else{print}}' \
      | awk 'BEGIN{OFS="\t"}{$2=sprintf("%.0f",$2); $3=sprintf("%.0f",$3); print}' \
      | sort -k1,1 -k2,2n > temp_file_001

    echo ${species}

    # intersect
    bedtools intersect -a ${file_001} -b temp_file_001 -wao > ${file_001}.TEanno.split.all.gaps

    # sum
    # wc -l ${file_001}.TEanno.split >> TEanno.split.within.gaps.all
  done < <(ls -htl *.bed.gaps | awk '{print $NF}')

