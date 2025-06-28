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
   echo "${file_000##*/}"

   test_rm_ith "filter.var"
   awk -F, 'BEGIN{OFS="\t"}{if(NR>1 && NF == 20){gsub(/"/,""); print $(NF-1), $NF}else if(NR>1 && NF == 19){print $NF, "NA"}}' ${file_000} \
     | sort \
     | uniq > filter.var

    test_rm_ith "${species}_all_centromeric_arrays.tsv"
    while read -r filt; do
      filt_1=$(echo ${filt} | awk '{print $1}')
      filt_2=$(echo ${filt} | awk '{print $2}')
      if ! grep -P -q "NA" filter.var; then
        awk -F, -v filt_1=${filt_1} -v filt_2=${filt_2} 'BEGIN{OFS="\t"}NR>1 && $(NF-1) == filt_1 && $NF == filt_2{gsub(/"/,""); print $1, $2, $3, $8, $4, filt_1, filt_2}' ${file_000} \
          | awk 'BEGIN{OFS="\t"}{if($2 < 0){$2 = 0; print}else{print}}' \
          | awk 'BEGIN{OFS="\t"}{$2=sprintf("%.0f",$2); $3=sprintf("%.0f",$3); print}' \
          | sort -k1,1 -k2,2n >> ${species}_all_centromeric_arrays.tsv
      else
        awk -F, -v filt_1=${filt_1} -v filt_2=${filt_2} 'BEGIN{OFS="\t"}NR>1 && $(NF) == filt_1{gsub(/"/,""); print $1, $2, $3, $8, $4, filt_1, filt_2}' ${file_000} \
          | awk 'BEGIN{OFS="\t"}{if($2 < 0){$2 = 0; print}else{print}}' \
          | awk 'BEGIN{OFS="\t"}{$2=sprintf("%.0f",$2); $3=sprintf("%.0f",$3); print}' \
          | sort -k1,1 -k2,2n >> ${species}_all_centromeric_arrays.tsv
      fi
    done < filter.var

  done < <(ls -htl *_centromeric_arrays.csv | awk '{print $NF}')
