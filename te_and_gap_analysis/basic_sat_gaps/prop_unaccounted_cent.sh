#!/bin/bash

set -e
set -u
set -o pipefail

  bin_path=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

  ## load lib
  source ${bin_path}/wfun_lib.sh

  # prop of centromere unaccounted for

  while read -r file_000; do
   # cent array bed file
   species=$(echo ${file_000##*/} | awk -F "_" '{print $1}')
   awk -F, 'BEGIN{OFS="\t"}{if(NR>1 && NF == 20){gsub(/"/,""); print $(NF-1), $NF}else if(NR>1 && NF == 19){print $NF, "NA"}}' ${file_000} \
     | sort \
     | uniq > filter.var
   while read -r filt; do
     filt_1=$(echo ${filt} | awk '{print $1}')
     filt_2=$(echo ${filt} | awk '{print $2}')
     if ! grep -P -q "NA" filter.var; then
       awk -F, -v filt_1=${filt_1} -v filt_2=${filt_2} 'BEGIN{OFS="\t"}NR>1 && $(NF-1) == filt_1 && $NF == filt_2{gsub(/"/,""); print $1, $2-1, $3, $8, $4}' ${file_000} \
         | awk 'BEGIN{OFS="\t"}{if($2 < 0){$2 = 0; print}else{print}}' \
         | awk 'BEGIN{OFS="\t"}{$2=sprintf("%.0f",$2); $3=sprintf("%.0f",$3); print}' \
         | sort -k1,1 -k2,2n > ${species}_centromeric_arrays_${filt_1}_${filt_2}.bed
     else
       awk -F, -v filt_1=${filt_1} -v filt_2=${filt_2} 'BEGIN{OFS="\t"}NR>1 && $(NF) == filt_1{gsub(/"/,""); print $1, $2-1, $3, $8, $4}' ${file_000} \
         | awk 'BEGIN{OFS="\t"}{if($2 < 0){$2 = 0; print}else{print}}' \
         | awk 'BEGIN{OFS="\t"}{$2=sprintf("%.0f",$2); $3=sprintf("%.0f",$3); print}' \
         | sort -k1,1 -k2,2n > ${species}_centromeric_arrays_${filt_1}_${filt_2}.bed
     fi

    for val in 0 50 100 150 200 250; do
      awk -v filt_1=${filt_1} -v filt_2=${filt_2} -v val=${val} 'BEGIN{OFS="\t"}{total_length += $3-$2}END{print total_length, filt_1"_"filt_2, val, "cent_nt"}' ${species}_centromeric_arrays_${filt_1}_${filt_2}.bed >> ${species}_centromere_nt_prop
      awk -v filt_1=${filt_1} -v filt_2=${filt_2} -v val=${val} 'BEGIN{OFS="\t"}$3-$2 < 100000 && $3-$2 > val{total_length += $3-$2}END{print total_length, filt_1"_"filt_2, val, "gap_nt"}' ${species}_centromeric_repeats_${filt_1}_${filt_2}.bed.gaps >> ${species}_centromere_nt_prop
    done
  done < filter.var
  rm -I "${species}_centromeric_arrays_${filt_1}_${filt_2}.bed"
  done < <(ls -htl *_centromeric_arrays.csv | awk '{print $NF}')
