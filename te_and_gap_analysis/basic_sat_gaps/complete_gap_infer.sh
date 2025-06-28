#!/bin/bash

set -e
set -u
set -o pipefail

  bin_path=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

  ## load lib
  source ${bin_path}/wfun_lib.sh


  while read -r file_000; do
   # cent array bed file
   species=$(echo ${file_000##*/} | awk -F "_" '{print $1}')
   echo "${file_000##*/}"

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

     if [ ! -f ${species}_centromeric_repeats.bed ]; then
     while read -r TRASH_class; do
       # cent repeats bed file
         file_001=$(ls ${species}*_repeats_filtered.csv)
         awk -F, -v TRASH_class=${TRASH_class} 'BEGIN{OFS="\t"}NR>1{gsub(/"/,""); print $1, $3, $4, $(NF-1), $5}' ${file_001} \
           | awk -v TRASH_class=${TRASH_class} '$4 == TRASH_class' \
           | awk 'BEGIN{OFS="\t"}{$2=sprintf("%.0f",$2); $3=sprintf("%.0f",$3); print}' \
           | awk 'BEGIN{OFS="\t"}{$2=$2-1; print}' \
           | awk 'BEGIN{OFS="\t"}{if($2 < 0){$2 = 0; print}else{print}}' \
           | sort -k1,1 -k2,2n >> ${species}_centromeric_repeats.bed
     done < <(cut -f4 ${species}_centromeric_arrays_${filt_1}_${filt_2}.bed | var_count - | awk 'NR==1{print $2}' | tr ';' '\n')
     fi


     if [ -f ${species}_centromeric_repeats.bed ]; then
       # intersect
       bedtools intersect -a ${species}_centromeric_arrays_${filt_1}_${filt_2}.bed -b ${species}_centromeric_repeats.bed > temp_file_000
       mv temp_file_000 ${species}_centromeric_repeats_${filt_1}_${filt_2}.bed
       # infer gaps
       Rscript gap_infer.R ${species}_centromeric_repeats_${filt_1}_${filt_2}.bed
     else
       echo "${species}_centromeric_repeats.bed" >> missing.repeats
     fi

    test_rm_ith "${species}_centromeric_repeats.bed"
    test_rm_ith "${species}_centromeric_arrays_${filt_1}_${filt_2}.bed"
    test_rm_ith "${species}_centromeric_repeats_${filt_1}_${filt_2}.bed"
   done < filter.var
  done < <(ls -htl *_centromeric_arrays.csv | awk '{print $NF}')
