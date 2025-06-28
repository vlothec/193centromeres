#!/bin/bash

set -e
set -u
set -o pipefail

  bin_path=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

  ## load lib
  source ${bin_path}/wfun_lib.sh

  min=${@:$OPTIND:1}
  step=${@:$OPTIND+1:1}

  mkdir temp_${min}_${step}

  while read -r dir_001; do

    file_001=$(ls -htl ${dir_001}/* | awk -F/ '{print $NF}' | grep -P "csv$")

    tr -d '\\"' < ${dir_001}/${file_001} \
      | awk -F, 'BEGIN{OFS="\t"}NR>1{print $2, $3, $20, $5, $6, $7, $8, $9, "ID="$10";Name="$11";Classification="$12";Sequence_ontology="$13";Identity="$14";Method="$15";TSD="$16";TIR="$17";motif="$18";tsd="$19";oldV3="$20";overlapping_bp="$21";width="$22";overlapping_percentage="$23}' > temp_${min}_${step}/temp_file_001_${min}.gff3

    cd temp_${min}_${step}
    ./repeat_region_rescue temp_file_001_${min}.gff3

    sed "s/temp_file_001_${min}/${file_001%%_*}/g" temp_file_001_${min}.gff3.reassigned.sum > ${file_001}.reassigned.sum
    awk 'BEGIN{OFS=","}{split($9,atr,";|="); print NR,$1,$2,$3,$4,$5,$6,$7,$8,atr[2],atr[4],atr[6],atr[8],atr[10],atr[12],atr[14],atr[16],atr[18],atr[20],atr[22],atr[24],atr[26],atr[28],$NF}' temp_file_001_${min}.gff3.reassigned > ${file_001}.reassigned

    test_rm_ith "temp_file_001_${min}.gff3"
    test_rm_ith "temp_file_001_${min}.gff3.reassigned"
    test_rm_ith "temp_file_001_${min}.gff3.reassigned.sum"

    cd ../
  done < <(ls -htl *filtered.csv | awk -F/ '$0 ~ /^\//{print "edta/"$NF}' | awk -v min=${min} -v step=${step} 'NR >= min && NR < min+step' | tr -d ":")
