#!/bin/bash

set -e
set -u
set -o pipefail

  bin_path=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

  ## load lib
  source ${bin_path}/wfun_lib.sh

  # final tab must include the following cols
  # c("nt.count","filter_edge","min.gap","side","feature","TE.cls","step")

  # following loop iterates and produces the prop of gap space within centromere land
  # the output will be produced in the wd
  # however the input files must be called using abs path, just in case

  while read -r file_000; do
    # cent array bed file
    species=$(echo ${file_000##*/} | awk -F "_" '{print $1}')
    echo ${species}
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

      step="1"
      for val in 0 50 100 150 200 250; do
        side="upper"
        awk -v filt_1=${filt_1} -v filt_2=${filt_2} -v val=${val} -v side=${side} -v step=${step} 'BEGIN{OFS="\t"}{total_length += $3-$2}END{if(total_length > 0){print total_length, filt_1"_"filt_2, val, side, "cent_nt", "NA", step}else{print 0, filt_1"_"filt_2, val, side, "cent_nt", "NA", step}}' ${species}_centromeric_arrays_${filt_1}_${filt_2}.bed >> ${species}_centromere_nt_prop_combined
        awk -v filt_1=${filt_1} -v filt_2=${filt_2} -v val=${val} -v side=${side} -v step=${step} 'BEGIN{OFS="\t"}$3-$2 < 100000 && $3-$2 >= val{total_length += $3-$2}END{if(total_length > 0){print total_length, filt_1"_"filt_2, val, side,"gap_nt", "NA", step}else{print 0, filt_1"_"filt_2, val, side,"gap_nt", "NA", step}}' ${species}_centromeric_repeats_${filt_1}_${filt_2}.bed.gaps >> ${species}_centromere_nt_prop_combined
      done

      for val in 0 50 100 150 200 250; do
        side="lower"
        awk -v filt_1=${filt_1} -v filt_2=${filt_2} -v val=${val} -v side=${side} -v step=${step} 'BEGIN{OFS="\t"}{total_length += $3-$2}END{if(total_length > 0){print total_length, filt_1"_"filt_2, val, side, "cent_nt", "NA", step}else{print 0, filt_1"_"filt_2, val, side, "cent_nt", "NA", step}}' ${species}_centromeric_arrays_${filt_1}_${filt_2}.bed >> ${species}_centromere_nt_prop_combined
        awk -v filt_1=${filt_1} -v filt_2=${filt_2} -v val=${val} -v side=${side} -v step=${step} 'BEGIN{OFS="\t"}$3-$2 < 100000 && $3-$2 < val{total_length += $3-$2}END{if(total_length > 0){print total_length, filt_1"_"filt_2, val, side,"gap_nt", "NA", step}else{print 0, filt_1"_"filt_2, val, side,"gap_nt", "NA", step}}' ${species}_centromeric_repeats_${filt_1}_${filt_2}.bed.gaps >> ${species}_centromere_nt_prop_combined
      done

    done < filter.var
    rm -I "${species}_centromeric_arrays_${filt_1}_${filt_2}.bed"
  done < <(ls -htl *_centromeric_arrays.csv | awk '{print $NF}')

  # what follows is the prop of TE space within all gap land
  while read -r file_001; do
    # cent array bed file
    species=$(echo ${file_001##*/} | awk -F "_" '{print $1}')
    echo ${species}
    awk -F, 'BEGIN{OFS="\t"}{if(NR>1 && NF == 20){gsub(/"/,""); print $(NF-1)"_"$NF}else if(NR>1 && NF == 19){print $NF"_""NA"}}' ${file_001} \
      | sort \
      | uniq > filter.var

    while read -r filt; do
      while read -r file_000; do

        step="2"
        for val in 0 50 100 150 200 250; do
          side="upper"
          awk -v val=${val} 'BEGIN{OFS="\t"}NR>1 && $3-$2 < 100000 && $3-$2 >= val{print $1, $2, $3, $4, $5, $6}' ${file_000} \
            | awk '!seen[$0]++' \
            | awk -v filt=${filt} -v val=${val} -v side=${side} -v step=${step} 'BEGIN{OFS="\t"}{total_length += $3-$2}END{if(total_length > 0){print total_length, filt, val, side, "gap_nt", "NA", step}else{print 0, filt, val, side, "gap_nt", "NA", step}}' >> ${species}_centromere_nt_prop_combined

          awk -v val=${val} 'BEGIN{OFS="\t"}NR>1 && $3-$2 < 100000 && $3-$2 >= val{print $7, $8, $9, $10, $11, $12, $13}' ${file_000} \
            | awk '!seen[$0]++' \
            | awk -v filt=${filt} -v val=${val} -v side=${side} -v step=${step} 'BEGIN{OFS="\t"}{total_length += $NF}END{if(total_length > 0){print total_length, filt, val, side, "TE_nt", "NA", step}else{print 0, filt, val, side, "TE_nt", "NA", step}}' >> ${species}_centromere_nt_prop_combined
        done

        for val in 0 50 100 150 200 250; do
          side="lower"
          awk -v val=${val} 'BEGIN{OFS="\t"}NR>1 && $3-$2 < 100000 && $3-$2 < val{print $1, $2, $3, $4, $5, $6}' ${file_000} \
            | awk '!seen[$0]++' \
            | awk -v filt=${filt} -v val=${val} -v side=${side} -v step=${step} 'BEGIN{OFS="\t"}{total_length += $3-$2}END{if(total_length > 0){print total_length, filt, val, side, "gap_nt", "NA", step}else{print 0, filt, val, side, "gap_nt", "NA", step}}' >> ${species}_centromere_nt_prop_combined

          awk -v val=${val} 'BEGIN{OFS="\t"}NR>1 && $3-$2 < 100000 && $3-$2 < val{print $7, $8, $9, $10, $11, $12, $13}' ${file_000} \
            | awk '!seen[$0]++' \
            | awk -v filt=${filt} -v val=${val} -v side=${side} -v step=${step} 'BEGIN{OFS="\t"}{total_length += $NF}END{if(total_length > 0){print total_length, filt, val, side, "TE_nt", "NA", step}else{print 0, filt, val, side, "TE_nt", "NA", step}}' >> ${species}_centromere_nt_prop_combined
        done

       done < <(ls -htl ${species}_centromeric_repeats_${filt}.bed.gaps.TEanno.split.all.gaps.parsed | awk '{print $NF}')
    done < filter.var
  done < <(ls -htl *_centromeric_arrays.csv | awk '{print $NF}')

  # split by TE.type
  # prop of gaps annotated as TEs
  while read -r file_001; do
    # cent array bed file
    species=$(echo ${file_001##*/} | awk -F "_" '{print $1}')
    echo ${species}
    awk -F, 'BEGIN{OFS="\t"}{if(NR>1 && NF == 20){gsub(/"/,""); print $(NF-1)"_"$NF}else if(NR>1 && NF == 19){print $NF"_""NA"}}' ${file_001} \
      | sort \
      | uniq > filter.var

    while read -r filt; do
      while read -r file_000; do
      species=$(echo ${file_000##*/} | awk -F "_" '{print $1}')
      echo ${species}

      sort -k10,10 ${file_000} \
        | join -1 10 -2 2 -a 1 -a 2 - <(tr ' ' '\t' < TE_general_cls | sort -k2,2) \
        | awk 'NF > 2' \
        | tr ' ' '\t' > temp_file_000

        while read -r cls; do
          step="3"
          # checking the upper side of the dat
          for val in 0 50 100 150 200 250; do
            side="upper"
            awk -v val=${val} 'BEGIN{OFS="\t"}NR>1 && $3-$2 < 100000 && $3-$2 >= val{print $1, $2, $3, $4, $5, $6}' ${file_000} \
              | awk '!seen[$0]++' \
              | awk -v filt=${filt} -v val=${val} -v cls=${cls} -v side=${side} -v step=${step} 'BEGIN{OFS="\t"}{total_length += $3-$2}END{if(total_length > 0){print total_length, filt, val, side, "gap_nt", cls, step}else{print 0, filt, val, side, "gap_nt", cls, step}}' >> ${species}_centromere_nt_prop_combined

            awk -v val=${val} -v cls=${cls} 'BEGIN{OFS="\t"}$4-$3 < 100000 && $4-$3 >= val && $NF == cls{print $8, $9, $10, $11, $12, $1, $13}' temp_file_000 \
              | awk '!seen[$0]++' \
              | awk -v filt=${filt} -v val=${val} -v cls=${cls} -v side=${side} -v step=${step} 'BEGIN{OFS="\t"}{total_length += $NF}END{if(total_length > 0){print total_length, filt, val, side, "TE_nt", cls, step}else{print 0, filt, val, side, "TE_nt", cls, step}}' >> ${species}_centromere_nt_prop_combined
          done

          # checking the lower side of the data
          for val in 0 50 100 150 200 250; do
            side="lower"
            awk -v val=${val} 'BEGIN{OFS="\t"}NR>1 && $3-$2 < 100000 && $3-$2 < val{print $1, $2, $3, $4, $5, $6}' ${file_000} \
              | awk '!seen[$0]++' \
              | awk -v filt=${filt} -v val=${val} -v cls=${cls} -v side=${side} -v step=${step} 'BEGIN{OFS="\t"}{total_length += $3-$2}END{if(total_length > 0){print total_length, filt, val, side, "gap_nt", cls, step}else{print 0, filt, val, side, "gap_nt", cls, step}}' >> ${species}_centromere_nt_prop_combined

            awk -v val=${val} -v cls=${cls} 'BEGIN{OFS="\t"}$4-$3 < 100000 && $4-$3 < val && $NF == cls{print $8, $9, $10, $11, $12, $1, $13}' temp_file_000 \
              | awk '!seen[$0]++' \
              | awk -v filt=${filt} -v val=${val} -v cls=${cls} -v side=${side} -v step=${step} 'BEGIN{OFS="\t"}{total_length += $NF}END{if(total_length > 0){print total_length, filt, val, side, "TE_nt", cls, step}else{print 0, filt, val, side, "TE_nt", cls, step}}' >> ${species}_centromere_nt_prop_combined
          done
        done < <(awk '{print $1}' TE_general_cls | awk '!seen[$0]++')

        if [ -f temp_file_000 ]; then
          rm -I temp_file_000
        fi

       done < <(ls -htl ${species}_centromeric_repeats_${filt}.bed.gaps.TEanno.split.all.gaps.parsed | awk '{print $NF}')
    done < filter.var
  done < <(ls -htl *_centromeric_arrays.csv | awk '{print $NF}')

  # split by EDTA.cls
  # prop of gaps annotated as TEs
  while read -r file_001; do
    # cent array bed file
    species=$(echo ${file_001##*/} | awk -F "_" '{print $1}')
    echo ${species}
    awk -F, 'BEGIN{OFS="\t"}{if(NR>1 && NF == 20){gsub(/"/,""); print $(NF-1)"_"$NF}else if(NR>1 && NF == 19){print $NF"_""NA"}}' ${file_001} \
      | sort \
      | uniq > filter.var

    while read -r filt; do
      while read -r file_000; do
      species=$(echo ${file_000##*/} | awk -F "_" '{print $1}')
      echo ${species}

      sort -k10,10 ${file_000} \
        | join -1 10 -2 2 -a 1 -a 2 - <(tr ' ' '\t' < TE_general_cls | sort -k2,2) \
        | awk 'NF > 2' \
        | tr ' ' '\t' > temp_file_000

      awk 'NR>1{print $10}' ${file_000} \
        | sort | uniq | awk 'length($0) > 1' > cls_list

        while read -r cls; do
          step="4"
          # checking the upper side of the dat
          for val in 0 50 100 150 200 250; do
            side="upper"
            awk -v val=${val} 'BEGIN{OFS="\t"}NR>1 && $3-$2 < 100000 && $3-$2 >= val{print $1, $2, $3, $4, $5, $6}' ${file_000} \
              | awk '!seen[$0]++' \
              | awk -v filt=${filt} -v val=${val} -v cls=${cls} -v side=${side} -v step=${step} 'BEGIN{OFS="\t"}{total_length += $3-$2}END{if(total_length > 0){print total_length, filt, val, side, "gap_nt", cls, step}else{print 0, filt, val, side, "gap_nt", cls, step}}' >> ${species}_centromere_nt_prop_combined

            awk -v val=${val} -v cls=${cls} 'BEGIN{OFS="\t"}$4-$3 < 100000 && $4-$3 >= val && $1 == cls{print $8, $9, $10, $11, $12, $1, $13}' temp_file_000 \
              | awk '!seen[$0]++' \
              | awk -v filt=${filt} -v val=${val} -v cls=${cls} -v side=${side} -v step=${step} 'BEGIN{OFS="\t"}{total_length += $NF}END{if(total_length > 0){print total_length, filt, val, side, "TE_nt", cls, step}else{print 0, filt, val, side, "TE_nt", cls, step}}' >> ${species}_centromere_nt_prop_combined
          done

          # checking the lower side of the data
          for val in 0 50 100 150 200 250; do
            side="lower"
            awk -v val=${val} 'BEGIN{OFS="\t"}NR>1 && $3-$2 < 100000 && $3-$2 < val{print $1, $2, $3, $4, $5, $6}' ${file_000} \
              | awk '!seen[$0]++' \
              | awk -v filt=${filt} -v val=${val} -v cls=${cls} -v side=${side} -v step=${step} 'BEGIN{OFS="\t"}{total_length += $3-$2}END{if(total_length > 0){print total_length, filt, val, side, "gap_nt", cls, step}else{print 0, filt, val, side, "gap_nt", cls, step}}' >> ${species}_centromere_nt_prop_combined

            awk -v val=${val} -v cls=${cls} 'BEGIN{OFS="\t"}$4-$3 < 100000 && $4-$3 < val && $1 == cls{print $8, $9, $10, $11, $12, $1, $13}' temp_file_000 \
              | awk '!seen[$0]++' \
              | awk -v filt=${filt} -v val=${val} -v cls=${cls} -v side=${side} -v step=${step} 'BEGIN{OFS="\t"}{total_length += $NF}END{if(total_length > 0){print total_length, filt, val, side, "TE_nt", cls, step}else{print 0, filt, val, side, "TE_nt", cls, step}}' >> ${species}_centromere_nt_prop_combined
          done
        done < cls_list

        if [ -f temp_file_000 ]; then
          rm -I temp_file_000
          rm -I cls_list
        fi

       done < <(ls -htl ${species}_centromeric_repeats_${filt}.bed.gaps.TEanno.split.all.gaps.parsed | awk '{print $NF}')
    done < filter.var
  done < <(ls -htl *_centromeric_arrays.csv | awk '{print $NF}')
