#!/bin/bash

set -e
set -u
set -o pipefail


bin_path=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

## load lib
source ${bin_path}/wfun_lib.sh

# rm -I *list*
rm -I all_species_tally.txt

for species in drHedHeli1 daBalNigr1 daSheArve1 daMisOron1 daLinVulg1 ddEupPepu3 icHarAxyr1; do

  # gaps
  # gaps and overlapping TEs
  awk 'BEGIN{OFS="\t"}NR > 1 && $8 != -1{print $7, $8, $9, $6"#"$10}' ${species}.1.fa_centromeric_repeats_TRUE_FALSE.bed.gaps.TEanno.split.all.gaps.parsed \
   | sort -k1,1 -k2,2n > ${species}.1.fa_centromeric_repeats_TRUE_FALSE.bed.gaps.TEanno.split.all.gaps.parsed.bed

  # just gapss
  awk 'BEGIN{OFS="\t"}$3-$2 <= 100000{print $1, $2, $3, $6}' ${species}.1.fa_centromeric_repeats_TRUE_FALSE.bed.gaps.TEanno.split.all.gaps \
    | awk '!a[$0]++' \
    | sort -k1,1 -k2,2n > ${species}.1.fa_centromeric_repeats_TRUE_FALSE.bed.gaps.TEanno.split.all.gaps.bed

  # soloLTR
  awk 'BEGIN{OFS="\t"}NR > 1 && $(NF-1) == "IN"{print $1, $4, $5, $3"#"$NF}' ${species}.1.fa_wga_soloLTR.gff3.cent_occurrence \
    | sort -k1,1 -k2,2n > ${species}.1_wga_soloLTR.gff3.cent_occurrence.bed

  # intact elements
  sed 's/ /:/g' ${species}.1.fa.mod.EDTA.intact.fa.rexdb-*.cls.tsv \
    | awk 'BEGIN{OFS="\t"}NR > 1 && $1 ~ $2 && $7 != "none"{split($1,coord,":|\\.|-"); print coord[1], coord[2], coord[4],$1"/"$3"_"$4}' \
    | sort -k1,1 -k2,2n > ${species}.1.fa.mod.EDTA.intact.fa.rexdb-plant.cls.tsv.bed

  #  soloLTR
  #  minimum soloLTR length == 150 bp

  if [ "$(wc -l < ${species}.1_wga_soloLTR.gff3.cent_occurrence.bed)" -gt 0 ]; then

    bedtools intersect -a ${species}.1.fa_centromeric_repeats_TRUE_FALSE.bed.gaps.TEanno.split.all.gaps.bed -b ${species}.1_wga_soloLTR.gff3.cent_occurrence.bed -wao \
      | awk '$6 != -1{split($4,true_gap_length,"_|#");
                        true_solo_length = $7-$6;
                        split($4,gap_id,"#");
                        print gap_id[1], $8, int(($7-$6)/100)*100"_"int($NF/100)*100"#"$(NF-1), true_gap_length[2], true_solo_length, true_gap_length[2]-true_solo_length, $5":"$6".."$7}' \
      | sort -k5,5n \
      | awk '$(NF-2) >= 150' \
      | awk '{print $1, $2, $NF, $(NF-3), $(NF-2), $(NF-1)}' \
      | awk '!a[$0]++' > ${species}_all_in_gaps_solo_list.txt

    wc -l ${species}_all_in_gaps_solo_list.txt \
      | awk -v sp=${species} 'BEGIN{OFS="\t"}{print sp, "soloLTR", "all_in_gaps", $1}' >> all_species_tally.txt

    bedtools intersect -a ${species}.1.fa_centromeric_repeats_TRUE_FALSE.bed.gaps.TEanno.split.all.gaps.parsed.bed -b ${species}.1_wga_soloLTR.gff3.cent_occurrence.bed -wao \
      | awk '$6 != -1{split($4,true_gap_length,"_|#");
                    true_solo_length = $7-$6;
                    split($4,gap_id,"#");
                    print gap_id[1], $8, int(($7-$6)/100)*100"_"int($NF/100)*100"#"$(NF-1), true_gap_length[2], true_solo_length, true_gap_length[2]-true_solo_length, $5":"$6".."$7}' \
      | sort -k5,5n \
      | awk '!a[$NF]++' \
      | awk '$(NF-2) >= 150' \
      | awk '{print $1, $2, $NF, $(NF-3), $(NF-2), $(NF-1)}' \
      | awk '!a[$0]++' > ${species}_all_solo_list.txt

    wc -l ${species}_all_solo_list.txt \
      | awk -v sp=${species} 'BEGIN{OFS="\t"}{print sp, "soloLTR", "all_in_overlap", $1}' >> all_species_tally.txt

    bedtools intersect -a ${species}.1.fa_centromeric_repeats_TRUE_FALSE.bed.gaps.TEanno.split.all.gaps.parsed.bed -b ${species}.1_wga_soloLTR.gff3.cent_occurrence.bed -wao \
      | awk '$6 != -1{split($4,true_gap_length,"_|#");
                  true_solo_length = $7-$6;
                  split($4,gap_id,"#");
                  print gap_id[1], $8, int(($7-$6)/100)*100"_"int($NF/100)*100"#"$(NF-1), true_gap_length[2], true_solo_length, true_gap_length[2]-true_solo_length, $5":"$6".."$7}' \
      | sort -k5,5n \
      | awk '!a[$NF]++' \
      | awk '$(NF-2) >= 150' \
      | awk '$(NF-1) >= -100 && $(NF-1) <= 100' \
      | awk '{print $1, $2, $NF, $(NF-3), $(NF-2), $(NF-1)}' >> ${species}_solo_list_100_bp.txt

    wc -l ${species}_solo_list_100_bp.txt \
      | awk -v sp=${species} 'BEGIN{OFS="\t"}{print sp, "soloLTR", "all_in_overlap_100", $1}' >> all_species_tally.txt
  fi

  # intact TEs

  bedtools intersect -a ${species}.1.fa_centromeric_repeats_TRUE_FALSE.bed.gaps.TEanno.split.all.gaps.bed -b ${species}.1.fa.mod.EDTA.intact.fa.rexdb-plant.cls.tsv.bed -wao \
    | awk '$6 != -1 && $NF > 1{$NF = int($NF/100)*100;
                                       gap_length = int(($3-$2)/100)*100; 
                                       split($4,true_gap_length,"_|#");
                                       true_edta_length = $7-$6;
                                       split($(NF-1),fam,"/"); 
                                       split($4,gap_id,"#");
                                       print gap_id[1], $NF"_"fam[2], fam[2], true_gap_length[2], true_edta_length, true_gap_length[2]-true_edta_length, $5":"$6".."$7}' \
    | sort -k3,3n \
    | awk '{print $1, $3, $NF, $(NF-3), $(NF-2), $(NF-1)}'\
    | awk '!a[$0]++' > ${species}_all_in_gaps_intact_list.txt

  wc -l ${species}_all_in_gaps_intact_list.txt \
    | awk -v sp=${species} 'BEGIN{OFS="\t"}{print sp, "intact", "all_in_gaps", $1}' >> all_species_tally.txt


  bedtools intersect -a ${species}.1.fa_centromeric_repeats_TRUE_FALSE.bed.gaps.TEanno.split.all.gaps.parsed.bed -b ${species}.1.fa.mod.EDTA.intact.fa.rexdb-plant.cls.tsv.bed -wao \
    | awk '$6 != -1 && $NF > 1{$NF = int($NF/100)*100;
                                     gap_length = int(($3-$2)/100)*100; 
                                     split($4,true_gap_length,"_|#");
                                     true_edta_length = $7-$6;
                                     split($(NF-1),fam,"/"); 
                                     split($4,gap_id,"#");
                                     print gap_id[1], $NF"_"fam[2], fam[2], true_gap_length[2], true_edta_length, true_gap_length[2]-true_edta_length, $5":"$6".."$7}' \
    | sort -k3,3n \
    | awk '!a[$NF]++' \
    | awk '{print $1, $3, $NF, $(NF-3), $(NF-2), $(NF-1)}'\
    | awk '!a[$0]++' > ${species}_all_intact_list.txt

  wc -l ${species}_all_intact_list.txt \
    | awk -v sp=${species} 'BEGIN{OFS="\t"}{print sp, "intact", "all_in_overlap", $1}' >> all_species_tally.txt

  # minimum intact length == 2000 bp

  bedtools intersect -a ${species}.1.fa_centromeric_repeats_TRUE_FALSE.bed.gaps.TEanno.split.all.gaps.parsed.bed -b ${species}.1.fa.mod.EDTA.intact.fa.rexdb-plant.cls.tsv.bed -wao \
    | awk '$6 != -1 && $NF > 1{$NF = int($NF/100)*100;
                                   gap_length = int(($3-$2)/100)*100; 
                                   split($4,true_gap_length,"_|#");
                                   true_edta_length = $7-$6;
                                   split($(NF-1),fam,"/"); 
                                   split($4,gap_id,"#");
                                   print gap_id[1], $NF"_"fam[2], fam[2], true_gap_length[2], true_edta_length, true_gap_length[2]-true_edta_length, $5":"$6".."$7}' \
    | awk '!a[$NF]++' \
    | awk '$(NF-1) <= 200 && $(NF-1) >= -200' \
    | sort -k3,3n \
    | awk '{print $1, $3, $NF, $(NF-3), $(NF-2), $(NF-1)}' \
    | awk '$(NF-2) >= 2000' \
    | awk '!a[$1]++' > ${species}_intact_list_200_bp.txt

  wc -l ${species}_intact_list_200_bp.txt \
    | awk -v sp=${species} 'BEGIN{OFS="\t"}{print sp, "intact", "all_in_overlap_200", $1}' >> all_species_tally.txt

done
