# Analysis of 193 genomes and their centromeres
## Requirements
Genome assembly data (example Arabidopsis thaliana chromosome 1 coordinates 1 : 20,000,000 attached in _/test_data/ddAraThal4.1_chr1_20Mb.fa_)
Installed 3rd party software according to its documentation:
* R version 4.1.3 or newer with libraries and their dependancies that include (each script lists required libraries):
	* seqinr
	* msa
	* Biostrings
	* dplyr
	* Matrix
	* pheatmap
	* GenomicRanges
* TRASH2 (https://github.com/vlothec/TRASH_2) which also requires:
	* MAFFT 7.526
	* HMMER 3.4
* EDTA
* TEsorter
* Helixer 0.3.4
## Putative centromeric repeat identification
1. Run annotation software:
	* TRASH2
	* EDTA
	* TEsorter
2. Filter annotations with scripts found in /repeat_te_gene_annotations_parsing, in order they appear. Modify paths at the beginning of the script to match the outputs from the step above
3. Calculate initial predictor scores for individual repeat families with _/centromere_identification/6_find_centromeric_repeat.R_ script
4. Create a genomic landscape plot to visualise lcoations and predictor scores of top scoring families and decide on the putative centromeric repeat with _/centromere_identification/10_global_plot.R_
## Centromere coordinates, gaps and genomic landscape replotting (Supplemental Data 1)
Most downstream analyses require a metadata .csv table with information on the previously identified centromeric families and classification of all genomes into a A. satellite based monocentric, B. satellite based holocentric, C. transposon based monocentric and D. unknown
Another metadata file containing information (source genome name, repeat TRASH2 name and repeat custom name) on each inentified centromeric family is also required.
* _/centromere_identification/11_array_cen_pericen_coordinates.R_ Identifies centromeric and centromere-proximal regions for genomes with identified centromeric repeats
* _/centromere_identification/12_find_gaps.R_ Identifies gaps in the centromeric arrays
* _/centromere_identification/10.1_global_plot_plus_cen_coordinates.R_ replots the genomic landscape plot with additional information from the _/11_array_cen_pericen_coordinates.R_ script and highlighting identified centromeric repeat
## HOR Analysis
* _/hor_scripts/make_hor_scripts.R_ Extracts which chromosomes contain centromeric satellites and what are their names in order to create SLURM submission commands that will analyse them individually.
* _/hor_scrips/9_HOR_periods.R_ Creates visualisations of identified HORs in order to assess runs completeness and inform on the HOR content of individual genomes
* _/hor_scripts/0.01_HOR_score.R_ Calculates a HOR score for each individual repeat that had them identified 
## Downstream Analysis
/main_downstream_analysis Contains scripts used for remaining analyses, with script names describing their purpose and comments section at the beginning of many scripts adding details on their functionality and algorithms when needed.
While numbered, these scripts do not have prerequisites of finishing previous scripts in this directory in order to run them. 
## Summary of genomes and centromeric satellites (Supplemental Table 1)
All scripts mentioned in the **Putative centromere repeat identification**, **Centromere coordinates...** and **HOR anlaysis** sections above need to be done. 3 additional scripts also need to be run and finished for the full table to be constructed, those calculate the most computationally expensive steps of the summary:
* _/main_summary/13.1_table_S1_HOR_stats.R_ Reads in and summarises results of the HOR analyses
* _/main_summary/13.2_table_S1_GC.R_ Calculates GC values for genomes, chromosomes and genomic intervals (repeats, transposons, their subsets etc.)
* _/main_downstream_analysis/37_similarity_values_within_between_chr.r_ Calculates similarity values within and between chromosomes for centromeric repeats
_/main_summary/13_table_S1_whole_summary.R_ Generates the full table for all analysed genomes
## dN/dS analysis
/dNdS/analyze_dnds.ipnyb describes and visualises the CENH3 protein site dN/dS Analysis
## Transposable element and gap Analysis
_/te_and_gap_analysis/_ directory contains scripts for transposon reannotation and rescue steps and those required for individual figures creation in subdirectories:
* _/te_and_gap_analysis/Figure_3/_
* _/te_and_gap_analysis/Figure_4/_
* _/te_and_gap_analysis/supp_data/_

## Other
* _/FISH_oligo-pools/_ contains pools of oligonucleotides used in FISH experiments
* [Figshare S1 figure](https://figshare.com/articles/dataset/193centromeres_S1/29436083) contains plots made with _/centromere_identification/10.1_global_plot_plus_cen_coordinates.R_ script that make the Supplemental Data 1 figure
* [Figshare additional data](https://figshare.com/articles/dataset/193centromeres/29412917) contains all main data created by the scripts in this repository

