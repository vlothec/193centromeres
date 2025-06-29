directrory	type	rep_name	obs	main_input_files
figure_3	script	figure_3_panel_A.R	bar plots showing proportion of centromeric features	data from cent_prop_combined.sh (*_centromere_nt_prop_combined)
figure_3	script	figure_3_panel_B.R	bar plots showing centromere length in Mb	data from cent_prop_combined.sh (*_centromere_nt_prop_combined)
figure_3	script	figure_3_panel_C.R	bar plots showing proportion of intact/fragmented transposon annotation	data from structural_homology_ratio.R (*_structural_homology_ratio)
figure_3	script	figure_3_panel_D.R	dot plots with LTR identity values	data from te_pident_stats_intact.R (*_tally_filt_TRUE_FALSE_pident)
figure_3	script	figure_3_panel_E.R	dot plots showing proportion of transposon classes	data from te_features_tally.R (*.fa_tally_filt_TRUE_FALSE)
figure_3	script	stats_figure_3_panel_C.R	stats (two sample independent t-test) for intact/fragmented transposon annotation	data from structural_homology_ratio.R (*_structural_homology_ratio)
figure_3	script	stats_figure_3_panel_D.R	stats (two sample independent t-test) for LTR identity values	data from te_pident_stats_intact.R (*.fa_tally_filt_TRUE_FALSE)
figure_3	script	figure_3_panel_H.R	gaps in satellite arrays	data from gap_infer_prop_all_gaps.R (*all.gaps.parsed)
figure_4	script	figure_4_panel_A.R	ATHILA tree	.
figure_4	script	figure_4_panel_B.R	drGeuUrba1 centromere space	.
figure_4	txt	data_for_figure_4_panel_a.tsv	data for ATHILA tree	.
figure_4	txt	plant.Gypsy_Athila.retree2.max100.mafft.fasttree	data for ATHILA tree	.
figure_4/panle_C	script	figure_4_panel_C.R	boxplot in figure_4_panel_C	.
figure_4/panle_C	txt	nest_1_SUPER_15_18115606_18180343.txt	coordinates of the sequences used for the dotplot shown in figure_4_panel_C	.
figure_4/panle_C/dotplot_seq	fasta	nest_1_strata_1.fa	dotplot sequences shown in figure_4_panel_C	.
figure_4/panle_C/dotplot_seq	fasta	nest_1_strata_2.fa	dotplot sequences shown in figure_4_panel_C	.
figure_4/panle_C/dotplot_seq	fasta	nest_1_strata_3.fa	dotplot sequences shown in figure_4_panel_C	.
figure_4/panle_C/dotplot_seq	fasta	nest_1_original_extended_region.fa	dotplot sequences shown in figure_4_panel_C	.
figure_4/panle_C/p_ident	txt	drGeuUrba1.1_stitched_strata_1.anno	extended data for TEs considered in figure_4_panel_C boxplot strata1	.
figure_4/panle_C/p_ident	txt	drGeuUrba1.1_stitched_strata_2.anno	extended data for TEs considered in figure_4_panel_C boxplot strata2	.
figure_4/panle_C/p_ident	txt	drGeuUrba1.1_stitched_strata_3.anno	extended data for TEs considered in figure_4_panel_C boxplot strata3	.
supp_data	script	figure_s7.R	bar plots showing centromere content in satellite-based species	data from cent_prop_combined.sh (*_centromere_nt_prop_combined)
supp_data	script	figure_s9.R	pie charts showing centrophilic LTRRT in plants	data from sat_invasion_tally.R (*.intact.fa.rexdb-plant.cls.tsv.tally)
supp_data	script	figure_s11_panel_B.R	frequency of gaps along a scaled chromosome	.
supp_data	script	figure_s11_panel_C.R	gaps in satellite arrays ordered by physical position	data from gap_infer_prop_all_gaps.R (*.split.all.gaps.parsed)
supp_data	script	figure_s12.R	multipage pdf with ATHILA tree	.
supp_data	script	figure_s13_panel_A.R	drGeuUrba1 centromere space by chromosome	data from drGeuUrba1_cent_prop_step_3.R and drGeuUrba1_cent_prop_step_5.R
supp_data	script	table_s2.R	transposon centromere content in satellite-based species	data from cent_prop_combined.sh (*_centromere_nt_prop_combined)
basic_te_anno	script	wfun_lib.sh	basic functions	.
basic_te_anno	script	TEsorter_with_EDTA_intact.sh	TEsorter classification	.
basic_te_anno	script	TEsorter_with_filtered_rescued_files.sh	TEsorter classification	.
basic_te_anno	script	ann_ltrrt_structure.sh	TEsorter classification	.
basic_te_anno	script	TEsorter_diagnostic.sh	TEsorter classification	.
basic_te_anno	script	concat_hmm_domains.sh	phylogenetic tree	.
basic_te_anno	script	subset_concat_hmm_domains.sh	phylogenetic tree	.
basic_te_anno	script	mafft_plus_fasttree_alone.sh	phylogenetic tree	.
basic_te_anno/rescue_step	script	dtol_rescue_step.sh	reassign repeat_region TEs when possible	.
basic_te_anno/rescue_step	script	repeat_region_rescue	reassign repeat_region TEs when possible	.
basic_te_anno/rescue_step	txt	step1.map	reassign repeat_region TEs when possible	.
basic_te_anno/rescue_step	txt	step2.map	reassign repeat_region TEs when possible	.
basic_te_anno/rescue_step	txt	step3.map	reassign repeat_region TEs when possible	.
basic_te_anno	txt	basic_te_readme.txt	information for basic_te_anno scripts	.
basic_sat_gaps	script	centromeric_arrays_tsv_files.sh	get relevant fields from *_centromeric_arrays.csv file	*_centromeric_arrays.csv
basic_sat_gaps	script	complete_gap_infer.sh	get gap coordinates 	*_repeats_filtered.csv and *_centromeric_arrays.csv 
basic_sat_gaps	script	gap_infer.R	get gap coordinates	.
basic_sat_gaps	script	gap_isec_with_split_anno.sh	intersect gap coordinates and the split_filtered_reassigned TE annotation 	data from complete_gap_infer.sh (*.bed.gaps)
basic_sat_gaps	script	gap_infer_prop_all_gaps_loop.sh	get proportion of mapped TEs per gap 	data from gap_isec_with_split_anno.sh (*.split.all.gaps)
basic_sat_gaps	script	gap_infer_prop_all_gaps.R	get proportion of mapped TEs per gap	.
basic_sat_gaps	script	cent_prop_combined.sh	get proportion of different features at different levels 	data from gap_isec_with_split_anno.sh (*.split.all.gaps) and *_centromeric_arrays.csv and TE_general_cls
basic_sat_gaps	script	prop_unaccounted_cent.sh	get proportion of different features at different levels	.
basic_sat_gaps	script	array_boundaries_annotation_loop.sh	annotate split_filtered_reassigned and intact TE files with array boundaries 	data from centromeric_arrays_tsv_files.sh (*all_centromeric_arrays.tsv) and either *.EDTA.intact.gff3 or *_edta_filtered.csv.reassigned
basic_sat_gaps	script	array_boundaries_annotation.R	annotate split_filtered_reassigned and intact TE files with array boundaries	.
basic_sat_gaps	script	structural_homology_ratio_loop.sh	structural/homology data for figure_3_panel_C 	data from array_boundaries_annotation_loop.sh (*.fa_edta_filtered.csv.reassigned.array_boundaries)
basic_sat_gaps	script	structural_homology_ratio.R	structural/homology data for figure_3_panel_C	.
basic_sat_gaps	script	te_pident_stats_loop.sh	LTR identity data for figure_3_panel_D 	data from array_boundaries_annotation_loop.sh (*.EDTA.intact.gff3.array_boundaries)
basic_sat_gaps	script	te_pident_stats_intact.R	LTR identity data for figure_3_panel_D	.
basic_sat_gaps	script	te_features_tally_loop.sh	TE counts for figure_3_panel_E 	data from array_boundaries_annotation_loop.sh (*.fa_edta_filtered.csv.reassigned.array_boundaries)
basic_sat_gaps	script	te_features_tally.R	TE counts for figure_3_panel_E	.
basic_sat_gaps	script	soloLTRseeker_run_sat_species.sh	soloLTR and intact data for the 6 species shown in figure_3_panel_h 	*.EDTA.intact.gff3 and *.fa
basic_sat_gaps	script	soloLTR_cent_occ_loop.sh	soloLTR and intact data for the 6 species shown in figure_3_panel_h 	data from soloLTRseeker_run_sat_species.sh (*_wga_soloLTR.gff3)
basic_sat_gaps	script	soloLTR_cent_occ.R	soloLTR and intact data for the 6 species shown in figure_3_panel_h 	data from gap_infer_prop_all_gaps.R (*.split.all.gaps.parsed)
basic_sat_gaps	script	all_intact_and_solo_gaps.sh	soloLTR and intact data for the 6 species shown in figure_3_panel_h 	"ata from TEsorter_with_EDTA_intact.sh (*.EDTA.intact.fa.rexdb-plant.cls.tsv), gap_isec_with_split_anno.sh (*.split.all.gap), gap_infer_prop_all_gaps.R (*.split.all.gaps.parsed) and soloLTR_cent_occ_loop.sh (*_wga_soloLTR.gff3.cent_occurrence)"
basic_sat_gaps	script	all_fam_tally_loop.sh	soloLTR and intact data for the 6 species shown in figure_3_panel_h	.
basic_sat_gaps	script	all_fam_tally.R	soloLTR and intact data for the 6 species shown in figure_3_panel_h 	data from all_intact_and_solo_gaps.sh (*_lists*)
basic_sat_gaps	script	sat_invasion_tally_loop.sh	LTRRT lineage distribution in plant species with satellite-based centromeres 	data from TEsorter_with_EDTA_intact.sh (*.EDTA.intact.fa.rexdb-plant.cls.tsv) and array_boundaries_annotation_loop.sh (*.EDTA.intact.gff3.array_boundaries)
basic_sat_gaps	script	sat_invasion_tally.R	LTRRT lineage distribution in plant species with satellite-based centromeres 	.
additional_info	txt	drAilAlti1_cent.coord	manually determined centromere coordinates for transposon-based species	.
additional_info	txt	drGeuUrba1_cent.coord	manually determined centromere coordinates for transposon-based species	.
additional_info	txt	drMalDome5_cent.coord	manually determined centromere coordinates for transposon-based species	.
additional_info	txt	drMalSylv7_cent.coord	manually determined centromere coordinates for transposon-based species	.
additional_info	txt	dtol_plant_attributes.txt	relevant metadata of dtol species	.
additional_info	txt	TE_general_cls	reference table for TE groups used in some sections	.
additional_info	txt	TE_general_cls_extended	reference table for TE groups used in some sections	.
additional_info	fasta	sce_outliers.concat_dom.faa	outlier used in figure_4_panel_A	.
additional_info	txt	proper_chr.txt	list of chr considered	.