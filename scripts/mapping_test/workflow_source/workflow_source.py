#!/bin/env python3
from gwf import Workflow
from gwf.executors import Conda
from gwf.workflow import collect
import os, yaml, glob, sys
from workflow_templates import *

def mapping_resequencing_data_population_genetics_workflow(config_file: str = glob.glob('*config.y*ml')[0]):
	"""
	Workflow: Align resequencing data to reference genome and do basic filtering.
	
	:param str config_file:
		Configuration file containing pre-defined set of variables
	"""
	# --------------------------------------------------
	#                  Configuration
	# --------------------------------------------------
	
	CONFIG = yaml.safe_load(open(config_file))
	ACCOUNT: str = CONFIG['account']
	TAXONOMY: str = CONFIG['taxonomic_group'].lower()
	SPECIES_NAME: str = CONFIG['species_name']
	REFERENCE_GENOME: str = CONFIG['reference_genome_path']
	WORK_DIR: str = CONFIG['working_directory_path'][:len(CONFIG['working_directory_path']) - 1] if CONFIG['working_directory_path'].endswith('/') else CONFIG['working_directory_path']
	OUTPUT_DIR: str = CONFIG['output_directory_path']
	SAMPLE_LISTS: list = CONFIG['sample_lists']
	ADAPTERREMOVAL_SETTINGS: dict = CONFIG['adapterremoval_settings']
	AR_MIN_QUAL: int | None = ADAPTERREMOVAL_SETTINGS['min_quality'] if ADAPTERREMOVAL_SETTINGS['min_quality'] else 25
	AR_MIN_LENGTH: int | None = ADAPTERREMOVAL_SETTINGS['min_length'] if ADAPTERREMOVAL_SETTINGS['min_length'] else 20
	AR_SEQUENCE1: str | None = ADAPTERREMOVAL_SETTINGS['adaptersequence1'] if ADAPTERREMOVAL_SETTINGS['adaptersequence1'] else 'AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA'
	AR_SEQUENCE2: str | None = ADAPTERREMOVAL_SETTINGS['adaptersequence2'] if ADAPTERREMOVAL_SETTINGS['adaptersequence2'] else 'AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG'
	SAMTOOLS_FILTER_SETTINGS: dict = CONFIG['filter_settings']
	STF_EXCLUDE: int | None = SAMTOOLS_FILTER_SETTINGS['flags_excluded'] if SAMTOOLS_FILTER_SETTINGS['flags_excluded'] else 3844
	STF_REQUIRED: int | None = SAMTOOLS_FILTER_SETTINGS['flags_required'] if SAMTOOLS_FILTER_SETTINGS['flags_required'] else None
	STF_MIN_MQ: int | None = SAMTOOLS_FILTER_SETTINGS['min_mq'] if SAMTOOLS_FILTER_SETTINGS['min_mq'] else 20
	
	# --------------------------------------------------
	#                  Workflow
	# --------------------------------------------------
	
	gwf = Workflow(
		defaults={'account': ACCOUNT}
	)
	
	top_dir = f'{WORK_DIR}/{TAXONOMY.replace(" ", "_")}/{SPECIES_NAME.replace(" ", "_")}/mapping'
	top_out = f'{OUTPUT_DIR}/alignments/{TAXONOMY.replace(" ", "_")}/{SPECIES_NAME.replace(" ", "_")}'
	
	index = gwf.target_from_template(
		name=f'{SPECIES_NAME.replace(" ", "_")}_index_reference',
		template=index_reference_genome(
			reference_genome_file=REFERENCE_GENOME,
			output_directory=top_dir)
	)

	between_group_qualimap = []

	for GROUP in SAMPLE_LISTS:
		if not GROUP['sample_folder_list']:
			continue
		GROUP_NAME: str = GROUP['group_name'].lower()
		SAMPLE_FOLDER_LIST: list = GROUP['sample_folder_list']
		sample_list = [get_sample_data(path, data_type=2) for path in SAMPLE_FOLDER_LIST]
		within_group_qualimap = []
		preFilterAlignmentsMerge = []
		postFilterAlignmentsMerge = []
		preFilterAlignmentsNoMerge = []
		postFilterAlignmentsNoMerge = []
		preFilterAlignmentsNoCollapse = []
		postFilterAlignmentsNoCollapse = []
		for sample in sample_list:	
			adapterremoval = gwf.target_from_template(
				name=f'{sample["sample_name"].replace("-", "_")}_adapterremoval',
				template=adapterremoval_pairedend(
					sample_name=sample['sample_name'],
					read1_files=sample['read1_files'],
					read2_files=sample['read2_files'],
					output_directory=f'{top_dir}/{GROUP_NAME}',
					min_quality=AR_MIN_QUAL,
					min_length=AR_MIN_LENGTH,
					adapter1=AR_SEQUENCE1,
					adapter2=AR_SEQUENCE2
				)
			)

			adapterremoval_no_collapse = gwf.target_from_template(
				name=f'{sample["sample_name"].replace("-", "_")}_adapterremoval_no_collapse',
				template=adapterremoval_pairedend_no_collapse(
					sample_name=sample['sample_name'],
					read1_files=sample['read1_files'],
					read2_files=sample['read2_files'],
					output_directory=f'{top_dir}/{GROUP_NAME}',
					min_quality=AR_MIN_QUAL,
					min_length=AR_MIN_LENGTH,
					adapter1=AR_SEQUENCE1,
					adapter2=AR_SEQUENCE2
				)
			)

			for i in range(0,3):
				if i == 0:
					align_paired0 = gwf.target_from_template(
						name=f'{sample["sample_name"].replace("-", "_")}_paired_alignment0',
						template=alignment_pairedend(
							read1_file=adapterremoval.outputs['pair1'],
							read2_file=adapterremoval.outputs['pair2'],
							reference_genome_file=index.outputs['symlink'],
							sample_name=sample['sample_name'],
							output_directory=f'{top_dir}/{GROUP_NAME}/merge'
						)
					)

					align_collapsed0 = gwf.target_from_template(
						name=f'{sample["sample_name"].replace("-", "_")}_collapsed_alignment0',
						template=alignment_collapsed(
							collapsed_read_files=adapterremoval.outputs['collapsed'],
							reference_genome_file=index.outputs['symlink'],
							sample_name=sample['sample_name'],
							output_directory=f'{top_dir}/{GROUP_NAME}/merge'
						)
					)

					merge0 = gwf.target_from_template(
						name=f'{sample["sample_name"].replace("-", "_")}_merge_alignments0',
						template=merge_alignments(
							alignment_files=[align_paired0.outputs['alignment'],
											 align_collapsed0.outputs['alignment']],
							sample_name=sample['sample_name'],
							output_directory=f'{top_dir}/{GROUP_NAME}/merge'
						)
					)

					mark_duplicates0 = gwf.target_from_template(
						name=f'{sample["sample_name"].replace("-", "_")}_markdup0',
						template=mark_duplicates_samtools(
							alignment_file=merge0.outputs['merged'],
							sample_name=sample['sample_name'],
							output_directory=f'{top_dir}/{GROUP_NAME}/merge'
						)
					)

					preFilterAlignmentsMerge.append(mark_duplicates0.outputs['markdup'])

					insertSizePre0 = gwf.target_from_template(
						name=f'{sample["sample_name"].replace("-", "_")}_pre_filter_insert_size0',
						template=insert_size(
							bamFile=mark_duplicates0.outputs['markdup'],
							outputName=f'{sample["sample_name"]}.preFilter',
							outputDirectory=f'{top_dir}/{GROUP_NAME}/merge'
						)
					)

					depthDistributionPre0 = gwf.target_from_template(
						name=f'{sample["sample_name"].replace("-", "_")}_pre_filter_depth_distribution0',
						template=depth_distribution(
							bamFile=mark_duplicates0.outputs['markdup'],
							lowerThreshold=50,
							outputDirectory=f'{top_dir}/{GROUP_NAME}/merge/depth',
							sampleName=f'{sample["sample_name"]}.preFilter.depthDist'
						)
					)

					depthAlongReferencePre0 = gwf.target_from_template(
						name=f'{sample["sample_name"].replace("-", "_")}_pre_filter_depth_along_reference0',
						template=depth_along_reference(
							bamFile=mark_duplicates0.outputs['markdup'],
							outputDirectory=f'{top_dir}/{GROUP_NAME}/merge/depth',
							sampleName=f'{sample["sample_name"]}.preFilter.depthAlong'
						)
					)

					insertDistributionPre0 = gwf.target_from_template(
						name=f'{sample["sample_name"].replace("-", "_")}_pre_filter_insert_distribution0',
						template=insert_size_distribution(
							insertSizeFile=insertSizePre0.outputs['insert'],
							outputDirectory=f'{top_dir}/{GROUP_NAME}/merge/insert_size',
							sampleName=f'{sample["sample_name"]}.preFilter.insertDist'
						)
					)

					unmapped0 = gwf.target_from_template(
						name=f'{sample["sample_name"].replace("-", "_")}_extract_unmapped_reads0',
						template=extract_unmapped_reads(
							alignment_file=mark_duplicates0.outputs['markdup'],
							sample_name=sample['sample_name'],
							output_directory=f'{top_out}/{GROUP_NAME}/merge' if OUTPUT_DIR else f'{top_dir}/{GROUP_NAME}/merge/unmapped'
						)
					)

					pre_filter_stats0 = gwf.target_from_template(
						name=f'{sample["sample_name"].replace("-", "_")}_pre_filter_stats0',
						template=samtools_stats(
							alignment_file=mark_duplicates0.outputs['markdup'],
							output_directory=f'{top_out}/{GROUP_NAME}/merge/{sample['sample_name']}/qc/pre_filtering_stats' if OUTPUT_DIR else f'{top_dir}/{GROUP_NAME}/merge/pre_filtering_stats/{sample["sample_name"]}'
						)
					)

					filter_alignment0 = gwf.target_from_template(
						name=f'{sample["sample_name"].replace("-", "_")}_filtering0',
						template=samtools_filter(
							alignment_file=mark_duplicates0.outputs['markdup'],
							sample_name=sample['sample_name'],
							output_directory=f'{top_out}/{GROUP_NAME}/merge' if OUTPUT_DIR else f'{top_dir}/{GROUP_NAME}/merge/filtered_alignment',
							flags_excluded=STF_EXCLUDE,
							flags_required=STF_REQUIRED,
							min_mq=STF_MIN_MQ
						)
					)

					insertSizePost0 = gwf.target_from_template(
						name=f'{sample["sample_name"].replace("-", "_")}_post_filter_insert_size0',
						template=insert_size(
							bamFile=filter_alignment0.outputs['filtered'],
							outputName=f'{sample["sample_name"]}.postFilter',
							outputDirectory=f'{top_dir}/{GROUP_NAME}/merge'
						)
					)

					depthDistributionPost0 = gwf.target_from_template(
						name=f'{sample["sample_name"].replace("-", "_")}_post_filter_depth_distribution0',
						template=depth_distribution(
							bamFile=filter_alignment0.outputs['filtered'],
							lowerThreshold=50,
							outputDirectory=f'{top_dir}/{GROUP_NAME}/merge/depth',
							sampleName=f'{sample["sample_name"]}.postFilter.depthDist'
						)
					)

					depthAlongReferencePost0 = gwf.target_from_template(
						name=f'{sample["sample_name"].replace("-", "_")}_post_filter_depth_along_reference0',
						template=depth_along_reference(
							bamFile=filter_alignment0.outputs['filtered'],
							outputDirectory=f'{top_dir}/{GROUP_NAME}/merge/depth',
							sampleName=f'{sample["sample_name"]}.postFilter.depthAlong'
						)
					)

					insertDistributionPost0 = gwf.target_from_template(
						name=f'{sample["sample_name"].replace("-", "_")}_post_filter_insert_distribution0',
						template=insert_size_distribution(
							insertSizeFile=insertSizePost0.outputs['insert'],
							outputDirectory=f'{top_dir}/{GROUP_NAME}/merge/insert_size',
							sampleName=f'{sample["sample_name"]}.postFilter.insertDist'
						)
					)

					postFilterAlignmentsMerge.append(filter_alignment0.outputs['filtered'])

					post_filter_stats0 = gwf.target_from_template(
						name=f'{sample["sample_name"].replace("-", "_")}_post_filter_stats0',
						template=samtools_stats(
							alignment_file=filter_alignment0.outputs['filtered'],
							output_directory=f'{top_out}/{GROUP_NAME}/merge/{sample['sample_name']}/qc/pre_filtering_stats' if OUTPUT_DIR else f'{top_dir}/{GROUP_NAME}/merge/post_filtering_stats/{sample["sample_name"]}'
						)
					)

					qualimap0 = gwf.target_from_template(
						name=f'{sample["sample_name"].replace("-", "_")}_qualimap0',
						template=qc_qualimap(
							alignment_file=filter_alignment0.outputs['filtered'],
							output_directory=f'{top_out}/{GROUP_NAME}/merge/{sample['sample_name']}/qc/bamqc' if OUTPUT_DIR else f'{top_dir}/bamqc/individual/{GROUP_NAME}/merge/{sample["sample_name"]}'
						)
					)

				elif i == 1:
					align_paired1 = gwf.target_from_template(
						name=f'{sample["sample_name"].replace("-", "_")}_paired_alignment1',
						template=alignment_pairedend(
							read1_file=adapterremoval.outputs['pair1'],
							read2_file=adapterremoval.outputs['pair2'],
							reference_genome_file=index.outputs['symlink'],
							sample_name=sample['sample_name'],
							output_directory=f'{top_dir}/{GROUP_NAME}/no_merge'
						)
					)

					mark_duplicates1 = gwf.target_from_template(
						name=f'{sample["sample_name"].replace("-", "_")}_markdup1',
						template=mark_duplicates_samtools(
							alignment_file=align_paired1.outputs['alignment'],
							sample_name=sample['sample_name'],
							output_directory=f'{top_dir}/{GROUP_NAME}/no_merge'
						)
					)

					preFilterAlignmentsNoMerge.append(mark_duplicates1.outputs['markdup'])

					insertSizePre1 = gwf.target_from_template(
						name=f'{sample["sample_name"].replace("-", "_")}_pre_filter_insert_size1',
						template=insert_size(
							bamFile=mark_duplicates1.outputs['markdup'],
							outputName=f'{sample["sample_name"]}.preFilter',
							outputDirectory=f'{top_dir}/{GROUP_NAME}/no_merge'
						)
					)

					depthDistributionPre1 = gwf.target_from_template(
						name=f'{sample["sample_name"].replace("-", "_")}_pre_filter_depth_distribution1',
						template=depth_distribution(
							bamFile=mark_duplicates1.outputs['markdup'],
							lowerThreshold=50,
							outputDirectory=f'{top_dir}/{GROUP_NAME}/no_merge/depth',
							sampleName=f'{sample["sample_name"]}.preFilter.depthDist'
						)
					)

					depthAlongReferencePre1 = gwf.target_from_template(
						name=f'{sample["sample_name"].replace("-", "_")}_pre_filter_depth_along_reference1',
						template=depth_along_reference(
							bamFile=mark_duplicates1.outputs['markdup'],
							outputDirectory=f'{top_dir}/{GROUP_NAME}/no_merge/depth',
							sampleName=f'{sample["sample_name"]}.preFilter.depthAlong'
						)
					)

					insertDistributionPre1 = gwf.target_from_template(
						name=f'{sample["sample_name"].replace("-", "_")}_pre_filter_insert_distribution1',
						template=insert_size_distribution(
							insertSizeFile=insertSizePre1.outputs['insert'],
							outputDirectory=f'{top_dir}/{GROUP_NAME}/no_merge/insert_size',
							sampleName=f'{sample["sample_name"]}.preFilter.insertDist'
						)
					)

					unmapped1 = gwf.target_from_template(
						name=f'{sample["sample_name"].replace("-", "_")}_extract_unmapped_reads1',
						template=extract_unmapped_reads(
							alignment_file=mark_duplicates1.outputs['markdup'],
							sample_name=sample['sample_name'],
							output_directory=f'{top_out}/{GROUP_NAME}/no_merge' if OUTPUT_DIR else f'{top_dir}/{GROUP_NAME}/no_merge/unmapped'
						)
					)

					pre_filter_stats1 = gwf.target_from_template(
						name=f'{sample["sample_name"].replace("-", "_")}_pre_filter_stats1',
						template=samtools_stats(
							alignment_file=mark_duplicates1.outputs['markdup'],
							output_directory=f'{top_out}/{GROUP_NAME}/no_merge/{sample['sample_name']}/qc/pre_filtering_stats' if OUTPUT_DIR else f'{top_dir}/{GROUP_NAME}/no_merge/pre_filtering_stats/{sample["sample_name"]}'
						)
					)

					filter_alignment1 = gwf.target_from_template(
						name=f'{sample["sample_name"].replace("-", "_")}_filtering1',
						template=samtools_filter(
							alignment_file=mark_duplicates1.outputs['markdup'],
							sample_name=sample['sample_name'],
							output_directory=f'{top_out}/{GROUP_NAME}/no_merge' if OUTPUT_DIR else f'{top_dir}/{GROUP_NAME}/no_merge/filtered_alignment',
							flags_excluded=STF_EXCLUDE,
							flags_required=STF_REQUIRED,
							min_mq=STF_MIN_MQ
						)
					)

					insertSizePost1 = gwf.target_from_template(
						name=f'{sample["sample_name"].replace("-", "_")}_post_filter_insert_size1',
						template=insert_size(
							bamFile=filter_alignment1.outputs['filtered'],
							outputName=f'{sample["sample_name"]}.postFilter',
							outputDirectory=f'{top_dir}/{GROUP_NAME}/no_merge'
						)
					)

					depthDistributionPost1 = gwf.target_from_template(
						name=f'{sample["sample_name"].replace("-", "_")}_post_filter_depth_distribution1',
						template=depth_distribution(
							bamFile=filter_alignment1.outputs['filtered'],
							lowerThreshold=50,
							outputDirectory=f'{top_dir}/{GROUP_NAME}/no_merge/depth',
							sampleName=f'{sample["sample_name"]}.postFilter.depthDist'
						)
					)

					depthAlongReferencePost1 = gwf.target_from_template(
						name=f'{sample["sample_name"].replace("-", "_")}_post_filter_depth_along_reference1',
						template=depth_along_reference(
							bamFile=filter_alignment1.outputs['filtered'],
							outputDirectory=f'{top_dir}/{GROUP_NAME}/no_merge/depth',
							sampleName=f'{sample["sample_name"]}.postFilter.depthAlong'
						)
					)

					insertDistributionPost1 = gwf.target_from_template(
						name=f'{sample["sample_name"].replace("-", "_")}_post_filter_insert_distribution1',
						template=insert_size_distribution(
							insertSizeFile=insertSizePost1.outputs['insert'],
							outputDirectory=f'{top_dir}/{GROUP_NAME}/no_merge/insert_size',
							sampleName=f'{sample["sample_name"]}.postFilter.insertDist'
						)
					)

					postFilterAlignmentsNoMerge.append(filter_alignment1.outputs['filtered'])

					post_filter_stats1 = gwf.target_from_template(
						name=f'{sample["sample_name"].replace("-", "_")}_post_filter_stats1',
						template=samtools_stats(
							alignment_file=filter_alignment1.outputs['filtered'],
							output_directory=f'{top_out}/{GROUP_NAME}/no_merge/{sample['sample_name']}/qc/pre_filtering_stats' if OUTPUT_DIR else f'{top_dir}/{GROUP_NAME}/no_merge/post_filtering_stats/{sample["sample_name"]}'
						)
					)

					qualimap1 = gwf.target_from_template(
						name=f'{sample["sample_name"].replace("-", "_")}_qualimap1',
						template=qc_qualimap(
							alignment_file=filter_alignment1.outputs['filtered'],
							output_directory=f'{top_out}/{GROUP_NAME}/no_merge/{sample['sample_name']}/qc/bamqc' if OUTPUT_DIR else f'{top_dir}/bamqc/individual/{GROUP_NAME}/no_merge/{sample["sample_name"]}'
						)
					)

				elif i == 2:
					align_paired2 = gwf.target_from_template(
						name=f'{sample["sample_name"].replace("-", "_")}_paired_alignment2',
						template=alignment_pairedend(
							read1_file=adapterremoval_no_collapse.outputs['pair1'],
							read2_file=adapterremoval_no_collapse.outputs['pair2'],
							reference_genome_file=index.outputs['symlink'],
							sample_name=sample['sample_name'],
							output_directory=f'{top_dir}/{GROUP_NAME}/no_collapse'
						)
					)

					mark_duplicates2 = gwf.target_from_template(
						name=f'{sample["sample_name"].replace("-", "_")}_markdup2',
						template=mark_duplicates_samtools(
							alignment_file=align_paired2.outputs['alignment'],
							sample_name=sample['sample_name'],
							output_directory=f'{top_dir}/{GROUP_NAME}/no_collapse'
						)
					)

					preFilterAlignmentsNoCollapse.append(mark_duplicates2.outputs['markdup'])

					insertSizePre2 = gwf.target_from_template(
						name=f'{sample["sample_name"].replace("-", "_")}_pre_filter_insert_size2',
						template=insert_size(
							bamFile=mark_duplicates2.outputs['markdup'],
							outputName=f'{sample["sample_name"]}.preFilter',
							outputDirectory=f'{top_dir}/{GROUP_NAME}/no_collapse'
						)
					)

					depthDistributionPre2 = gwf.target_from_template(
						name=f'{sample["sample_name"].replace("-", "_")}_pre_filter_depth_distribution2',
						template=depth_distribution(
							bamFile=mark_duplicates2.outputs['markdup'],
							lowerThreshold=50,
							outputDirectory=f'{top_dir}/{GROUP_NAME}/no_collapse/depth',
							sampleName=f'{sample["sample_name"]}.preFilter.depthDist'
						)
					)

					depthDistributionPreOverlapOnce2 = gwf.target_from_template(
						name=f'{sample["sample_name"].replace("-", "_")}_pre_filter_depth_distribution_overlap_once2',
						template=depth_distribution(
							bamFile=mark_duplicates2.outputs['markdup'],
							lowerThreshold=50,
							outputDirectory=f'{top_dir}/{GROUP_NAME}/no_collapse/depth',
							sampleName=f'{sample["sample_name"]}.preFilter.depthDist.overlapOnce',
							additionalOptions='-s'
						)
					)

					depthAlongReferencePre2 = gwf.target_from_template(
						name=f'{sample["sample_name"].replace("-", "_")}_pre_filter_depth_along_reference2',
						template=depth_along_reference(
							bamFile=mark_duplicates2.outputs['markdup'],
							outputDirectory=f'{top_dir}/{GROUP_NAME}/no_collapse/depth',
							sampleName=f'{sample["sample_name"]}.preFilter.depthAlong'
						)
					)

					insertDistributionPre2 = gwf.target_from_template(
						name=f'{sample["sample_name"].replace("-", "_")}_pre_filter_insert_distribution2',
						template=insert_size_distribution(
							insertSizeFile=insertSizePre2.outputs['insert'],
							outputDirectory=f'{top_dir}/{GROUP_NAME}/no_collapse/insert_size',
							sampleName=f'{sample["sample_name"]}.preFilter.insertDist'
						)
					)

					unmapped2 = gwf.target_from_template(
						name=f'{sample["sample_name"].replace("-", "_")}_extract_unmapped_reads2',
						template=extract_unmapped_reads(
							alignment_file=mark_duplicates2.outputs['markdup'],
							sample_name=sample['sample_name'],
							output_directory=f'{top_out}/{GROUP_NAME}/no_collapse' if OUTPUT_DIR else f'{top_dir}/{GROUP_NAME}/no_collapse/unmapped'
						)
					)

					pre_filter_stats2 = gwf.target_from_template(
						name=f'{sample["sample_name"].replace("-", "_")}_pre_filter_stats2',
						template=samtools_stats(
							alignment_file=mark_duplicates2.outputs['markdup'],
							output_directory=f'{top_out}/{GROUP_NAME}/no_collapse/{sample['sample_name']}/qc/pre_filtering_stats' if OUTPUT_DIR else f'{top_dir}/{GROUP_NAME}/no_collapse/pre_filtering_stats/{sample["sample_name"]}'
						)
					)

					filter_alignment2 = gwf.target_from_template(
						name=f'{sample["sample_name"].replace("-", "_")}_filtering2',
						template=samtools_filter(
							alignment_file=mark_duplicates2.outputs['markdup'],
							sample_name=sample['sample_name'],
							output_directory=f'{top_out}/{GROUP_NAME}/no_collapse' if OUTPUT_DIR else f'{top_dir}/{GROUP_NAME}/no_collapse/filtered_alignment',
							flags_excluded=STF_EXCLUDE,
							flags_required=STF_REQUIRED,
							min_mq=STF_MIN_MQ
						)
					)

					postFilterAlignmentsNoCollapse.append(filter_alignment2.outputs['filtered'])

					insertSizePost2 = gwf.target_from_template(
						name=f'{sample["sample_name"].replace("-", "_")}_post_filter_insert_size2',
						template=insert_size(
							bamFile=filter_alignment2.outputs['filtered'],
							outputName=f'{sample["sample_name"]}.postFilter',
							outputDirectory=f'{top_dir}/{GROUP_NAME}/no_collapse'
						)
					)

					depthDistributionPost2 = gwf.target_from_template(
						name=f'{sample["sample_name"].replace("-", "_")}_post_filter_depth_distribution2',
						template=depth_distribution(
							bamFile=filter_alignment2.outputs['filtered'],
							lowerThreshold=50,
							outputDirectory=f'{top_dir}/{GROUP_NAME}/no_collapse/depth',
							sampleName=f'{sample["sample_name"]}.postFilter.depthDist'
						)
					)

					depthDistributionPostOverlapOnce2 = gwf.target_from_template(
						name=f'{sample["sample_name"].replace("-", "_")}_post_filter_depth_distribution_overlap_once2',
						template=depth_distribution(
							bamFile=filter_alignment2.outputs['filtered'],
							lowerThreshold=50,
							outputDirectory=f'{top_dir}/{GROUP_NAME}/no_collapse/depth',
							sampleName=f'{sample["sample_name"]}.postFilter.depthDist.overlapOnce',
							additionalOptions='-s'
						)
					)

					depthAlongReferencePost2 = gwf.target_from_template(
						name=f'{sample["sample_name"].replace("-", "_")}_post_filter_depth_along_reference2',
						template=depth_along_reference(
							bamFile=filter_alignment2.outputs['filtered'],
							outputDirectory=f'{top_dir}/{GROUP_NAME}/no_collapse/depth',
							sampleName=f'{sample["sample_name"]}.postFilter.depthAlong'
						)
					)

					insertDistributionPost2 = gwf.target_from_template(
						name=f'{sample["sample_name"].replace("-", "_")}_post_filter_insert_distribution2',
						template=insert_size_distribution(
							insertSizeFile=insertSizePost2.outputs['insert'],
							outputDirectory=f'{top_dir}/{GROUP_NAME}/no_collapse/insert_size',
							sampleName=f'{sample["sample_name"]}.postFilter.insertDist'
						)
					)

					post_filter_stats2 = gwf.target_from_template(
						name=f'{sample["sample_name"].replace("-", "_")}_post_filter_stats2',
						template=samtools_stats(
							alignment_file=filter_alignment2.outputs['filtered'],
							output_directory=f'{top_out}/{GROUP_NAME}/no_collapse/{sample['sample_name']}/qc/pre_filtering_stats' if OUTPUT_DIR else f'{top_dir}/{GROUP_NAME}/no_collapse/post_filtering_stats/{sample["sample_name"]}'
						)
					)

					qualimap2 = gwf.target_from_template(
						name=f'{sample["sample_name"].replace("-", "_")}_qualimap2',
						template=qc_qualimap(
							alignment_file=filter_alignment2.outputs['filtered'],
							output_directory=f'{top_out}/{GROUP_NAME}/no_collapse/{sample['sample_name']}/qc/bamqc' if OUTPUT_DIR else f'{top_dir}/bamqc/individual/{GROUP_NAME}/no_collapse/{sample["sample_name"]}'
						)
					)

			resultReport = gwf.target_from_template(
				name=f'result_report_{sample["sample_name"].replace("-", "_")}',
				template=result_report(
					sampleName=sample['sample_name'],
					preFlagstatCollapse=pre_filter_stats0.outputs['stats'][1],
					preStatsCollapse=pre_filter_stats0.outputs['stats'][3],
					preDepthDistributionCollapse=depthDistributionPre0.outputs['png'],
					preDepthAlongReferenceCollapse=depthAlongReferencePre0.outputs['png'],
					postFlagstatCollapse=post_filter_stats0.outputs['stats'][1],
					postStatsCollapse=post_filter_stats0.outputs['stats'][3],
					postDepthDistributionCollapse=depthDistributionPost0.outputs['png'],
					postDepthAlongReferenceCollapse=depthAlongReferencePost0.outputs['png'],
					preFlagstatDiscard=pre_filter_stats1.outputs['stats'][1],
					preStatsDiscard=pre_filter_stats1.outputs['stats'][3],
					preDepthDistributionDiscard=depthDistributionPre1.outputs['png'],
					preDepthAlongReferenceDiscard=depthAlongReferencePre1.outputs['png'],
					postFlagstatDiscard=post_filter_stats1.outputs['stats'][1],
					postStatsDiscard=post_filter_stats1.outputs['stats'][3],
					postDepthDistributionDiscard=depthDistributionPost1.outputs['png'],
					postDepthAlongReferenceDiscard=depthAlongReferencePost1.outputs['png'],
					preFlagstatIgnore=pre_filter_stats2.outputs['stats'][1],
					preStatsIgnore=pre_filter_stats2.outputs['stats'][3],
					preDepthDistributionIgnore=depthDistributionPreOverlapOnce2.outputs['png'],
					preDepthAlongReferenceIgnore=depthAlongReferencePre2.outputs['png'],
					postFlagstatIgnore=post_filter_stats2.outputs['stats'][1],
					postStatsIgnore=post_filter_stats2.outputs['stats'][3],
					postDepthDistributionIgnore=depthDistributionPostOverlapOnce2.outputs['png'],
					postDepthAlongReferenceIgnore=depthAlongReferencePost2.outputs['png'],
					outputDirectory=f'{top_dir}/{GROUP_NAME}/{sample["sample_name"]}/report'
				)
			)
		
		# perSiteCoverageMergePre = gwf.target_from_template(
		# 	name=f'per_site_coverage_merge_pre',
		# 	template=per_site_coverage(
		# 		bamFiles=preFilterAlignmentsMerge,
		# 		outputName=f'{species_abbreviation(SPECIES_NAME)}.preFilter',
		# 		outputDirectory=f'{top_dir}/{GROUP_NAME}/merge/depth'
		# 	)
		# )

		# perSiteCoverageMergePost = gwf.target_from_template(
		# 	name=f'per_site_coverage_merge_post',
		# 	template=per_site_coverage(
		# 		bamFiles=postFilterAlignmentsMerge,
		# 		outputName=f'{species_abbreviation(SPECIES_NAME)}.postFilter',
		# 		outputDirectory=f'{top_dir}/{GROUP_NAME}/merge/depth'
		# 	)
		# )

		# perSiteCoverageNoMergePre = gwf.target_from_template(
		# 	name=f'per_site_coverage_no_merge_pre',
		# 	template=per_site_coverage(
		# 		bamFiles=preFilterAlignmentsNoMerge,
		# 		outputName=f'{species_abbreviation(SPECIES_NAME)}.preFilter',
		# 		outputDirectory=f'{top_dir}/{GROUP_NAME}/no_merge/depth'
		# 	)
		# )

		# perSiteCoverageNoMergePost = gwf.target_from_template(
		# 	name=f'per_site_coverage_no_merge_post',
		# 	template=per_site_coverage(
		# 		bamFiles=postFilterAlignmentsNoMerge,
		# 		outputName=f'{species_abbreviation(SPECIES_NAME)}.postFilter',
		# 		outputDirectory=f'{top_dir}/{GROUP_NAME}/no_merge/depth'
		# 	)
		# )

		# perSiteCoverageNoCollapsePre = gwf.target_from_template(
		# 	name=f'per_site_coverage_no_collapse_pre',
		# 	template=per_site_coverage(
		# 		bamFiles=preFilterAlignmentsNoCollapse,
		# 		outputName=f'{species_abbreviation(SPECIES_NAME)}.preFilter',
		# 		outputDirectory=f'{top_dir}/{GROUP_NAME}/no_collapse/depth'
		# 	)
		# )

		# perSiteCoverageNoCollapsePost = gwf.target_from_template(
		# 	name=f'per_site_coverage_no_collapse_post',
		# 	template=per_site_coverage(
		# 		bamFiles=postFilterAlignmentsNoCollapse,
		# 		outputName=f'{species_abbreviation(SPECIES_NAME)}.postFilter',
		# 		outputDirectory=f'{top_dir}/{GROUP_NAME}/no_collapse/depth'
		# 	)
		# )

	return gwf