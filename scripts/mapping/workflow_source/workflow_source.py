#!/bin/env python3
from gwf import Workflow
from gwf.workflow import collect
import os, yaml, glob, sys
from workflow_templates import *

def mapping_resequencing_data_workflow(config_file: str = glob.glob('*config.y*ml')[0]):
	"""
	Workflow: Align resequencing data to reference genome and do basic filtering.
	
	:param str config_file:
		Configuration file containing pre-defined set of variables
	"""
	# --------------------------------------------------
	#                  Configuration
	# --------------------------------------------------
	
	config = yaml.safe_load(open(config_file))
	ACCOUNT: str = config['account']
	SPECIES_NAME: str = config['species_name']
	REFERENCE_GENOME: str = config['reference_genome_path']
	SAMPLE_GROUP1: dict = config['sample_group1']
	# SAMPLE_DIR1: list = config['sample_folder_list1']
	# SAMPLE_DIR2_NAME: str = config['sample_folder_list2_name']
	# SAMPLE_DIR2: list = config['sample_folder_list2']
	# SAMPLE_DIR3_NAME: str = config['sample_folder_list3_name']
	# SAMPLE_DIR3: list = config['sample_folder_list3']
	WORK_DIR: str = config['working_directory_path']
	OUTPUT_DIR: str = config['output_directory_path']
	
	# --------------------------------------------------
	#                  Workflow
	# --------------------------------------------------
	
	gwf = Workflow(
		defaults={'account': ACCOUNT}
	)
	
	top_dir = f'{WORK_DIR}/{SPECIES_NAME.replace(" ", "_")}/mapping'
	
	
	
	print(SAMPLE_GROUP1)

	# all_samples = [i for j in [SAMPLE_DIR1, SAMPLE_DIR2, SAMPLE_DIR3] if j for i in j]
	# sample_list = [get_sample_data(path) for path in all_samples]
	# qualimap_list = []

	# index = gwf.target_from_template(
	# 	name=f'{SPECIES_NAME.replace(" ", "_")}_index_reference',
	# 	template=index_reference_genome(
	# 		reference_genome_file=REFERENCE_GENOME,
	# 		output_directory=top_dir)
	# )

	# for sample in sample_list:	
	# 	if OUTPUT_DIR == WORK_DIR:
	# 		OUTPUT_DIR = f'{top_dir}/alignment/{sample["sample_name"]}'
	# 	else:
	# 		OUTPUT_DIR = f'{OUTPUT_DIR}/{sample["sample_name"]}'

	# 	adapterremoval = gwf.target_from_template(
	# 		name=f'{sample["sample_name"].replace("-", "_")}_adapterremoval',
	# 		template=adapterremoval_pairedend(
	# 			sample_name=sample['sample_name'],
	# 			read1_files=sample['read1_files'],
	# 			read2_files=sample['read2_files'],
	# 			output_directory=top_dir
	# 		)
	# 	)

	# 	align_paired = gwf.target_from_template(
	# 		name=f'{sample["sample_name"].replace("-", "_")}_paired_alignment',
	# 		template=alignment_pairedend(
	# 			read1_file=adapterremoval.outputs['pair1'],
	# 			read2_file=adapterremoval.outputs['pair2'],
	# 			reference_genome_file=index.outputs['symlink'],
	# 			sample_name=sample['sample_name'],
	# 			output_directory=top_dir
	# 		)
	# 	)

	# 	align_collapsed = gwf.target_from_template(
	# 		name=f'{sample["sample_name"].replace("-", "_")}_collapsed_alignment',
	# 		template=alignment_collapsed(
	# 			collapsed_read_files=adapterremoval.outputs['collapsed'],
	# 			reference_genome_file=index.outputs['symlink'],
	# 			sample_name=sample['sample_name'],
	# 			output_directory=top_dir
	# 		)
	# 	)

	# 	merge = gwf.target_from_template(
	# 		name=f'{sample["sample_name"].replace("-", "_")}_merge_alignments',
	# 		template=merge_alignments(
	# 			alignment_files=[align_paired.outputs['alignment'],
	# 				 			 align_collapsed.outputs['alignment']],
	# 			sample_name=sample['sample_name'],
	# 			output_directory=top_dir
	# 		)
	# 	)

	# 	mark_duplicates = gwf.target_from_template(
	# 		name=f'{sample["sample_name"].replace("-", "_")}_markdup',
	# 		template=mark_duplicates_samtools(
	# 			alignment_file=merge.outputs['merged'],
	# 			sample_name=sample['sample_name'],
	# 			output_directory=top_dir
	# 		)
	# 	)

	# 	unmapped = gwf.target_from_template(
	# 		name=f'{sample["sample_name"].replace("-", "_")}_extract_unmapped_reads',
	# 		template=extract_unmapped_reads(
	# 			alignment_file=mark_duplicates.outputs['markdup'],
	# 			sample_name=sample['sample_name'],
	# 			output_directory=OUTPUT_DIR
	# 		)
	# 	)

	# 	pre_filter_stats = gwf.target_from_template(
	# 		name=f'{sample["sample_name"].replace("-", "_")}_pre_filter_stats',
	# 		template=samtools_stats(
	# 			alignment_file=mark_duplicates.outputs['markdup'],
	# 			output_directory=f'{top_dir}/alignment/{sample["sample_name"]}/pre_filtering_stats'
	# 		)
	# 	)

	# 	filter_alignment = gwf.target_from_template(
	# 		name=f'{sample["sample_name"].replace("-", "_")}_filtering',
	# 		template=samtools_filter(
	# 			alignment_file=mark_duplicates.outputs['markdup'],
	# 			sample_name=sample['sample_name'],
	# 			output_directory=OUTPUT_DIR
	# 		)
	# 	)

	# 	post_filter_stats = gwf.target_from_template(
	# 		name=f'{sample["sample_name"].replace("-", "_")}_post_filter_stats',
	# 		template=samtools_stats(
	# 			alignment_file=filter_alignment.outputs['filtered'],
	# 			output_directory=f'{top_dir}/alignment/{sample["sample_name"]}/post_filtering_stats'
	# 		)
	# 	)

	# 	qualimap = gwf.target_from_template(
	# 		name=f'{sample["sample_name"].replace("-", "_")}_qualimap',
	# 		template=qc_qualimap(
	# 			alignment_file=filter_alignment.outputs['filtered'],
	# 			sample_name=sample['sample_name'],
	# 			output_directory=OUTPUT_DIR
	# 		)
	# 	)

	# 	if 
	# 	qualimap_list.append({'sample_name': sample['sample_name'], 'path': qualimap.outputs['raw'], 'group': sample['group']})

	return gwf