#!/bin/env python3
from gwf import Workflow
from gwf.workflow import collect
import os, yaml, glob, sys
from workflow_templates import *

def freebayes_vcf_workflow(config_file: str = glob.glob('*config.y*ml')[0]):
	"""
	Workflow: Create :format:`VCF` using :script:`freebayes`.
	
	:param str config_file:
		Configuration file containing pre-defined set of variables
	"""
	# --------------------------------------------------
	#                  Configuration
	# --------------------------------------------------
	
	CONFIG = yaml.safe_load(open(config_file))
	ACCOUNT: str = CONFIG['account']
	SPECIES_NAME: str = CONFIG['species_name']
	SAMPLE_LIST: list = CONFIG['sample_list']
	REFERENCE_GENOME: str = CONFIG['reference_genome_path']
	OUTPUT_DIR: str = CONFIG['output_directory_path']
	WORK_DIR: str = CONFIG['working_directory_path']
	PARTITION_SIZE: int = CONFIG['partition_size']
	PLOIDY: int = CONFIG['sample_ploidy']
	BESTN: int = CONFIG['best_n_alleles']
	ALT_FRACTION: float | int = CONFIG['min_alternate_fraction']
	ALT_COUNT: int = CONFIG['min_alternate_count']

	# --------------------------------------------------
	#                  Workflow
	# --------------------------------------------------
	
	gwf = Workflow(
		defaults={'account': ACCOUNT}
	)
	
	if os.path.exists(f'reference_partitions.{PARTITION_SIZE}bp.txt') and os.path.exists('reference_sequences.txt'):
		# If files exists reads data directly from files
		# Loads reference genome partitioning
		with open(f'reference_partitions.{PARTITION_SIZE}bp.txt', 'r') as infile:
			partitions = [{'num': entry.split(sep='\t')[0].strip(), 'region': entry.split(sep='\t')[1].strip(), 'start': entry.split(sep='\t')[2].strip(), 'end': entry.split(sep='\t')[3].strip()} for entry in infile]
			npadding = len(str(sum(1 for line in partitions)))
		# Loads list of contigs in reference genome
		with open('reference_sequences.txt', 'r') as infile:
			contigs = [{'contig': entry.split(sep='\t')[0].strip()} for entry in infile]
	else:
		# If files don't exist, generate data and write files
		sequences = parse_fasta(REFERENCE_GENOME)
		with open(f'reference_sequences.txt', 'w') as outfile:
			outfile.write('\n'.join('\t'.join(str(i) for i in entry.values()) for entry in sequences))
		# Partitions reference genome
		npadding = padding_calculator(parse_fasta=sequences, size=PARTITION_SIZE)
		partitions = partition_chrom(parse_fasta=sequences, size=PARTITION_SIZE, npad=npadding)
		with open(f'reference_partitions.{PARTITION_SIZE}bp.txt', 'w') as outfile:
			outfile.write('\n'.join('\t'.join(str(i) for i in entry.values()) for entry in partitions))
		# Creates list of contigs in reference genome
		contigs = [{'contig': contig['sequence_name']} for contig in sequences]

	top_dir = f'{WORK_DIR}/{SPECIES_NAME.replace(" ", "_")}/freebayes'
	output_dir = f'{OUTPUT_DIR}/{SPECIES_NAME.replace(" ", "_")}/freebayes'

	freebayes_parts = gwf.map(
		template_func=freebayes_chrom,
		inputs=partitions,
		name=name_freebayes_chrom,
		extra={'reference_genome_file': REFERENCE_GENOME,
			   'bam_file_list': SAMPLE_LIST,
			   'output_directory': top_dir,
			   'species_name': SPECIES_NAME,
			   'ploidy': PLOIDY,
			   'best_n_alleles': BESTN,
			   'min_alternate_fraction': ALT_FRACTION,
			   'min_alternate_count': ALT_COUNT}
		)
	
	concat_freebayes = gwf.target_from_template(
			name='concatenate_freebayes',
			template=concat_vcf(
				files=collect(freebayes_parts.outputs, ['vcf'])['vcfs'],
				output_name=f'{species_abbreviation(SPECIES_NAME)}.freebayes_n{BESTN}_p{PLOIDY}_minaltfrc{ALT_FRACTION}_minaltcnt{ALT_COUNT}',
				output_directory=output_dir,
				compress=True
			)
		)

	return gwf

def freebayes_population_set_workflow(config_file: str = glob.glob('*config.y*ml')[0]):
	"""
	Workflow: Create :format:`VCF` file for each sample in configuration.
	
	:param str config_file:
		Configuration file containing pre-defined set of variables
	"""
	# --------------------------------------------------
	#                  Configuration
	# --------------------------------------------------
	
	CONFIG = yaml.safe_load(open(config_file))
	ACCOUNT: str = CONFIG['account']
	TAXONOMY: str = CONFIG['taxonomic_group'].lower() if CONFIG['taxonomic_group'] else None
	SPECIES_NAME: str = CONFIG['species_name']
	REFERENCE_GENOME: str = CONFIG['reference_genome_path']
	WORK_DIR: str = CONFIG['working_directory_path'][:len(CONFIG['working_directory_path']) - 1] if CONFIG['working_directory_path'].endswith('/') else CONFIG['working_directory_path']
	OUTPUT_DIR: str = (CONFIG['output_directory_path'][:len(CONFIG['output_directory_path']) - 1] if CONFIG['output_directory_path'].endswith('/') else CONFIG['output_directory_path']) if CONFIG['output_directory_path'] else None
	PARTITION_SIZE: int | None = CONFIG['partition_size'] if CONFIG['partition_size'] else 500000
	FREEBAYES_SETTINGS: dict = CONFIG['freebayes_settings']
	FREEBAYES_PLOIDY: int | None = FREEBAYES_SETTINGS['sample_ploidy'] if FREEBAYES_SETTINGS['sample_ploidy'] else 100
	FREEBAYES_BESTN: int | None = FREEBAYES_SETTINGS['best_n_alleles'] if FREEBAYES_SETTINGS['best_n_alleles'] else 3
	FREEBAYES_MINALTFRC: float | int | None = FREEBAYES_SETTINGS['min_alternate_fraction'] if FREEBAYES_SETTINGS['min_alternate_fraction'] else 0
	FREEBAYES_MINALTCNT: int | None = FREEBAYES_SETTINGS['min_alternate_count'] if FREEBAYES_SETTINGS['min_alternate_count'] else 2
	FILTERING: dict = CONFIG['filtering']
	FILTERING_MINDP: int | None = FILTERING['minimum_depth'] if FILTERING['minimum_depth'] else 300
	FILTERING_MAXDP: int | None = FILTERING['maximum_depth'] if FILTERING['maximum_depth'] else 600
	MODE: int = CONFIG['mode']
	SAMPLE_LIST: list = CONFIG['sample_list']

	# --------------------------------------------------
	#                  Workflow
	# --------------------------------------------------
	
	gwf = Workflow(
		defaults={'account': ACCOUNT}
	)
	
	if os.path.exists(f'reference_partitions.{PARTITION_SIZE}bp.txt') and os.path.exists('reference_sequences.txt'):
		# If files exists reads data directly from files
		# Loads reference genome partitioning
		with open(f'reference_partitions.{PARTITION_SIZE}bp.txt', 'r') as infile:
			partitions = [{'num': entry.split(sep='\t')[0].strip(), 'region': entry.split(sep='\t')[1].strip(), 'start': entry.split(sep='\t')[2].strip(), 'end': entry.split(sep='\t')[3].strip()} for entry in infile]
			npadding = len(str(sum(1 for line in partitions)))
		# Loads list of contigs in reference genome
		with open('reference_sequences.txt', 'r') as infile:
			contigs = [{'contig': entry.split(sep='\t')[0].strip()} for entry in infile]
	else:
		# If files don't exist, generate data and write files
		sequences = parse_fasta(REFERENCE_GENOME)
		with open(f'reference_sequences.txt', 'w') as outfile:
			outfile.write('\n'.join('\t'.join(str(i) for i in entry.values()) for entry in sequences))
		# Partitions reference genome
		npadding = padding_calculator(parse_fasta=sequences, size=PARTITION_SIZE)
		partitions = partition_chrom(parse_fasta=sequences, size=PARTITION_SIZE, npad=npadding)
		with open(f'reference_partitions.{PARTITION_SIZE}bp.txt', 'w') as outfile:
			outfile.write('\n'.join('\t'.join(str(i) for i in entry.values()) for entry in partitions))
		# Creates list of contigs in reference genome
		contigs = [{'contig': contig['sequence_name']} for contig in sequences]

	top_dir = f'{WORK_DIR}/{TAXONOMY.replace(" ", "_")}/{SPECIES_NAME.replace(" ", "_")}/vcf' if TAXONOMY else f'{WORK_DIR}/{SPECIES_NAME.replace(" ", "_")}/vcf'
	top_out = f'{OUTPUT_DIR}/vcf/{TAXONOMY.replace(" ", "_")}/{SPECIES_NAME.replace(" ", "_")}' if TAXONOMY else f'{OUTPUT_DIR}/vcf/{SPECIES_NAME.replace(" ", "_")}'
	full_bam_list = []

	for GROUP in SAMPLE_LIST:
		if not GROUP['bam_file_list']:
			continue
		GROUP_NAME: str = GROUP['group_name'].lower()

		for SAMPLE in GROUP['bam_file_list']:
			SAMPLE_NAME: str = os.path.basename(os.path.dirname(SAMPLE))
			
			full_bam_list.append(SAMPLE)

			# One VCF file per sample file.
			if MODE == 1 or MODE == 4:
				freebayes_parts_single = gwf.map(
					name=name_freebayes_partition_single,
					template_func=freebayes_partition_single,
					inputs=partitions,
					extra={'reference_genome_file': REFERENCE_GENOME,
						   'bam_file': SAMPLE,
						   'output_directory': top_dir,
						   'group_name': GROUP_NAME,
						   'sample_name': SAMPLE_NAME,
						   'ploidy': FREEBAYES_PLOIDY,
						   'best_n_alleles': FREEBAYES_BESTN,
						   'min_alternate_fraction': FREEBAYES_MINALTFRC,
						   'min_alternate_count': FREEBAYES_MINALTCNT}
				)

				concat_freebayes_single = gwf.target_from_template(
						name=f'concatenate_freebayes_vcf_{GROUP_NAME}_{SAMPLE_NAME.replace("-", "_")}',
						template=concat_vcf(
							files=collect(freebayes_parts_single.outputs, ['vcf'])['vcfs'],
							output_name=f'{SAMPLE_NAME}.freebayes_n{FREEBAYES_BESTN}_p{FREEBAYES_PLOIDY}_minaltfrc{FREEBAYES_MINALTFRC}_minaltcnt{FREEBAYES_MINALTCNT}',
							output_directory=f'{top_out}/{GROUP_NAME}/{SAMPLE_NAME}/vcf' if OUTPUT_DIR else f'{top_dir}/raw_vcf/{GROUP_NAME}/{SAMPLE_NAME}',
							compress=True
						)
				)

				# filtering = gwf.target_from_template(
				# 	name=f'filter_vcf_{GROUP_NAME}_{SAMPLE_NAME.replace("-", "_")}',
				# 	template=filter_vcf(
				# 		vcf_file=concat_freebayes_single.outputs['concat_file'],
				# 		output_directory=f'{OUTPUT_DIR}' if OUTPUT_DIR else f'{top_dir}/filtered_vcf',
				# 		sample_group=GROUP_NAME,
				# 		sample_name=SAMPLE_NAME,
				# 		min_depth=FILTERING_MINDP,
				# 		max_depth=FILTERING_MAXDP
				# 	)
				# )
		
		# One VCF file per sample group.
		if MODE == 2 or MODE == 4:
			freebayes_parts_group = gwf.map(
				name=name_freebayes_partition_group,
				template_func=freebayes_partition_group,
				inputs=partitions,
				extra={'reference_genome_file': REFERENCE_GENOME,
		  		   	   'bam_files': GROUP['bam_file_list'],
				   	   'output_directory': top_dir,
					   'species_name': SPECIES_NAME,
					   'group_name': GROUP_NAME,
				   	   'ploidy': FREEBAYES_PLOIDY,
				   	   'best_n_alleles': FREEBAYES_BESTN,
				   	   'min_alternate_fraction': FREEBAYES_MINALTFRC,
				   	   'min_alternate_count': FREEBAYES_MINALTCNT}
			)
		
			concat_freebayes_group = gwf.target_from_template(
				name=f'cocatenate_freebayes_vcf_group',
				template=concat_vcf(
					files=collect(freebayes_parts_group.outputs, ['vcf'])['vcfs'],
					output_name=f'{species_abbreviation(SPECIES_NAME)}_{GROUP_NAME}.freebayes_n{FREEBAYES_BESTN}_p{FREEBAYES_PLOIDY}_minaltfrc{FREEBAYES_MINALTFRC}_minaltcnt{FREEBAYES_MINALTCNT}',
					output_directory=f'{top_out}/{GROUP_NAME}/vcf' if OUTPUT_DIR else f'{top_dir}/raw_vcf/{GROUP_NAME}'
				)
			)
	
	# One VCF file containing all sample files.
	if MODE == 3 or MODE == 4:
		freebayes_parts_all = gwf.map(
			name=name_freebayes_partition_all,
			template_func=freebayes_partition_all,
			inputs=partitions,
			extra={'reference_genome_file': REFERENCE_GENOME,
		  		   'bam_files': full_bam_list,
				   'output_directory': top_dir,
				   'species_name': SPECIES_NAME,
				   'ploidy': FREEBAYES_PLOIDY,
				   'best_n_alleles': FREEBAYES_BESTN,
				   'min_alternate_fraction': FREEBAYES_MINALTFRC,
				   'min_alternate_count': FREEBAYES_MINALTCNT}
		)

		concat_freebayes_all = gwf.target_from_template(
			name=f'concatenate_freebayes_vcf_all',
			template=concat_vcf(
				files=collect(freebayes_parts_all.outputs, ['vcf'])['vcfs'],
				output_name=f'{species_abbreviation(SPECIES_NAME)}.freebayes_n{FREEBAYES_BESTN}_p{FREEBAYES_PLOIDY}_minaltfrc{FREEBAYES_MINALTFRC}_minaltcnt{FREEBAYES_MINALTCNT}',
				output_directory=f'{top_out}/vcf' if OUTPUT_DIR else f'{top_dir}/raw_vcf',
				compress=True
			)
		)
	
	return gwf
