#!/bin/env python3
from gwf import Workflow
from gwf.workflow import collect
import os, yaml, glob, sys, gzip
from workflow_templates import *

def genetic_load_workflow(config_file: str = glob.glob('*config.y*ml')[0]):
	"""
	Workflow: description
	
	:param str config_file:
		Configuration file containing pre-defined set of variables
	"""
	# --------------------------------------------------
	#                  Configuration
	# --------------------------------------------------
	
	CONFIG = yaml.safe_load(open(config_file))
	ACCOUNT: str = CONFIG['account']
	SPECIES_NAME: str = CONFIG['species_name']
	REFERENCE: str = CONFIG['reference_genome_file']
	GTF: str = CONFIG['gtf_annotation_file']
	WORK_DIR: str = CONFIG['working_directory_path']
	OUTPUT_DIR: str = CONFIG['output_directory_path']
	SNPGENIE_SETTINGS: dict = CONFIG['snpgenie_settings']
	SNPGENIE_MINFREQ: float | int = SNPGENIE_SETTINGS['minfreq'] if SNPGENIE_SETTINGS['minfreq'] else 0
	SNPGENIE_SLIDINGWINDOW: int = SNPGENIE_SETTINGS['slidingwindow'] if SNPGENIE_SETTINGS['slidingwindow'] else 9
	SAMPLES: list = CONFIG['sample_list']
	
	SNPEFF_DIR = f'{os.path.dirname(os.path.realpath(__file__))}/software/snpeff'

	# --------------------------------------------------
	#                  Workflow
	# --------------------------------------------------
	
	gwf = Workflow(
		defaults={'account': ACCOUNT}
	)

	top_dir = f'{WORK_DIR}/{SPECIES_NAME.replace(" ", "_")}/genetic_load'

	if os.path.exists(f'reference_sequences_{os.path.splitext(os.path.basename(REFERENCE))[0]}.txt'):
		with open(f'reference_sequences_{os.path.splitext(os.path.basename(REFERENCE))[0]}.txt', 'r') as infile:
			sequences = [{'sequence_name': entry.split(sep='\t')[0].strip(), 'sequence_length': entry.split(sep='\t')[1].strip()} for entry in infile]
	else:
		sequences = parse_fasta(REFERENCE)
		with open(f'reference_sequences_{os.path.splitext(os.path.basename(REFERENCE))[0]}.txt', 'w') as outfile:
			outfile.write('\n'.join('\t'.join(str(i) for i in entry.values()) for entry in sequences))

	database_entry = gwf.target_from_template(
		name=f'{species_abbreviation(SPECIES_NAME)}_snpeff_database_entry',
		template=snpeff_database_build(
			gtf_annotation_file=GTF,
			reference_genome_file=REFERENCE,
			species_name=SPECIES_NAME,
			snpeff_directory=SNPEFF_DIR
		)
	)
	
	snpeff_results_list = []
	snpgenie_results_list = []

	for sample in SAMPLES:
		SAMPLE_NAME = sample['sample_name']
		SAMPLE_GROUP = sample['sample_group']
		VCF = sample['vcf_file']
		BAM = sample['bam_file']
		variant_annotation = gwf.target_from_template(
			name=f'snpeff_annotation_{SAMPLE_GROUP}_{SAMPLE_NAME.replace("-", "_")}',
			template=snpeff_annotation(
				vcf_file=VCF,
				snpeff_predictor_file=database_entry.outputs['predictor'],
				snpeff_config_file=f'{SNPEFF_DIR}/snpEff.config',
				output_directory=top_dir,
				sample_group=SAMPLE_GROUP,
				sample_name=SAMPLE_NAME
			)
		)

		annotation_frequencies = gwf.target_from_template(
			name=f'snpeff_frequencies_{SAMPLE_GROUP}_{SAMPLE_NAME.replace("-", "_")}',
			template=snpeff_freqs(
				ann=variant_annotation.outputs['ann']
			)
		)

		site_count = gwf.target_from_template(
			name=f'cds_site_count_{SAMPLE_GROUP}_{SAMPLE_NAME.replace("-", "_")}',
			template=cds_site_count(
				bam_file=BAM,
				gtf_annotation_file=GTF,
				output_directory=top_dir,
				sample_group=SAMPLE_GROUP,
				sample_name=SAMPLE_NAME,
				min_coverage=300,
				max_coverage=600
			)
		)

		snpeff_results = gwf.target_from_template(
			name=f'snpeff_results_{SAMPLE_GROUP}_{SAMPLE_NAME.replace("-", "_")}',
			template=snpeff_result(
				effectsummary_file=annotation_frequencies.outputs['csv'],
				sitecount_file=site_count.outputs['sites'],
				output_directory=top_dir,
				sample_group=SAMPLE_GROUP,
				sample_name=SAMPLE_NAME)
		)

		snpeff_results_list.append(snpeff_results.outputs['tsv'])

		for sequence in sequences:
			REGION = sequence['sequence_name']
			REGION_LENGTH = sequence['sequence_length']

			snpgenie_pi = gwf.target_from_template(
				name=f'snpgenie_{SAMPLE_GROUP}_{SAMPLE_NAME.replace("-", "_")}_{REGION}',
				template=snpgenie_withinpool(
					reference_genome_file=REFERENCE,
					gtf_annotation_file=GTF,
					vcf_file=VCF,
					sample_group=SAMPLE_GROUP,
					sample_name=SAMPLE_NAME,
					region=REGION,
					region_length=REGION_LENGTH,
					output_directory=top_dir,
					min_allele_frequency=SNPGENIE_MINFREQ,
					sliding_window_size=SNPGENIE_SLIDINGWINDOW
				)
			)

			snpgenie_results_list.append(snpgenie_pi.outputs['plus']['summary'])
			snpgenie_results_list.append(snpgenie_pi.outputs['minus']['summary'])

	concat_snpeff_results = gwf.target_from_template(
		name=f'snpeff_concatenate_{SPECIES_NAME.replace(" ", "_")}',
		template=concatenate_snpeff_results(
			files=snpeff_results_list,
			output_name=f'{species_abbreviation(SPECIES_NAME)}.snpeff_results',
			output_directory=f'{top_dir}/snpEff'
		)
	)

	summarize_snpgenie = gwf.target_from_template(
		name=f'snpgenie_summarize_result_{SPECIES_NAME.replace(" ", "_")}',
		template=snpgenie_summarize_results_population(
			population_summary_files=snpgenie_results_list,
			output_directory=top_dir,
			species_name=SPECIES_NAME
		)
	)

	return gwf