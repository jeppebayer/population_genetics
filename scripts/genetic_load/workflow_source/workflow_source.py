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

	database_entry = gwf.target_from_template(
		name=f'{species_abbreviation(SPECIES_NAME)}_snpeff_database_entry',
		template=snpeff_database_build(
			gtf_annotation_file=GTF,
			reference_genome_file=REFERENCE,
			species_name=SPECIES_NAME,
			snpeff_directory=SNPEFF_DIR
		)
	)
	
	variant_annotation = gwf.map(
		name=name_snpeff_annotation,
		template_func=snpeff_annotation,
		inputs=SAMPLES,
		extra={'snpeff_predictor_file': database_entry.outputs['predictor'],
		 	   'snpeff_config_file': f'{SNPEFF_DIR}/snpEff.config',
			   'output_directory': top_dir}
	)

	snpgenie_pi = gwf.map(
		name=name_snpgenie,
		template_func=snpgenie_withinpool,
		inputs=SAMPLES,
		extra={'reference_genome_file': REFERENCE,
			   'gtf_annotation_file': GTF,
			   'output_directory': top_dir,
			   'min_allele_frequency': SNPGENIE_MINFREQ,
			   'sliding_window_size': SNPGENIE_SLIDINGWINDOW}
	)

	return gwf