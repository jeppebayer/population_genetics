#!/bin/env python3
from gwf import Workflow
from gwf.workflow import collect
import os, yaml, glob, sys
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
	
	config = yaml.safe_load(open(config_file))
	ACCOUNT: str = config['account']
	SPECIES_NAME: str = config['species_name']
	VCF: str = config['vcf_file']
	REFERENCE: str = config['reference_genome_path']
	GTF: str = config['gtf_annotation_file']
	WORK_DIR: str = config['working_directory_path']
	OUTPUT_DIR: str = config['output_directory_path']
	
	snpeff_directory = f'{os.path.dirname(os.path.realpath(__file__))}/software/snpeff'

	# --------------------------------------------------
	#                  Workflow
	# --------------------------------------------------
	
	gwf = Workflow(
		defaults={'account': ACCOUNT}
	)

	top_dir = f'{WORK_DIR}/{SPECIES_NAME.replace(" ", "_")}/genetic_load'
	if not OUTPUT_DIR:
		OUTPUT_DIR = top_dir

	database_entry = gwf.target_from_template(
		name=f'{species_abbreviation(SPECIES_NAME)}_snpeff_database_entry',
		template=snpeff_database_build(
			gtf_annotation_file=GTF,
			reference_genome_file=REFERENCE,
			species_name=SPECIES_NAME,
			snpeff_directory=snpeff_directory
		)
	)
	
	variant_annotation = gwf.target_from_template(
		name=f'{species_abbreviation(SPECIES_NAME)}_snpeff_annotation',
		template=snpeff_annotation(
			vcf_file=VCF,
			snpeff_predictor_file=database_entry.outputs['predictor'],
			snpeff_config_file=f'{snpeff_directory}/snpEff.config',
			output_directory=top_dir,
			species_name=SPECIES_NAME
		)
	)

	return gwf