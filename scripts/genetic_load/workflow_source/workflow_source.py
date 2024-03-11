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
	REFERENCE: str = config['reference_genome_path']
	GTF: str = config['gtf_annotation_file']
	WORK_DIR: str = config['working_directory_path']
	OUTPUT_DIR: str = config['output_directory_path']
	
	# --------------------------------------------------
	#                  Workflow
	# --------------------------------------------------
	
	gwf = Workflow(
		defaults={'account': ACCOUNT}
	)
	
	database_entry = gwf.target_from_template(
		name=f'{species_abbreviation(SPECIES_NAME)}_snpeff_database_entry',
		template=snpeff_database_build(
			gtf_annotation_file=GTF,
			reference_genome_file=REFERENCE,
			species_name=SPECIES_NAME
		)
	)
	
	return gwf