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
	SAMPLE_DIRS: list = config['sample_folder_list']
	WORK_DIR: str = config['working_directory_path']
	OUTPUT_DIR: str = config['output_directory_path']
	
	# --------------------------------------------------
	#                  Workflow
	# --------------------------------------------------
	
	gwf = Workflow(
		defaults={'account': ACCOUNT}
	)
	
	top_dir = f'{WORK_DIR}/{SPECIES_NAME.replace(" ", "_")}/mapping'
	
	sample_list = [get_sample_data(path) for path in SAMPLE_DIRS]
	print(sample_list)

	adapter_removal = gwf.map(
	)

	return gwf