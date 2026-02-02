#!/bin/env python3
from gwf import Workflow
from gwf.workflow import collect
from gwf.executors import Conda
import os, sys, glob, yaml
from workflow_templates import *

def geneticLoad_workflow(configFile: str = glob.glob('*config.y*ml')[0]):
	"""
	Workflow: description
	
	:param str configFile:
		Configuration file containing pre-defined set of variables
	"""
	# --------------------------------------------------
	#                  Configuration
	# --------------------------------------------------
	
	CONFIG = yaml.safe_load(open(configFile))
	ACCOUNT: str = CONFIG['account']
	CONDA_ENV_01: str = CONFIG['condaEnvironment01']
	CONDA_ENV_02: str = CONFIG['condaEnvironment02'] # Environment containing snpgenie
	TAXONOMY: str | None = CONFIG['taxonomicGroup'].lower() if CONFIG['taxonomicGroup'] else None
	SPECIES_NAME: str = CONFIG['speciesName']
	WORK_DIR: str =  CONFIG['workingDirectoryPath'][:len(CONFIG['workingDirectoryPath']) - 1] if CONFIG['workingDirectoryPath'].endswith('/') else CONFIG['workingDirectoryPath']
	OUTPUT_DIR: str | None = (CONFIG['outputDirectoryPath'][:len(CONFIG['outputDirectoryPath']) - 1] if CONFIG['outputDirectoryPath'].endswith('/') else CONFIG['outputDirectoryPath']) if CONFIG['outputDirectoryPath'] else None
	REFERENCE_GENOME: str = CONFIG['referenceGenome']
	ANNOTATION_FILE: str = CONFIG['annotationFile']
	SAMPLE_SETUP: list = CONFIG['sampleSetup']

	# --------------------------------------------------
	#                  Workflow
	# --------------------------------------------------
	
	gwf = Workflow(
		defaults={'account': ACCOUNT},
		executor=Conda(CONDA_ENV_01)
	)
	
	sequenceNameList = sequenceNamesFasta(REFERENCE_GENOME)

	topDir = f'{WORK_DIR}/{TAXONOMY.replace(" ", "_")}/{SPECIES_NAME.replace(" ", "_")}/genetic_load' if TAXONOMY else f'{WORK_DIR}/{SPECIES_NAME.replace(" ", "_")}/genetic_load'
	topOut = f'{OUTPUT_DIR}/genetic_load/{TAXONOMY.replace(" ", "_")}/{SPECIES_NAME.replace(" ", "_")}' if TAXONOMY else f'{OUTPUT_DIR}/genetic_load/{SPECIES_NAME.replace(" ", "_")}'
	
	for sample in SAMPLE_SETUP:
		GROUP_NAME: str = sample['groupName'].lower().replace(' ', '_')
		VCF_FILE: str = sample['vcfFile']

		sampleNameList = sampleNamesVCF(VCF_FILE)

		snpeffResults = gwf.target_from_template(
			name=f'snpeff_results_{GROUP_NAME}',
			template=snpeff_results(
				vcfFile=VCF_FILE,
				gtfAnnotationFile=ANNOTATION_FILE,
				speciesName=SPECIES_NAME,
				outputDirectory=f'{topDir}/{GROUP_NAME}',
				environment=CONDA_ENV_01,
				group=GROUP_NAME
			)
		)

		vcfReformat = gwf.target_from_template(
			name=f'vcf_reformat_{GROUP_NAME}',
			template=vcf_reformat(
				vcfFile=VCF_FILE,
				outputDirectory=F'{topDir}/{GROUP_NAME}',
				environment=CONDA_ENV_01,
				group=GROUP_NAME
			)
		)

		# for sampleName in sampleNameList:
		# 	for sequence in sequenceNameList:
		# 		snpgenieSequence = gwf.target_from_template(
		# 			name=f'snpgenie_{sequence.replace("-", "_")}_{sampleName.replace("-", "_")}_{GROUP_NAME}',
		# 			template=snpgenie_sequence(
		# 				referenceGenome=REFERENCE_GENOME,
		# 				gtfFile=ANNOTATION_FILE,
		# 				vcfFile=vcfReformat.outputs['vcf'],
		# 				sampleName=sampleName,
		# 				region=sequence,
		# 				outputDirectory=f'{topDir}/{GROUP_NAME}',
		# 				environment=CONDA_ENV_02,
		# 				group=GROUP_NAME
		# 			)
		# 		)
	
	print(f'Intermediary files will be place at: {topDir}')
	print(f'Output files will be placed at: {topOut if OUTPUT_DIR else topDir}')

	return gwf