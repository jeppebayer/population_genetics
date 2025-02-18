#!/bin/env python3
from gwf import Workflow
from gwf.workflow import collect
import os, yaml, glob, sys
from workflow_templates import *

def vcfFilterAndAnnotation_workflow(configFile: str = glob.glob('*config.y*ml')[0]):
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
	TAXONOMY: str | None = CONFIG['taxonomicGroup'].lower() if CONFIG['taxonomicGroup'] else None
	SPECIES_NAME: str = CONFIG['speciesName']
	WORK_DIR: str =  CONFIG['workingDirectoryPath'][:len(CONFIG['workingDirectoryPath']) - 1] if CONFIG['workingDirectoryPath'].endswith('/') else CONFIG['workingDirectoryPath']
	OUTPUT_DIR: str | None = (CONFIG['outputDirectoryPath'][:len(CONFIG['outputDirectoryPath']) - 1] if CONFIG['outputDirectoryPath'].endswith('/') else CONFIG['outputDirectoryPath']) if CONFIG['outputDirectoryPath'] else None
	OUTGROUP_VCF: str | None = CONFIG['outgroupVcfFile'] if CONFIG['outgroupVcfFile'] else None
	VCF_GROUP_LIST: list = CONFIG['vcfGroupList']
	
	# --------------------------------------------------
	#                  Workflow
	# --------------------------------------------------
	
	gwf = Workflow(
		defaults={'account': ACCOUNT}
	)
	
	topDir = f'{WORK_DIR}/{TAXONOMY.replace(" ", "_")}/{SPECIES_NAME.replace(" ", "_")}/vcf' if TAXONOMY else f'{WORK_DIR}/{SPECIES_NAME.replace(" ", "_")}/vcf'
	topOut = f'{OUTPUT_DIR}/vcf/{TAXONOMY.replace(" ", "_")}/{SPECIES_NAME.replace(" ", "_")}' if TAXONOMY else f'{OUTPUT_DIR}/vcf/{SPECIES_NAME.replace(" ", "_")}'
	
	setupDict = {group['groupName'].lower().replace(' ', '_'): {'name': group['groupName'].lower().replace(' ', '_'),
				  												'vcfFile': os.path.abspath(group['vcfFile']),
																'minDP': group['minimumCoverage'],
																'maxDP': group['maximumCoverage'],
																'bedFile': os.path.abspath(group['withinThresholdBedFile'])}
				for group in VCF_GROUP_LIST if group['vcfFile']}

	# indexReferenceGenome = gwf.target_from_template(
	# 	name=f'index_reference_genome_{SPECIES_NAME.replace(" ", "_")}',
	# 	template=index_reference_genome(
	# 		referenceGenomeFile=REFERENCE_GENOME,
	# 		outputDirectory=topDir
	# 	)
	# )

	for group in setupDict:
		filterVcf = gwf.target_from_template(
			name=f'filter_vcf_{speciesAbbreviation(SPECIES_NAME)}_{setupDict[group]['name']}',
			template=filter_vcf(
				vcfFile=setupDict[group]['vcfFile'],
				depthThresholdBed=setupDict[group]['bedFile'],
				minDepth=setupDict[group]['minDP'],
				maxDepth=setupDict[group]['maxDP'],
				outputDirectory=f'{topOut}/{setupDict[group]['name']}' if OUTPUT_DIR else f'{topDir}/filtered/{setupDict[group]['name']}'
			)
		)

	print(f'Intermediary files will be place at: {topDir}/')
	print(f'Output files will be placed at: {topOut if OUTPUT_DIR else topDir}/')
	
	return gwf