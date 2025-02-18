#!/bin/env python3
from gwf import Workflow
from gwf.workflow import collect
import os, yaml, glob, sys
from workflow_templates import *

def depthFractions_workflow(configFile: str = glob.glob('*config.y*ml')[0]):
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
	SAMPLE_LIST: list = CONFIG['sampleList']

	# --------------------------------------------------
	#                  Workflow
	# --------------------------------------------------
	
	gwf = Workflow(
		defaults={'account': ACCOUNT}
	)
	
	topDir = f'{WORK_DIR}/{TAXONOMY.replace(" ", "_")}/{SPECIES_NAME.replace(" ", "_")}/mapping' if TAXONOMY else f'{WORK_DIR}/{SPECIES_NAME.replace(" ", "_")}/mapping'
	topOut = f'{OUTPUT_DIR}/alignments/{TAXONOMY.replace(" ", "_")}/{SPECIES_NAME.replace(" ", "_")}' if TAXONOMY else f'{OUTPUT_DIR}/alignments/{SPECIES_NAME.replace(" ", "_")}'
	
	for GROUP in SAMPLE_LIST:
		if not GROUP['bamFileList']:
			continue
		GROUP_NAME: str = GROUP['groupName'].lower().replace(" ", "_")
		FILE_LIST: list = GROUP['bamFileList']
		mapListOfDict = [{'bamFile': file, 'sampleName': os.path.basename(file).split('.')[0]} for file in FILE_LIST]

		depthFractions = gwf.map(
			name=name_depth_fractions,
			template_func=depth_fractions,
			inputs=mapListOfDict,
			extra={
				'lowerThreshold': GROUP['lowerCoverageThreshold'],
				'outputDirectory': f'{topDir}/{GROUP_NAME}/depth'
				}
		)

		concatDepthFractionsTsv = gwf.target_from_template(
			name=f'concat_depth_fractions_tsv_{GROUP_NAME}',
			template=concat_depth_fractions_tsv(
				tsvFiles=collect(depthFractions.outputs, ['tsv'])['tsvs'],
				outputDirectory=f'{topOut}/{GROUP_NAME}/depth' if OUTPUT_DIR else f'{topDir}/{GROUP_NAME}/depth',
				outputName=f'{speciesAbbreviation(SPECIES_NAME)}.{GROUP_NAME}.depthFractions'
			)
		)

		concatDepthFractionsPlots = gwf.target_from_template(
			name=f'concat_depth_fractions_plots_{GROUP_NAME}',
			template=concat_depth_fractions_plots(
				plotFiles=collect(depthFractions.outputs, ['png'])['pngs'],
				outputDirectory=f'{topOut}/{GROUP_NAME}/depth' if OUTPUT_DIR else f'{topDir}/{GROUP_NAME}/depth',
				outputName=f'{speciesAbbreviation(SPECIES_NAME)}.{GROUP_NAME}.depthFractions'
			)
		)

	return gwf