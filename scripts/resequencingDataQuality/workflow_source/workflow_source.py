#!/bin/env python3
from gwf import Workflow
from gwf.workflow import collect
from gwf.executors import Conda
import os, yaml, glob, sys
from workflow_templates import *

def resequencing_data_quality_workflow(config_file: str = glob.glob('*config.y*ml')[0]):
	"""
	
	
	:param str config_file:
		Configuration file containing pre-defined set of variables
	"""
	# --------------------------------------------------
	#                  Configuration
	# --------------------------------------------------
	
	CONFIG = yaml.safe_load(open(config_file))
	ACCOUNT: str = CONFIG['account']
	ENVIRONMENT: str = CONFIG['condaEnvironment']
	WORK_DIR: str = CONFIG['outputDirectory']
	SETUP: list = CONFIG['setup']
	
	# --------------------------------------------------
	#                  Workflow
	# --------------------------------------------------
	
	gwf = Workflow(
		defaults={'account': ACCOUNT},
		executor=Conda(ENVIRONMENT)
	)

	summaryStatsFiles = []

	for species in SETUP:
		SPECIES_NAME: str = species['speciesName']
		TAXONOMY: str | None = species['taxonomicGroup'].lower().replace(' ', '_') if species['taxonomicGroup'] else None
		PLATFORM: str = species['sequencingPlatform']
		ALIGNMENT_FILES: list = species['alignmentFiles']

		topDir = f'{WORK_DIR}/{TAXONOMY}/{SPECIES_NAME.replace(" ", "_")}/alignment' if TAXONOMY else f'{WORK_DIR}/{SPECIES_NAME.replace(" ", "_")}/alignment'

		if PLATFORM.lower().startswith('novaseq'):
			opticalDistance = 2500
		elif PLATFORM.lower().startswith('hiseq'):
			opticalDistance = 100
		else:
			opticalDistance = 0

		for alignmentFile in ALIGNMENT_FILES:
			sampleName = os.path.basename(os.path.dirname(alignmentFile))
			currentDir = f'{topDir}/{sampleName}'

			markDuplicates = gwf.target_from_template(
				name=f'mark_duplicates_{sampleName.replace("-", "_")}',
				template=mark_duplicates(
					alignmentFile=alignmentFile,
					sampleName=sampleName,
					outputDirectory=currentDir,
					environment=ENVIRONMENT,
					opticalDistance=opticalDistance,
					group=species_abbreviation(SPECIES_NAME),
				)
			)

			alignmentStats = gwf.target_from_template(
				name=f'alignment_stat_{sampleName.replace("-", "_")}',
				template=alignment_stats(
					alignmentFile=markDuplicates.outputs['markdup'],
					outputDirectory=currentDir,
					environment=ENVIRONMENT,
					group=species_abbreviation(SPECIES_NAME)
				)
			)

			summaryStats = gwf.target_from_template(
				name=f'summary_stats_{sampleName.replace("-", "_")}',
				template=summary_stats(
					markdupstatsFile=markDuplicates.outputs['stats'],
					statsFile=alignmentStats.outputs['stats'],
					flagstatFile=alignmentStats.outputs['flagstat'],
					sampleName=sampleName,
					outputDirectory=currentDir,
					environment=ENVIRONMENT,
					group=species_abbreviation(SPECIES_NAME)
				)
			)

			summaryStatsFiles.append(summaryStats.outputs['stats'])

	collateStats = gwf.target_from_template(
		name=f'collate_summary_stats',
		template=collate_stats(
			summaryStatsFiles=summaryStatsFiles,
			outputDirectory=WORK_DIR,
			environment=ENVIRONMENT
		)
	)

	return gwf