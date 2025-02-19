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
	REFERENCE_GENOME: str = CONFIG['referenceGenome']
	GENOME_ANNOTATION: str = CONFIG['gtfAnnotation']
	OUTGROUP_SETTINGS: dict = CONFIG['outgroupSettings']
	OUTGROUP_VCF: str | None = OUTGROUP_SETTINGS['vcfFile']
	OUTGROUP_MINDP: int | None = OUTGROUP_SETTINGS['minimumCoverage']
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

	indexReferenceGenome = gwf.target_from_template(
		name=f'index_reference_genome_{SPECIES_NAME.replace(" ", "_")}',
		template=index_reference_genome(
			referenceGenomeFile=REFERENCE_GENOME,
			outputDirectory=topDir
		)
	)

	if GENOME_ANNOTATION:
		snpeffBuildDatabase = gwf.target_from_template(
			name=f'snpeff_build_database_{speciesAbbreviation(SPECIES_NAME)}',
			template=snpeff_build_database(
				referenceGenome=indexReferenceGenome.outputs['symlink'],
				gtfAnnotation=GENOME_ANNOTATION,
				outputDirectory=topDir,
				speciesName=SPECIES_NAME
			)
		)

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

		variableSiteCount = gwf.target_from_template(
			name=f'variable_site_count_{speciesAbbreviation(SPECIES_NAME)}_{setupDict[group]['name']}',
			template=variable_site_count(
				vcfFileBefore=setupDict[group]['vcfFile'],
				vcfFileAfter=filterVcf.outputs['vcf'],
				outputDirectory=topOut if OUTPUT_DIR else topDir,
				outputName=f'{speciesAbbreviation(SPECIES_NAME)}.{setupDict[group]['name']}.variable.sitetable'
			)
		)

		if GENOME_ANNOTATION:
			snpeffAnnotation = gwf.target_from_template(
				name=f'snpeff_annotation_{speciesAbbreviation(SPECIES_NAME)}_{setupDict[group]['name']}',
				template=snpeff_annotation(
					vcfFile=filterVcf.outputs['vcf'],
					snpeffPredictorFile=snpeffBuildDatabase.outputs['predictor'],
					outputDirectory=topOut
				)
			)



	print(f'Intermediary files will be place at: {topDir}/')
	print(f'Output files will be placed at: {topOut if OUTPUT_DIR else topDir}/')
	
	return gwf