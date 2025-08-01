#!/bin/env python3
from gwf import Workflow
from gwf.workflow import collect
from gwf.executors import Conda
import os, yaml, glob, sys, subprocess
from workflow_templates import *

def sfs_qc_workflow(configFile: str = glob.glob('*config.y*ml')[0]):
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
	CONDA_ENV: str = CONFIG['condaEnvironment01']
	SPECIES_NAME: str = CONFIG['speciesName']
	TAXONOMY: str | None = CONFIG['taxonomicGroup'].lower() if CONFIG['taxonomicGroup'] else None
	WORK_DIR: str =  CONFIG['workingDirectoryPath'][:len(CONFIG['workingDirectoryPath']) - 1] if CONFIG['workingDirectoryPath'].endswith('/') else CONFIG['workingDirectoryPath']
	OUTPUT_DIR: str | None = (CONFIG['outputDirectoryPath'][:len(CONFIG['outputDirectoryPath']) - 1] if CONFIG['outputDirectoryPath'].endswith('/') else CONFIG['outputDirectoryPath']) if CONFIG['outputDirectoryPath'] else None
	INTERGENIC_BED: str = CONFIG['intergenicBed']
	REPEATS_BED: str = CONFIG['repeatsBed']
	VCF_FILES: list = CONFIG['vcfFiles']
	
	# --------------------------------------------------
	#                  Workflow
	# --------------------------------------------------
	
	gwf = Workflow(
		defaults={'account': ACCOUNT},
        executor=Conda(CONDA_ENV)
	)

	topDir = f'{WORK_DIR}/{TAXONOMY.replace(" ", "_")}/{SPECIES_NAME.replace(" ", "_")}/vcf/sfs_qc' if TAXONOMY else f'{WORK_DIR}/{SPECIES_NAME.replace(" ", "_")}/vcf/sfs_qc'
	topOut = f'{OUTPUT_DIR}/vcf/sfs_qc/{TAXONOMY.replace(" ", "_")}/{SPECIES_NAME.replace(" ", "_")}' if TAXONOMY else f'{OUTPUT_DIR}/vcf/sfs_qc/{SPECIES_NAME.replace(" ", "_")}'
	
	for vcf in VCF_FILES:
		sampleName=os.path.basename(vcf).replace("-", "_")
		
		sfsCount = gwf.target_from_template(
			name=f'sfs_count_{sampleName}',
			template=sfs_count(
				vcfFile=vcf,
				intergenicBed=INTERGENIC_BED,
				repeatsBed=REPEATS_BED,
				outputDirectory=topDir,
				outputName=os.path.basename(os.path.splitext(os.path.splitext(vcf)[0])[0]) if vcf.endswith('.gz') else os.path.basename(os.path.splitext(vcf)[0]),
				environment=CONDA_ENV,
				group=sampleName
			)
		)

		sfsPlot = gwf.target_from_template(
			name=f'sfs_plot_{sampleName}',
			template=sfs_plot(
				sfsFile=sfsCount.outputs['sfs'],
				outputDirectory=topOut if OUTPUT_DIR else topDir,
				environment=CONDA_ENV,
				group=sampleName
			)
		)
	
	return gwf