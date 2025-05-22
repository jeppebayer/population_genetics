#!/bin/env python3
from gwf import Workflow
from gwf.workflow import collect
import os, yaml, glob, sys
from workflow_templates import *

def ancestral_allele_workflow(configFile: str = glob.glob('*config.y*ml')[0]):
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
	VCF_BED: str = CONFIG['ingroupsBedFile']
	REFERENCE_GENOME: str = CONFIG['referenceGenomeFile']
	OUTGROUP_SETUP: list = CONFIG['outgroupSetup']
	INGROUP_SETUP: list = CONFIG['ingroupSetup']
	FREEBAYES_SETTINGS: dict = CONFIG['freebayesSettings']
	FREEBAYES_PLOIDY: int | None = FREEBAYES_SETTINGS['samplePloidy'] if FREEBAYES_SETTINGS['samplePloidy'] else 100
	FREEBAYES_BESTN: int | None = FREEBAYES_SETTINGS['bestNAlleles'] if FREEBAYES_SETTINGS['bestNAlleles'] else 3
	FREEBAYES_MINALTFRC: float | int | None = FREEBAYES_SETTINGS['minAlternateFraction'] if FREEBAYES_SETTINGS['minAlternateFraction'] else 0
	FREEBAYES_MINALTCNT: int | None = FREEBAYES_SETTINGS['minAlternateCount'] if FREEBAYES_SETTINGS['minAlternateCount'] else 2
	FREEBAYES_MEM: int | None = FREEBAYES_SETTINGS['memory'] if FREEBAYES_SETTINGS['memory'] else 50
	FREEBAYES_TIME: str | None = FREEBAYES_SETTINGS['time'] if FREEBAYES_SETTINGS['time'] else '24:00:00'
	
	# --------------------------------------------------
	#                  Workflow
	# --------------------------------------------------
	
	gwf = Workflow(
		defaults={'account': ACCOUNT}
	)
	
	topDir = f'{WORK_DIR}/{TAXONOMY.replace(" ", "_")}/{SPECIES_NAME.replace(" ", "_")}/vcf' if TAXONOMY else f'{WORK_DIR}/{SPECIES_NAME.replace(" ", "_")}/vcf'
	topOut = f'{OUTPUT_DIR}/vcf/{TAXONOMY.replace(" ", "_")}/{SPECIES_NAME.replace(" ", "_")}' if TAXONOMY else f'{OUTPUT_DIR}/vcf/{SPECIES_NAME.replace(" ", "_")}'
	
	if not VCF_BED:
		print(f"A bed file matching the variant position of the vcf file must be supplied.\nIf you do not have one it can be created running the following command:\n\tbcftools query --format '%CHROM\\t%POS0\\t%END' <file.vcf> | bedtools merge -i - > <file.vcf.bed>")
		sys.exit()
	
	partitions = partition_bed(parse_bed(VCF_BED))
	outgroupVcfList = []

	indexReferenceGenome = gwf.target_from_template(
		name=f'index_reference_genome',
		template=index_reference_genome(
			referenceGenomeFile=REFERENCE_GENOME,
			outputDirectory=topDir
		)
	)

	for sample in OUTGROUP_SETUP:
		sampleName = sample['sampleName'].replace(' ',' _')
		freebayesPartitionOutgroup = gwf.map(
			name=name_freebayes_partition_outgroup,
			template_func=freebayes_partition_outgroup,
			inputs=partitions,
			extra={
				'referenceGenomeFile': indexReferenceGenome.outputs['symlink'],
				'bamFile': sample['bamFile'],
				'outputDirectory': topDir,
				'sampleName': sampleName,
				'bedFile': VCF_BED,
				'ploidy': FREEBAYES_PLOIDY,
				'bestNAlleles': FREEBAYES_BESTN,
				'minAlternateFraction': FREEBAYES_MINALTFRC,
				'minAlternateCount': FREEBAYES_MINALTCNT,
				'memory': FREEBAYES_MEM,
				'time': FREEBAYES_TIME
			}
		)

		concatenateFreebayesPartitionOutgroup = gwf.target_from_template(
			name=f'concatenate_outgroup_vcf_{sampleName.replace("-", "_")}',
			template=concat_vcf(
				files=collect(freebayesPartitionOutgroup.outputs, ['vcf'])['vcfs'],
				outputName=f'{sample['sampleName']}.outgroups.freebayes',
				outputDirectory=f'{topOut}/outgroups/{sampleName}' if OUTPUT_DIR else f'{topDir}/outgroups/{sampleName}',
				compress=True
			)
		)

		outgroupVcfList.append(concatenateFreebayesPartitionOutgroup.outputs['concat_file'])

	if len(outgroupVcfList) > 1:
		mergeOutgroups = gwf.target_from_template(
			name=f'outgroups_vcf_merge',
			template=merge_vcf(
				vcfFiles=outgroupVcfList,
				outputName=f'{speciesAbbreviation(SPECIES_NAME)}.outgroups.freebayes',
				outputDirectory=f'{topDir}/outgroups'
			)
		)

	else:
		cpRenameOutgroups = gwf.target_from_template(
			name=f'outgroups_vcf_cp_rename',
			template=cp_rename_vcf(
				vcfFile=outgroupVcfList[0],
				outputName=f'{speciesAbbreviation(SPECIES_NAME)}.outgroups.freebayes',
				outputDirectory=f'{topDir}/outgroups'
			)
		)

	normalizeVcf = gwf.target_from_template(
		name=f'normalize_vcf_outgroups',
		template=normalize_vcf(
			vcfFile=mergeOutgroups.outputs['vcf'] if len(outgroupVcfList) > 1 else cpRenameOutgroups.outputs['vcf'],
			referenceGenomeFile=indexReferenceGenome.outputs['symlink'],
			outputDirectory=f'{topDir}/outgroups'
		)
	)

	vcfStatsOutgroups = gwf.target_from_template(
		name=f'vcf_stats_outgroups',
		template=vcf_stats(
			vcfFile=normalizeVcf.outputs['vcf'],
			referenceGenomeFile=indexReferenceGenome.outputs['symlink'],
			outputDirectory=f'{topDir}/outgroups/vcfStats'
		)
	)

	for num, group in enumerate(INGROUP_SETUP):
		groupName = group['groupName'].lower().replace(' ', '_')
		mergeIngroupOutgroup = gwf.target_from_template(
			name=f'merge_ingroup_outgroup_{groupName}_{num}',
			template=merge_ingroup_outgroup(
				ingroupVcfFile=group['vcfFile'],
				outgroupVcfFile=normalizeVcf.outputs['vcf'],
				outputDirectory=f'{topDir}/withOutgroups/{groupName}'
			)
		)

	return gwf