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
	OUTGROUP_BEDFILE: int | None = OUTGROUP_SETTINGS['withinThresholdBedFile']
	VCF_GROUP_LIST: list = CONFIG['vcfGroupList']
	INTERGENIC_BED: str | None = CONFIG['intergenicBedFile']
	REPEATS_BED: str | None = CONFIG['repeatsBedFile']
	
	# --------------------------------------------------
	#                  Workflow
	# --------------------------------------------------
	
	gwf = Workflow(
		defaults={'account': ACCOUNT}
	)
	
	topDir = f'{WORK_DIR}/{TAXONOMY.replace(" ", "_")}/{SPECIES_NAME.replace(" ", "_")}/vcf' if TAXONOMY else f'{WORK_DIR}/{SPECIES_NAME.replace(" ", "_")}/vcf'
	topOut = f'{OUTPUT_DIR}/vcf/{TAXONOMY.replace(" ", "_")}/{SPECIES_NAME.replace(" ", "_")}' if TAXONOMY else f'{OUTPUT_DIR}/vcf/{SPECIES_NAME.replace(" ", "_")}'
	
	setupDict = {index: {'name': group['groupName'].lower().replace(' ', '_'),
				  												'vcfFile': os.path.abspath(group['vcfFile']),
																'minDP': group['minimumCoverage'],
																'maxDP': group['maximumCoverage'],
																'bedFile': os.path.abspath(group['withinThresholdBedFile'])}
				for index, group in enumerate(VCF_GROUP_LIST) if group['vcfFile']}

	bedFileList = []

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

	# if OUTGROUP_VCF and OUTGROUP_BEDFILE:
	# 	ancestralAlleleInferenceReference = gwf.target_from_template(
	# 		name=f'ancestral_allele_inference_reference_{os.path.basename(OUTGROUP_VCF)}',
	# 		template=ancestral_allele_inference_reference(
	# 			referenceGenome=indexReferenceGenome.outputs['symlink'],
	# 			outgroupSitesWithinCoverageBed=OUTGROUP_BEDFILE,
	# 			outputDirectory=topDir,
	# 			outputName=f'{speciesAbbreviation(SPECIES_NAME)}'
	# 		)
	# 	)

	# 	ancestralAlleleInferenceVariant = gwf.target_from_template(
	# 		name=f'ancestral_allele_inference_variant_{os.path.basename(OUTGROUP_VCF)}',
	# 		template=ancestral_allele_inference_variant(
	# 			outgroupVcf=OUTGROUP_VCF,
	# 			outgroupSitesWithinCoverageBed=OUTGROUP_BEDFILE,
	# 			outputDirectory=topDir,
	# 			outputName=f'{speciesAbbreviation(SPECIES_NAME)}'
	# 		)
	# 	)

	# 	ancestralAlleleInferenceMerge = gwf.target_from_template(
	# 		name=f'ancestral_allele_inference_merge_{os.path.basename(OUTGROUP_VCF)}',
	# 		template=ancestral_allele_inference_merge(
	# 			variantAnnotation=ancestralAlleleInferenceVariant.outputs['annotation'],
	# 			referenceAnnotation=ancestralAlleleInferenceReference.outputs['annotation'],
	# 			outputDirectory=topDir,
	# 			outputName=f'{speciesAbbreviation(SPECIES_NAME)}'
	# 		)
	# 	)

	for group in setupDict:
		samples = [{'sampleName': i} for i in sampleNamesVcf(setupDict[group]['vcfFile'])]

		normalizeVcf = gwf.target_from_template(
			name=f'normalize_vcf_{group}_{speciesAbbreviation(SPECIES_NAME)}_{setupDict[group]['name']}',
			template=normalize_vcf(
				vcfFile=setupDict[group]['vcfFile'],
				referenceGenomeFile=indexReferenceGenome.outputs['symlink'],
				outputDirectory=f'{topOut}/{setupDict[group]['name']}' if OUTPUT_DIR else f'{topDir}/normalized/{setupDict[group]['name']}'
			)
		)

		normalizedVcfStats = gwf.target_from_template(
			name=f'normalized_vcf_stats_{group}_{speciesAbbreviation(SPECIES_NAME)}_{setupDict[group]['name']}',
			template=vcf_stats(
				vcfFile=normalizeVcf.outputs['vcf'],
				referenceGenomeFile=indexReferenceGenome.outputs['symlink'],
				outputDirectory=f'{topOut}/{setupDict[group]['name']}' if OUTPUT_DIR else f'{topDir}/normalized/{setupDict[group]['name']}'
			)
		)

		filterVcf = gwf.target_from_template(
			name=f'filter_vcf_{group}_{speciesAbbreviation(SPECIES_NAME)}_{setupDict[group]['name']}',
			template=filter_vcf(
				vcfFile=normalizeVcf.outputs['vcf'],
				referenceGenomeFile=indexReferenceGenome.outputs['symlink'],
				depthThresholdBed=setupDict[group]['bedFile'],
				minDepth=setupDict[group]['minDP'],
				maxDepth=setupDict[group]['maxDP'],
				outputDirectory=f'{topOut}/{setupDict[group]['name']}' if OUTPUT_DIR else f'{topDir}/filtered/{setupDict[group]['name']}'
			)
		)

		filteredVcfStats = gwf.target_from_template(
			name=f'filtered_vcf_stats_{group}_{speciesAbbreviation(SPECIES_NAME)}_{setupDict[group]['name']}',
			template=vcf_stats(
				vcfFile=filterVcf.outputs['vcf'],
				referenceGenomeFile=indexReferenceGenome.outputs['symlink'],
				outputDirectory=f'{topOut}/{setupDict[group]['name']}' if OUTPUT_DIR else f'{topDir}/filtered/{setupDict[group]['name']}'
			)
		)

		# variableSiteCount = gwf.target_from_template(
		# 	name=f'variable_site_count_{group}_{speciesAbbreviation(SPECIES_NAME)}_{setupDict[group]['name']}',
		# 	template=variable_site_count(
		# 		vcfFileBefore=setupDict[group]['vcfFile'],
		# 		vcfFileAfter=filterVcf.outputs['vcf'],
		# 		outputDirectory=topOut if OUTPUT_DIR else topDir
		# 	)
		# )

		vcfToBed = gwf.target_from_template(
			name=f'vcf_to_bed_{os.path.basename(filterVcf.outputs['vcf']).replace("-", "_")}',
			template=vcf_to_bed(
				vcfFile=filterVcf.outputs['vcf'],
				outputDirectory=f'{topDir}/outgroups'
			)
		)

		bedFileList.append(vcfToBed.outputs['bed'])

		if GENOME_ANNOTATION:
			snpeffAnnotation = gwf.target_from_template(
				name=f'snpeff_annotation_{group}_{speciesAbbreviation(SPECIES_NAME)}_{setupDict[group]['name']}',
				template=snpeff_annotation(
					vcfFile=filterVcf.outputs['vcf'],
					snpeffPredictorFile=snpeffBuildDatabase.outputs['predictor'],
					vcfOutputDirectory=f'{topOut}/{setupDict[group]['name']}' if OUTPUT_DIR else f'{topDir}/filtered/{setupDict[group]['name']}',
					snpeffOutputDirectory=topDir
				)
			)
		
		if INTERGENIC_BED and REPEATS_BED:
			sfsNeutral = gwf.map(
				name=f'sfs_neutral_{group}_{speciesAbbreviation(SPECIES_NAME)}_{setupDict[group]['name']}',
				template_func=sfs_neutral,
				inputs=samples,
				extra={'vcfFile': filterVcf.outputs['vcf'],
		   			   'intergenicBed': INTERGENIC_BED,
					   'repeatsBed': REPEATS_BED,
					   'outputDirectory': f'{topDir}/filtered/{setupDict[group]['name']}/sfs/tmp'}
			)

			sfsMerge = gwf.target_from_template(
				name=f'merge_sfs_{group}_{speciesAbbreviation(SPECIES_NAME)}_{setupDict[group]['name']}',
				template=sfs_merge(
					sfsFiles=collect(sfsNeutral.outputs, ['sfs'])['sfss'],
					outputName=os.path.basename(os.path.splitext(os.path.splitext(filterVcf.outputs['vcf'])[0])[0]) if filterVcf.outputs['vcf'].endswith('.gz') else os.path.basename(os.path.splitext(filterVcf.outputs['vcf'])[0]),
					outputDirectory=f'{topDir}/filtered/{setupDict[group]['name']}/sfs'
				)
			)

			sfsPlot = gwf.target_from_template(
				name=f'plot_sfs_{group}_{speciesAbbreviation(SPECIES_NAME)}_{setupDict[group]['name']}',
				template=sfs_plot(
					sfsFile=sfsMerge.outputs['sfs'],
					outputDirectory=f'{topOut}/{setupDict[group]['name']}/vcfStats' if OUTPUT_DIR else f'{topDir}/filtered/{setupDict[group]['name']}/sfs'
				)
			)

		# 	if OUTGROUP_VCF and OUTGROUP_BEDFILE:
		# 		updateAncetralAlleleInformation = gwf.target_from_template(
		# 			name=f'update_ancestral_allele_information_{group}_{speciesAbbreviation(SPECIES_NAME)}_{setupDict[group]['name']}',
		# 			template=update_ancetral_allele_information(
		# 				vcfFile=snpeffAnnotation.outputs['vcf'],
		# 				ancestralAnnotationFile=ancestralAlleleInferenceMerge.outputs['annotation'],
		# 				outputDirectory=f'{topOut}/{setupDict[group]['name']}' if OUTPUT_DIR else f'{topDir}/filtered/{setupDict[group]['name']}'
		# 			)
		# 		)

	mergeBedFiles = gwf.target_from_template(
		name=f'merge_bed_files',
		template=merge_bed_files(
			bedFiles=bedFileList,
			outputName=f'{speciesAbbreviation(SPECIES_NAME)}.merge',
			outputDirectory=os.path.dirname(bedFileList[0])
		)
	)

	print(f'Intermediary files will be place at: {topDir}/')
	print(f'Output files will be placed at: {topOut if OUTPUT_DIR else topDir}/')
	
	return gwf