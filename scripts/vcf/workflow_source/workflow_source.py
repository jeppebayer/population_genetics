#!/bin/env python3
from gwf import Workflow
from gwf.workflow import collect
import os, yaml, glob, sys
from workflow_templates import *

def freebayes_population_set_workflow(configFile: str = glob.glob('*config.y*ml')[0]):
	"""
	Workflow: Create :format:`VCF` file for each sample in configuration.
	
	:param str config_file:
		Configuration file containing pre-defined set of variables
	"""
	# --------------------------------------------------
	#                  Configuration
	# --------------------------------------------------
	
	CONFIG = yaml.safe_load(open(configFile))
	ACCOUNT: str = CONFIG['account']
	TAXONOMY: str | None = CONFIG['taxonomicGroup'].lower() if CONFIG['taxonomicGroup'] else None
	SPECIES_NAME: str = CONFIG['speciesName']
	REFERENCE_GENOME: str = CONFIG['referenceGenomePath']
	WORK_DIR: str = CONFIG['workingDirectoryPath'][:len(CONFIG['workingDirectoryPath']) - 1] if CONFIG['workingDirectoryPath'].endswith('/') else CONFIG['workingDirectoryPath']
	OUTPUT_DIR: str | None = (CONFIG['outputDirectoryPath'][:len(CONFIG['outputDirectoryPath']) - 1] if CONFIG['outputDirectoryPath'].endswith('/') else CONFIG['outputDirectoryPath']) if CONFIG['outputDirectoryPath'] else None
	PARTITION_SIZE: int = CONFIG['partitionSize'] if CONFIG['partitionSize'] else 500000
	FREEBAYES_SETTINGS: dict = CONFIG['freebayesSettings']
	FREEBAYES_MODE: int = FREEBAYES_SETTINGS['mode']
	FREEBAYES_PLOIDY: int | None = FREEBAYES_SETTINGS['samplePloidy'] if FREEBAYES_SETTINGS['samplePloidy'] else 100
	FREEBAYES_BESTN: int | None = FREEBAYES_SETTINGS['bestNAlleles'] if FREEBAYES_SETTINGS['bestNAlleles'] else 3
	FREEBAYES_MINALTFRC: float | int | None = FREEBAYES_SETTINGS['minAlternateFraction'] if FREEBAYES_SETTINGS['minAlternateFraction'] else 0
	FREEBAYES_MINALTCNT: int | None = FREEBAYES_SETTINGS['minAlternateCount'] if FREEBAYES_SETTINGS['minAlternateCount'] else 2
	VCF_MEM: int | None = FREEBAYES_SETTINGS['memory'] if FREEBAYES_SETTINGS['memory'] else 80
	VCF_TIME: str | None = str(FREEBAYES_SETTINGS['time']) if FREEBAYES_SETTINGS['time'] else '48:00:00'
	FILTERING: dict = CONFIG['filtering']
	INGROUP_FILTERING_MINDP: int | None = FILTERING['ingroupMinimumDepth'] if FILTERING['ingroupMinimumDepth'] else 200
	OUTGROUP_FILTERING_MINDP: int | None = FILTERING['outgroupMinimumDepth'] if FILTERING['outgroupMinimumDepth'] else 5
	SAMPLE_LIST: list = CONFIG['sampleList']
	INTERGENIC_BED: str = CONFIG['intergenicBedFile']
	REPEATS_BED: str | None = CONFIG['repeatsBedFile']
	BATCHSETTINGS: dict = CONFIG['batchSettings']
	NBATCHES: int = BATCHSETTINGS['numberOfBatches'] if BATCHSETTINGS['numberOfBatches'] else 0
	BATCHNR: int | None = BATCHSETTINGS['currentBatchNumber'] if BATCHSETTINGS['currentBatchNumber'] else None

	# --------------------------------------------------
	#                  Workflow
	# --------------------------------------------------
	
	gwf = Workflow(
		defaults={'account': ACCOUNT}
	)
	
	if os.path.exists(f'reference_partitions.{PARTITION_SIZE}bp.txt') and os.path.exists('reference_sequences.txt'):
		# If files exists reads data directly from files
		# Loads reference genome partitioning
		with open(f'reference_partitions.{PARTITION_SIZE}bp.txt', 'r') as infile:
			partitions = [{'num': entry.split(sep='\t')[0].strip(), 'region': entry.split(sep='\t')[1].strip(), 'start': entry.split(sep='\t')[2].strip(), 'end': entry.split(sep='\t')[3].strip()} for entry in infile]
			nPadding = len(str(sum(1 for line in partitions)))
		# Loads list of contigs in reference genome
		with open('reference_sequences.txt', 'r') as infile:
			contigs = [{'contig': entry.split(sep='\t')[0].strip()} for entry in infile]
	else:
		# If files don't exist, generate data and write files
		sequences = parse_fasta(REFERENCE_GENOME)
		with open(f'reference_sequences.txt', 'w') as outfile:
			outfile.write('\n'.join('\t'.join(str(i) for i in entry.values()) for entry in sequences))
		# Partitions reference genome
		nPadding = padding_calculator(parseFasta=sequences, size=PARTITION_SIZE)
		partitions = partition_chrom(parseFasta=sequences, size=PARTITION_SIZE, nPad=nPadding)
		with open(f'reference_partitions.{PARTITION_SIZE}bp.txt', 'w') as outfile:
			outfile.write('\n'.join('\t'.join(str(i) for i in entry.values()) for entry in partitions))
		# Creates list of contigs in reference genome
		contigs = [{'contig': contig['sequence_name']} for contig in sequences]

	# When the argument line for any command becomes too long it cannot be executed. This can become an issue when jobs have too many dependencies.
	# For this workflow in can occur when cancatenating massively parallellised task, eg. create the VCF parts.
	# To alleviate the issue, multiple concatenation jobs will be created so that no concatenation has more than 5000 dependencies.
	nPartitions = len(partitions)
	segmentSize = 5000
	nSegments = int(round(nPartitions / segmentSize, 0) + 1) if (nPartitions / segmentSize > round(nPartitions / segmentSize, 0)) else int(round(nPartitions / segmentSize, 0))

	topDir = f'{WORK_DIR}/{TAXONOMY.replace(" ", "_")}/{SPECIES_NAME.replace(" ", "_")}/vcf' if TAXONOMY else f'{WORK_DIR}/{SPECIES_NAME.replace(" ", "_")}/vcf'
	topOut = f'{OUTPUT_DIR}/vcf/{TAXONOMY.replace(" ", "_")}/{SPECIES_NAME.replace(" ", "_")}' if TAXONOMY else f'{OUTPUT_DIR}/vcf/{SPECIES_NAME.replace(" ", "_")}'

	indexReferenceGenome = gwf.target_from_template(
		name=f'index_reference_genome_{SPECIES_NAME.replace(" ", "_")}',
		template=index_reference_genome(
			referenceGenomeFile=REFERENCE_GENOME,
			outputDirectory=topDir
		)
	)

	fullBamList = []

	# Small, half-baked job batch system in case the workflow would produce to many jobs to queue at once.
	# Batched branch. Includes only generation of VCF parts.
	if NBATCHES != 0 and BATCHNR:
		batchSize = int(len(partitions)/NBATCHES)
		firstJob = (BATCHNR - 1) * batchSize
		lastJob = BATCHNR * batchSize
		if NBATCHES == BATCHNR:
			if len(partitions) - batchSize * NBATCHES != 0:
				lastJob += len(partitions) - batchSize * NBATCHES
		
		for GROUP in SAMPLE_LIST:
			if not GROUP['bamFileList']:
				continue
			GROUP_NAME: str = GROUP['groupName'].lower()

			for SAMPLE in GROUP['bamFileList']:
				SAMPLE_NAME: str = os.path.basename(os.path.dirname(SAMPLE))
				
				fullBamList.append(SAMPLE)

				# One VCF file per sample file.
				if FREEBAYES_MODE == 1 or FREEBAYES_MODE == 0:
					freebayesPartitionSingle = gwf.map(
						name=name_freebayes_partition_single,
						template_func=freebayes_partition_single,
						inputs=partitions[firstJob:lastJob],
						extra={'referenceGenomeFile': indexReferenceGenome.outputs['symlink'],
							   'bamFile': SAMPLE,
							   'outputDirectory': topDir,
							   'groupName': GROUP_NAME,
							   'sampleName': SAMPLE_NAME,
							   'ploidy': FREEBAYES_PLOIDY,
							   'bestNAlleles': FREEBAYES_BESTN,
							   'minAlternateFraction': FREEBAYES_MINALTFRC,
							   'minAlternateCount': FREEBAYES_MINALTCNT,
							   'memory': VCF_MEM,
							   'time': VCF_TIME}
					)
		
		# One VCF file containing all sample files.
		if FREEBAYES_MODE == 2 or FREEBAYES_MODE == 0:
			freebayesPartitionAll = gwf.map(
				name=name_freebayes_partition_all,
				template_func=freebayes_partition_all,
				inputs=partitions[firstJob:lastJob],
				extra={'referenceGenomeFile': indexReferenceGenome.outputs['symlink'],
					   'bamFiles': fullBamList,
					   'outputDirectory': topDir,
					   'speciesName': SPECIES_NAME,
					   'ploidy': FREEBAYES_PLOIDY,
					   'bestNAlleles': FREEBAYES_BESTN,
					   'minAlternateFraction': FREEBAYES_MINALTFRC,
					   'minAlternateCount': FREEBAYES_MINALTCNT,
					   'memory': VCF_MEM,
					   'time': VCF_TIME}
			)

	# Non-batched branch. Includes all jobs.
	else:
		vcfSingleDict = {}
		groupStatusDict = {}

		if REPEATS_BED:
			bedExcludeOverlapRepeats = gwf.target_from_template(
				name=f'intergenic_exluding_repeats_bed',
				template=bed_exclude_overlap(
					mainBedFile=INTERGENIC_BED,
					subtractionBedFile=REPEATS_BED,
					outputDirectory=topDir,
					speciesName=SPECIES_NAME
				)
			)

		else:
			extractSoftmaskedIntervals = gwf.target_from_template(
				name=f'extract_repetitive_intervals',
				template=extract_softmasked_intervals(
					referenceGenomeFile=REFERENCE_GENOME,
					outputDirectory=topOut if OUTPUT_DIR else topDir
				)
			)

			bedExcludeOverlapRepeats = gwf.target_from_template(
				name=f'intergenic_exluding_repeats_bed',
				template=bed_exclude_overlap(
					mainBedFile=INTERGENIC_BED,
					subtractionBedFile=extractSoftmaskedIntervals.outputs['bed'],
					outputDirectory=topDir,
					speciesName=SPECIES_NAME
				)
			)

		for GROUP in SAMPLE_LIST:
			if not GROUP['bamFileList']:
				continue
			GROUP_STATUS: str = GROUP['groupStatus'].lower() if GROUP['groupStatus'] else 'i'
			GROUP_NAME: str = GROUP['groupName'].lower()
			GROUP_BAMS: list = GROUP['bamFileList']
			vcfSingleList = []

			vcfSingleDict[GROUP_NAME] = []
			groupStatusDict[GROUP_NAME] = GROUP_STATUS

			depthDistribution = gwf.target_from_template(
				name=f'depth_distribution_{GROUP_NAME}',
				template=depth_distribution(
					bamFiles=GROUP_BAMS,
					outputDirectory=topDir,
					outputName=f'{species_abbreviation(SPECIES_NAME)}_{GROUP_NAME}.{'ingroup' if GROUP_STATUS == 'i' else 'outgroup'}'
				)
			)

			if GROUP_STATUS == 'i':
				depthDistributionPlot = gwf.target_from_template(
					name=f'depth_distribution_plot_{GROUP_NAME}',
					template=depth_distribution_plot(
						depthDistributionFile=depthDistribution.outputs['depth'],
						minCoverageThreshold=INGROUP_FILTERING_MINDP,
						outputDirectory=topOut if OUTPUT_DIR else topDir,
						mode=0
					)
				)

			elif GROUP_STATUS == 'o':
				depthDistributionPlot = gwf.target_from_template(
					name=f'depth_distribution_plot_{GROUP_NAME}',
					template=depth_distribution_plot(
						depthDistributionFile=depthDistribution.outputs['depth'],
						minCoverageThreshold=OUTGROUP_FILTERING_MINDP,
						outputDirectory=topOut if OUTPUT_DIR else topDir,
						mode=1
					)
				)

			sharedSitesWithinThresholdBed = gwf.target_from_template(
				name=f'depth_threshold_bed_{GROUP_NAME}',
				template=shared_sites_within_threshold_bed(
					depthDistributionFile=depthDistribution.outputs['depth'],
					depthDistributionTsv=depthDistributionPlot.outputs['tsv'],
					outputDirectory=topOut if OUTPUT_DIR else topDir,
					outputName=f'{species_abbreviation(SPECIES_NAME)}_{GROUP_NAME}.{'ingroup' if GROUP_STATUS == 'i' else 'outgroup'}'
				)
			)

			siteCountAll = gwf.target_from_template(
				name=f'site_count_all_{GROUP_NAME}',
				template=site_count_region(
					bamFiles=GROUP_BAMS,
					depthDistributionTsv=depthDistributionPlot.outputs['tsv'],
					bedFile=None,
					siteType='all',
					outputDirectory=topDir,
					outputName=f'{species_abbreviation(SPECIES_NAME)}_{GROUP_NAME}.{'ingroup' if GROUP_STATUS == 'i' else 'outgroup'}'
				)
			)

			siteCountIntergenic = gwf.target_from_template(
				name=f'site_count_intergenic_{GROUP_NAME}',
				template=site_count_region(
					bamFiles=GROUP_BAMS,
					depthDistributionTsv=depthDistributionPlot.outputs['tsv'],
					siteType='intergenic',
					outputDirectory=topDir,
					outputName=f'{species_abbreviation(SPECIES_NAME)}_{GROUP_NAME}.{'ingroup' if GROUP_STATUS == 'i' else 'outgroup'}',
					bedFile=INTERGENIC_BED
				)
			)

			siteCountRegionIntergenicExclRepeats = gwf.target_from_template(
				name=f'site_count_intergenic_excl_repeats_{GROUP_NAME}',
				template=site_count_region(
					bamFiles=GROUP_BAMS,
					depthDistributionTsv=depthDistributionPlot.outputs['tsv'],
					siteType='intergenic_excl_repeats',
					outputDirectory=topDir,
					outputName=f'{species_abbreviation(SPECIES_NAME)}_{GROUP_NAME}.{'ingroup' if GROUP_STATUS == 'i' else 'outgroup'}',
					bedFile=bedExcludeOverlapRepeats.outputs['bed']
				)
			)

			mergeSiteTablesSingle = gwf.target_from_template(
				name=f'merge_site_tables_single_{GROUP_NAME}',
				template=merge_site_tables(
					siteTables=[siteCountAll.outputs['sitetable'],
								siteCountIntergenic.outputs['sitetable'],
								siteCountRegionIntergenicExclRepeats.outputs['sitetable']],
					outputName=f'{species_abbreviation(SPECIES_NAME)}_{GROUP_NAME}.{'ingroup' if GROUP_STATUS == 'i' else 'outgroup'}.singlecall' if FREEBAYES_MODE == 0 else f'{species_abbreviation(SPECIES_NAME)}_{GROUP_NAME}.{'ingroup' if GROUP_STATUS == 'i' else 'outgroup'}',
					outputDirectory=topOut if OUTPUT_DIR else topDir
				)
			)

			for SAMPLE in GROUP['bamFileList']:
				SAMPLE_NAME: str = os.path.basename(os.path.dirname(SAMPLE))
				
				fullBamList.append(SAMPLE)

				# One VCF file per sample file.
				if FREEBAYES_MODE == 1 or FREEBAYES_MODE == 0:
					freebayesPartitionSingle = gwf.map(
						name=name_freebayes_partition_single,
						template_func=freebayes_partition_single,
						inputs=partitions,
						extra={'referenceGenomeFile': indexReferenceGenome.outputs['symlink'],
							   'bamFile': SAMPLE,
							   'outputDirectory': topDir,
							   'groupName': GROUP_NAME,
							   'sampleName': SAMPLE_NAME,
							   'ploidy': FREEBAYES_PLOIDY,
							   'bestNAlleles': FREEBAYES_BESTN,
							   'minAlternateFraction': FREEBAYES_MINALTFRC,
							   'minAlternateCount': FREEBAYES_MINALTCNT,
							   'memory': VCF_MEM,
							   'time': VCF_TIME}
					)

					if nSegments <= 1:
						concatenateFreebayesSingle = gwf.target_from_template(
							name=f'concatenate_freebayes_vcf_{GROUP_NAME}_{SAMPLE_NAME.replace("-", "_")}',
							template=concat_vcf(
								files=collect(freebayesPartitionSingle.outputs, ['vcf'])['vcfs'],
								outputName=f'{SAMPLE_NAME}.freebayes_n{FREEBAYES_BESTN}_p{FREEBAYES_PLOIDY}_minaltfrc{FREEBAYES_MINALTFRC}_minaltcnt{FREEBAYES_MINALTCNT}_singlecall',
								outputDirectory=f'{topOut}/{GROUP_NAME}/{SAMPLE_NAME}' if OUTPUT_DIR else f'{topDir}/raw_vcf/{GROUP_NAME}/{SAMPLE_NAME}',
								compress=True
							)
						)

						normVcfSingle = gwf.target_from_template(
							name=f'normalize_vcf_{GROUP_NAME}_{SAMPLE_NAME.replace("-", "_")}',
							template=norm_vcf(
								vcfFile=concatenateFreebayesSingle.outputs['concat_file'],
								referenceGenomeFile=REFERENCE_GENOME,
								outputName=f'{SAMPLE_NAME}.freebayes_n{FREEBAYES_BESTN}_p{FREEBAYES_PLOIDY}_minaltfrc{FREEBAYES_MINALTFRC}_minaltcnt{FREEBAYES_MINALTCNT}_singlecall',
								outputDirectory=f'{topDir}/raw_vcf/{GROUP_NAME}/{SAMPLE_NAME}'
							)
						)

					else:
						segmentList = []
						start = 0
						end = segmentSize
						collection = collect(freebayesPartitionSingle.outputs, ['vcf'])['vcfs']

						for i in range(nSegments):
							concatenateFreebayesSingleSegment = gwf.target_from_template(
								name=f'concatenate_freebayes_vcf_{GROUP_NAME}_{SAMPLE_NAME.replace("-", "_")}_segment_{i+1}',
								template=concat_vcf(
									files=collection[start : end],
									outputName=f'{SAMPLE_NAME}.freebayes_n{FREEBAYES_BESTN}_p{FREEBAYES_PLOIDY}_minaltfrc{FREEBAYES_MINALTFRC}_minaltcnt{FREEBAYES_MINALTCNT}.segment{i+1}',
									outputDirectory=f'{topDir}/raw_vcf/{GROUP_NAME}/{SAMPLE_NAME}/tmp',
									compress=True
								)
							)

							segmentList.append(concatenateFreebayesSingleSegment.outputs['concat_file'])
							if i < nSegments - 1:
								start = end
								end += segmentSize
							elif i == nSegments - 1:
								start = end
								end = nPartitions

						concatenateFreebayesSingle = gwf.target_from_template(
							name=f'concatenate_freebayes_vcf_{GROUP_NAME}_{SAMPLE_NAME.replace("-", "_")}_complete',
							template=concat_vcf(
								files=segmentList,
								outputName=f'{SAMPLE_NAME}.freebayes_n{FREEBAYES_BESTN}_p{FREEBAYES_PLOIDY}_minaltfrc{FREEBAYES_MINALTFRC}_minaltcnt{FREEBAYES_MINALTCNT}_singlecall',
								outputDirectory=f'{topOut}/{GROUP_NAME}/{SAMPLE_NAME}' if OUTPUT_DIR else f'{topDir}/raw_vcf/{GROUP_NAME}/{SAMPLE_NAME}',
								compress=True
							)
						)

						normVcfSingle = gwf.target_from_template(
							name=f'normalize_vcf_{GROUP_NAME}_{SAMPLE_NAME.replace("-", "_")}',
							template=norm_vcf(
								vcfFile=concatenateFreebayesSingle.outputs['concat_file'],
								referenceGenomeFile=REFERENCE_GENOME,
								outputName=f'{SAMPLE_NAME}.freebayes_n{FREEBAYES_BESTN}_p{FREEBAYES_PLOIDY}_minaltfrc{FREEBAYES_MINALTFRC}_minaltcnt{FREEBAYES_MINALTCNT}_singlecall',
								outputDirectory=f'{topDir}/raw_vcf/{GROUP_NAME}/{SAMPLE_NAME}'
							)
						)
					
					# Collect concatenated single population VCF files
					vcfSingleList.append(normVcfSingle.outputs['vcf'])
					# Collect concatenated single population VCF files by group
					vcfSingleDict[GROUP_NAME].append(normVcfSingle.outputs['vcf'])
				
			if len(vcfSingleList) > 1:
				mergeVcfSingle = gwf.target_from_template(
					name=f'merge_vcf_single_{GROUP_NAME}',
					template=merge_vcf(
						vcfFiles=vcfSingleList,
						outputName=f'{species_abbreviation(SPECIES_NAME)}_{GROUP_NAME}.freebayes_n{FREEBAYES_BESTN}_p{FREEBAYES_PLOIDY}_minaltfrc{FREEBAYES_MINALTFRC}_minaltcnt{FREEBAYES_MINALTCNT}_singlecall.norm',
						outputDirectory=f'{topDir}/raw_vcf'
					)
				)

			filterVcfSingle = gwf.target_from_template(
				name=f'filter_vcf_single_{GROUP_NAME}',
				template=filter_vcf(
					vcfFile=normVcfSingle.outputs['vcf'] if len(vcfSingleList) > 1 else vcfSingleList[0],
					depthDistributionTsv=depthDistributionPlot.outputs['tsv'],
					outputDirectory=f'{topOut}' if OUTPUT_DIR else f'{topDir}/filtered_vcf',
					groupStatus=GROUP_STATUS
				)
			)
		
		# One VCF file containing all sample files.
		if FREEBAYES_MODE == 2 or FREEBAYES_MODE == 0:
			freebayesPartitionAll = gwf.map(
				name=name_freebayes_partition_all,
				template_func=freebayes_partition_all,
				inputs=partitions,
				extra={'referenceGenomeFile': indexReferenceGenome.outputs['symlink'],
					   'bamFiles': fullBamList,
					   'outputDirectory': topDir,
					   'speciesName': SPECIES_NAME,
					   'ploidy': FREEBAYES_PLOIDY,
					   'bestNAlleles': FREEBAYES_BESTN,
					   'minAlternateFraction': FREEBAYES_MINALTFRC,
					   'minAlternateCount': FREEBAYES_MINALTCNT,
					   'memory': VCF_MEM,
					   'time': VCF_TIME}
			)

			if nSegments <= 1:
				concatenateFreebayesAll = gwf.target_from_template(
					name=f'concatenate_freebayes_vcf_all',
					template=concat_vcf(
						files=collect(freebayesPartitionAll.outputs, ['vcf'])['vcfs'],
						outputName=f'{species_abbreviation(SPECIES_NAME)}.freebayes_n{FREEBAYES_BESTN}_p{FREEBAYES_PLOIDY}_minaltfrc{FREEBAYES_MINALTFRC}_minaltcnt{FREEBAYES_MINALTCNT}',
						outputDirectory=f'{topOut}' if OUTPUT_DIR else f'{topDir}/raw_vcf',
						compress=True
					)
				)
			
			else:
				segmentList = []
				start = 0
				end = segmentSize
				collection = collect(freebayesPartitionAll.outputs, ['vcf'])['vcfs']

				for i in range(nSegments):
					concatenateFreebayesAllSegment = gwf.target_from_template(
						name=f'concatenate_freebayes_vcf_all_segment_{i+1}',
						template=concat_vcf(
							files=collection[start : end],
							outputName=f'{species_abbreviation(SPECIES_NAME)}.freebayes_n{FREEBAYES_BESTN}_p{FREEBAYES_PLOIDY}_minaltfrc{FREEBAYES_MINALTFRC}_minaltcnt{FREEBAYES_MINALTCNT}.segment{i+1}',
							outputDirectory=f'{topDir}/raw_vcf/tmp',
							compress=True
						)
					)

					segmentList.append(concatenateFreebayesAllSegment.outputs['concat_file'])
					if i < nSegments - 1:
						start = end
						end += segmentSize
					elif i == nSegments - 1:
						start = end
						end = nPartitions

				concatenateFreebayesAll = gwf.target_from_template(
					name=f'concatenate_freebayes_vcf_all_complete',
					template=concat_vcf(
						files=segmentList,
						outputName=f'{species_abbreviation(SPECIES_NAME)}.freebayes_n{FREEBAYES_BESTN}_p{FREEBAYES_PLOIDY}_minaltfrc{FREEBAYES_MINALTFRC}_minaltcnt{FREEBAYES_MINALTCNT}',
						outputDirectory=f'{topOut}' if OUTPUT_DIR else f'{topDir}/raw_vcf',
						compress=True
					)
				)

		# if FREEBAYES_MODE == 2 or FREEBAYES_MODE == 0:
		# 	normVcfAll = gwf.target_from_template(
		# 		name=f'normalize_vcf_all',
		# 		template=norm_vcf(
		# 			vcfFile=concatenateFreebayesAll.outputs['concat_file'],
		# 			referenceGenomeFile=indexReferenceGenome.outputs['symlink'],
		# 			outputName=f'{species_abbreviation(SPECIES_NAME)}.freebayes_n{FREEBAYES_BESTN}_p{FREEBAYES_PLOIDY}_minaltfrc{FREEBAYES_MINALTFRC}_minaltcnt{FREEBAYES_MINALTCNT}_allcall',
		# 			outputDirectory=f'{topOut}' if OUTPUT_DIR else f'{topDir}/raw_vcf'
		# 		)
		# 	)

		# 	filterVcfAll = gwf.target_from_template(
		# 		name=f'filter_vcf_all',
		# 		template=filter_vcf(
		# 			normVcfAll.outputs['vcf'],
		# 			depthDistributionTsv=depthDistributionPlot.outputs['tsv'],
		# 			outputDirectory=topOut if OUTPUT_DIR else f'{topDir}/filtered_vcf',
		# 			speciesName=SPECIES_NAME,
		# 			minDepth=INGROUP_FILTERING_MINDP
		# 		)
		# 	)

		# 	mergeSiteTablesAll = gwf.target_from_template(
		# 		name=f'merge_site_tables_all',
		# 		template=merge_site_tables(
		# 			siteTables=[siteCountAll.outputs['sitetable'],
		# 						siteCountIntergenic.outputs['sitetable'],
		# 						siteCountRegionIntergenicExclRepeats.outputs['sitetable'],
		# 						filterVcfAll.outputs['sitetable']],
		# 			outputName=f'{species_abbreviation(SPECIES_NAME)}.allcall' if FREEBAYES_MODE == 0 else f'{species_abbreviation(SPECIES_NAME)}',
		# 			outputDirectory=topOut if OUTPUT_DIR else topDir
		# 		)
		# 	)
	
	return gwf
