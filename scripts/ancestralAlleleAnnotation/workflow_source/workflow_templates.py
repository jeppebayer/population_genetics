#!/bin/env python3
from gwf import AnonymousTarget
import os, glob

############################## Functions ##############################

def speciesAbbreviation(speciesName: str) -> str:
	"""Creates species abbreviation from species name.

	:param str speciesName:
		Species name written as *genus* *species*"""
	genus, species = speciesName.replace(' ', '_').split('_')
	genus = genus[0].upper() + genus[1:3]
	species = species[0].upper() + species[1:3]
	return genus + species

def parse_bed(bedFile: str):
	"""
	Parses :format:`BED` file returning a list of dictionaries
	containing chromsome name, number of sites within chromosome, and
	the start and end lines of the chromsome within the :format:`BED` file.

	::

		return [{'chrom': str, 'nSites': int, 'start': int, 'end': int}, ...]
	
	:param str bedFile:
		:format:`BED` file.
	"""
	bedList = []
	chromName = None
	siteCount = 0
	with open(bedFile, 'r') as bed:
		for lineNo, entry in enumerate(bed, start = 1):
			entry = entry.strip()
			entry = entry.split('\t')
			if chromName and chromName != entry[0]:
				end = lineNo - 1
				bedList.append({'chromName': chromName, 'nSites': siteCount, 'start': start, 'end': end})
				start = lineNo
				siteCount = 0
				chromName = entry[0]
			if not chromName:
				chromName = entry[0]
				start = lineNo
			siteCount += int(entry[2]) - int(entry[1])
		bedList.append({'chromName': chromName, 'nSites': siteCount, 'start': start, 'end': lineNo})
	return bedList

def partition_bed(parseBed: list, nLines: int = 250000):
	"""
	Takes the output from **parse_bed()** and splits each chromosome into chunks of nLines.

	::

		return [{'num': int, 'chromName': str, 'start' int, 'end': int}]

	:param list parseBed:
		Output from **parse_bed()** function. A list of dictionaries.
	:param int nLines:
		Number of lines to split each chromosome into.
	"""
	padding = 1
	for chrom in parseBed:
		wholeChunks = (chrom['end'] - chrom['start'] + 1) // nLines
		padding += (wholeChunks + 1)
	nPad = len(str(padding))
	bedPartition = []
	num = 1
	for chrom in parseBed:
		wholeChunks = (chrom['end'] - chrom['start'] + 1) // nLines
		partialChunk = (chrom['end'] - chrom['start'] + 1) - wholeChunks * nLines
		start = chrom['start']
		for chunk in range(wholeChunks):
			end = start + nLines - 1
			bedPartition.append({'num': f'{num:0{nPad}}', 'chromName': chrom['chromName'], 'start': start, 'end': end})
			start = end + 1
			num += 1
		if partialChunk:
			bedPartition.append({'num': f'{num:0{nPad}}', 'chromName': chrom['chromName'], 'start': start, 'end': start + partialChunk - 1})
			num += 1
	return bedPartition

############################## Templates ##############################

def index_reference_genome(referenceGenomeFile: str, outputDirectory: str):
	"""
	Template: Index reference genome with :script:`bwa index` and :script:`samtools faidx`.
	
	Template I/O::
	
		inputs = {'reference': referenceGenomeFile}
		outputs = {'symlink': *, 'bwa': [*.amb, *.ann, *.pac, *.bwt, *.sa], 'fai': *.fai}
	
	:param str referenceGenomeFile:
		Path to reference genome.
	:param str outputDirectory:
		Path to output directory.
	"""
	inputs = {'reference': referenceGenomeFile}
	outputs = {'symlink': f'{outputDirectory}/reference/{os.path.basename(referenceGenomeFile)}',
			   'bwa': [f'{outputDirectory}/reference/{os.path.basename(referenceGenomeFile)}.amb',
					   f'{outputDirectory}/reference/{os.path.basename(referenceGenomeFile)}.ann',
					   f'{outputDirectory}/reference/{os.path.basename(referenceGenomeFile)}.pac',
					   f'{outputDirectory}/reference/{os.path.basename(referenceGenomeFile)}.bwt',
					   f'{outputDirectory}/reference/{os.path.basename(referenceGenomeFile)}.sa'],
			   'fai': f'{outputDirectory}/reference/{os.path.basename(referenceGenomeFile)}.fai'}
	protect = [outputs['symlink'], outputs['bwa'][0], outputs['bwa'][1], outputs['bwa'][2], outputs['bwa'][3], outputs['bwa'][4], outputs['fai']]
	options = {
		'cores': 1,
		'memory': '10g',
		'walltime': '04:00:00'
	}
	spec = f"""
	# Sources environment
	if [ "$USER" == "jepe" ]; then
		source /home/"$USER"/.bashrc
		source activate popgen
	fi
	
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	
	[ -d {outputDirectory}/reference ] || mkdir -p {outputDirectory}/reference
	[ -e {outputDirectory}/reference/{os.path.basename(referenceGenomeFile)} ] && rm -f {outputDirectory}/reference/{os.path.basename(referenceGenomeFile)}
	ln -s {referenceGenomeFile} {outputDirectory}/reference/{os.path.basename(referenceGenomeFile)}
	
	bwa index \\
		-p {outputDirectory}/reference/{os.path.basename(referenceGenomeFile)} \\
		{outputDirectory}/reference/{os.path.basename(referenceGenomeFile)}
	
	samtools faidx \\
		-o {outputDirectory}/reference/{os.path.basename(referenceGenomeFile)}.prog.fai \\
		{outputDirectory}/reference/{os.path.basename(referenceGenomeFile)}
	
	mv {outputDirectory}/reference/{os.path.basename(referenceGenomeFile)}.prog.fai {outputs['fai']}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, protect=protect, options=options, spec=spec)

def name_freebayes_partition_outgroup(idx: str, target: AnonymousTarget) -> str:
	return f'freebayes_part_single_{os.path.basename(target.outputs["vcf"]).replace("-", "_").replace("|", "_")}'

def freebayes_partition_outgroup(referenceGenomeFile: str, bamFile: str, outputDirectory: str, sampleName: str, bedFile: str, num: int, chromName: str, start: int, end: int, ploidy: int = 100, bestNAlleles: int = 3, minAlternateFraction: float | int = 0, minAlternateCount: int = 2, memory: int = 80, time: str = '48:00:00'):
	"""
	Template: Create VCF file for each partition in a single pooled alignment.
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'reference': referenceGenomeFile,
		   	  'bam': bamFile}
	outputs = {'vcf': f'{outputDirectory}/outgroups/{sampleName}/tmp/{sampleName}.freebayes.{num}_{chromName.replace("|", "_")}.vcf.gz',
			   'index': f'{outputDirectory}/outgroups/{sampleName}/tmp/{sampleName}.freebayes.{num}_{chromName.replace("|", "_")}.vcf.gz.csi'}
	options = {
		'cores': 1,
		'memory': f'{memory}g',
		'walltime': time
	}
	spec = f"""
	# Sources environment
	if [ "$USER" == "jepe" ]; then
		source /home/"$USER"/.bashrc
		source activate popgen
	fi
	
	export _JAVA_OPTIONS="-Xmx{options['memory']}"

	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"

	[ -d {outputDirectory}/outgroups/{sampleName}/tmp ] || mkdir -p {outputDirectory}/outgroups/{sampleName}/tmp
	
	freebayes \\
		--fasta-reference {referenceGenomeFile} \\
		--use-best-n-alleles {bestNAlleles} \\
		--ploidy {ploidy} \\
		--targets <(awk \\
			'BEGIN{{
				FS = OFS = "\\t"
			}}
			{{
				if (NR >= {start} && NR <= {end})
				{{
					print $0
					next
				}}
				if (NR > {end})
				{{
					exit
				}}
			}}' \\
			{bedFile}) \\
		--min-alternate-fraction {minAlternateFraction} \\
		--min-alternate-count {minAlternateCount} \\
		--pooled-discrete \\
		--report-monomorphic \\
		-b {bamFile} \\
	| bcftools view \\
		--output-type z \\
		--output {outputDirectory}/outgroups/{sampleName}/tmp/{sampleName}.outgroups.freebayes.{num}_{chromName.replace("|", "_")}.prog.vcf.gz \\
		--write-index \\
		-
	
	mv {outputDirectory}/outgroups/{sampleName}/tmp/{sampleName}.freebayes.{num}_{chromName.replace("|", "_")}.prog.vcf.gz {outputs['vcf']}
	mv {outputDirectory}/outgroups/{sampleName}/tmp/{sampleName}.freebayes.{num}_{chromName.replace("|", "_")}.prog.vcf.gz.csi {outputs['index']}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def concat_vcf(files: list, outputName: str, outputDirectory: str| None = None, compress: bool = True):
	"""
	Template: Concatenates :format:`VCF` files. Optionally compresses output.
	
	Template I/O::
	
		inputs = {'files': list}
		outputs = {'concat_file': outputName.vcf | outputName.vcf.gz}
	
	:param list files:
		List containing :format:`VCF` files to concatenate. Can already be gzipped.
	:param str outputName:
		Desired name of output file, no extension.
	:param str outputDirectory:
		Path to output directory. Default is directory of 'files[0]'.
	:param bool compress:
		Bool indicating whether the output file should be compressed or not.
	"""
	if not outputDirectory:
		outputDirectory = os.path.dirname(files[0])
	inputs = {'files': files}
	if compress:
		outputs = {'concat_file': f'{outputDirectory}/{outputName}.vcf.gz',
			 	   'index': f'{outputDirectory}/{outputName}.vcf.gz.csi'}
	else:
		outputs = {'concat_file': f'{outputDirectory}/{outputName}.vcf'}
	options = {
		'cores': 32,
		'memory': '40g',
		'walltime': '24:00:00'
	}
	protect = outputs['concat_file']
	spec = f"""
	# Sources environment
	if [ "$USER" == "jepe" ]; then
		source /home/"$USER"/.bashrc
		source activate popgen
	fi
	
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	
	[ -d {outputDirectory} ] || mkdir -p {outputDirectory}
	
	if [ {compress} == 'False' ]; then
		bcftools concat \\
			--threads {options['cores']} \\
			--output-type v \\
			--output {outputDirectory}/{outputName}.prog.vcf \\
			{' '.join(files)}

			mv {outputDirectory}/{outputName}.prog.vcf {outputs['concat_file']}
	else
		bcftools concat \\
			--threads {options['cores']} \\
			--output-type z \\
			--output {outputDirectory}/{outputName}.prog.vcf.gz \\
			{' '.join(files)}
		
		bcftools index \\
			--threads {options['cores']} \\
			{outputDirectory}/{outputName}.prog.vcf.gz

			mv {outputDirectory}/{outputName}.prog.vcf.gz {outputs['concat_file']}
			mv {outputDirectory}/{outputName}.prog.vcf.gz.csi {outputs['index']}
	fi
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, protect=protect, spec=spec)

def merge_vcf(vcfFiles: list, outputName: str, outputDirectory: str):
	"""
	Template: template_description
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'vcfs': vcfFiles}
	outputs = {'vcf': f'{outputDirectory}/{outputName}.merge.vcf.gz',
			   'index': f'{outputDirectory}/{outputName}.merge.vcf.gz.csi'}
	protect = [outputs['vcf'], outputs['index']]
	options = {
		'cores': 30,
		'memory': '40g',
		'walltime': '24:00:00'
	}
	spec = f"""
	# Sources environment
	if [ "$USER" == "jepe" ]; then
		source /home/"$USER"/.bashrc
		source activate popgen
	fi
	
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	
	[ -d {outputDirectory} ] || mkdir -p {outputDirectory}
	
	bcftools merge \\
		--threads {options['cores']} \\
		--output-type z \\
		--output {outputDirectory}/{outputName}.merge.prog.vcf.gz \\
		--missing-to-ref \\
		--write-index \\
		{' '.join(vcfFiles)}
	
	mv {outputDirectory}/{outputName}.merge.prog.vcf.gz {outputs['vcf']}
	mv {outputDirectory}/{outputName}.merge.prog.vcf.gz.csi {outputs['index']}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, protect=protect, options=options, spec=spec)

def cp_rename_vcf(vcfFile: str, outputName: str, outputDirectory: str):
	"""
	Template: template_description
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'vcf': vcfFile}
	outputs = {'vcf': f'{outputDirectory}/{outputName}.vcf.gz',
			   'index': f'{outputDirectory}/{outputName}.vcf.gz.csi'}
	protect = [outputs['vcf'], outputs['index']]
	options = {
		'cores': 30,
		'memory': '20g',
		'walltime': '06:00:00'
	}
	spec = f"""
	# Sources environment
	if [ "$USER" == "jepe" ]; then
		source /home/"$USER"/.bashrc
		source activate popgen
	fi
	
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	
	[ -d {outputDirectory} ] || mkdir -p {outputDirectory}
	
	cp \\
		{vcfFile} \\
		{outputDirectory}/{outputName}.prog.vcf.gz

	bcftools index \\
		--threads {options['cores']} \\
		--csi \\
		{outputDirectory}/{outputName}.prog.vcf.gz
	
	mv {outputDirectory}/{outputName}.prog.vcf.gz {outputs['vcf']}
	mv {outputDirectory}/{outputName}.prog.vcf.gz.csi {outputs['index']}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, protect=protect, options=options, spec=spec)

def vcf_stats(vcfFile: str, referenceGenomeFile: str, outputDirectory: str):
	"""
	Template: template_description
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'vcf': vcfFile,
		   	  'reference': referenceGenomeFile}
	outputs = {'stats': f'{outputDirectory}/vcfStats/{os.path.basename(vcfFile)}.stats'}
	protect = [outputs['stats']]
	options = {
		'cores': 30,
		'memory': '20g',
		'walltime': '06:00:00'
	}
	spec = f"""
	# Sources environment
	if [ "$USER" == "jepe" ]; then
		source /home/"$USER"/.bashrc
		source activate popgen
	fi
	
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	
	[ -d {outputDirectory}/vcfStats ] || mkdir -p {outputDirectory}/vcfStats
	
	maxDepth="$(bcftools query -f '%INFO/DP' {vcfFile} | awk '{{if (max < $1) {{max = $1}}}} END{{print max}}' -)"

	bcftools stats \\
		--threads {options['cores']} \\
		--fasta-ref {referenceGenomeFile} \\
		--depth 0,"$maxDepth",1 \\
		--verbose \\
		{vcfFile} \\
		> {outputDirectory}/vcfStats/{os.path.basename(vcfFile)}.prog.stats
	
	mv {outputDirectory}/vcfStats/{os.path.basename(vcfFile)}.prog.stats {outputs['stats']}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, protect=protect, options=options, spec=spec)