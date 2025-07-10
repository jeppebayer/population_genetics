#!/bin/env python3
from gwf import AnonymousTarget
from gwf.executors import Conda
import os, yaml

########################## Functions ##########################

def species_abbreviation(species_name: str) -> str:
	"""Creates species abbreviation from species name.

	:param str species_name:
		Species name written as *genus* *species*"""
	genus, species = species_name.replace(' ', '_').split('_')
	genus = genus[0].upper() + genus[1:3]
	species = species[0].upper() + species[1:3]
	return genus + species

########################## Templates ##########################

def mark_duplicates(alignmentFile: str, sampleName: str, outputDirectory: str, environment: str, opticalDistance: int = 0, group: str | None = None):
	"""
	Template: Mark duplicate alignments using :script:`samtools markdup`.
	
	Template I/O::
	
		inputs = {'alignment': alignment_file}
		outputs = {'markdup': *.markdup.bam, 'bai': *.markdup.bam.bai, 'stats': *.markdup.bam.markdupstats}
	
	:param str alignment_file:
		Path to input alignment file.
	:param str sample_name:
		Name of sample.
	:param str output_directory:
		Path to output directory.
	"""
	outputDirectory = f'{outputDirectory}'
	inputs = {'alignment': alignmentFile}
	outputs = {'markdup': f'{outputDirectory}/{sampleName}.markdup.bam',
			   'bai': f'{outputDirectory}/{sampleName}.markdup.bam.bai',
			   'stats': f'{outputDirectory}/stats/{sampleName}.markdupstats.txt'}
	options = {
		'cores': 18,
		'memory': '60g',
		'walltime': '08:00:00'
	}
	spec = f"""
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	echo "Conda Environment Info:"
	conda env export --from-history
	
	[ -d {outputDirectory}/tmp ] || mkdir -p {outputDirectory}/tmp
	[ -d {outputDirectory}/stats ] || mkdir -p {outputDirectory}/stats

	samtools collate \\
		--threads {options['cores'] - 1} \\
		-O \\
		-u \\
		{alignmentFile} \\
	| samtools fixmate \\
		--threads {options['cores'] - 1} \\
		-m \\
		-u \\
		- \\
		- \\
	| samtools sort \\
		--threads {options['cores'] - 1} \\
		-T {outputDirectory}/tmp \\
		-u \\
		- \\
	| samtools markdup \\
		--threads {options['cores'] - 1} \\
		-T {outputDirectory}/tmp \\
		--include-fails \\
		-S \\
		-d {opticalDistance} \\
		-f {outputDirectory}/{sampleName}.markdupstats.txt \\
		- \\
		{outputDirectory}/{sampleName}.markdup.bam
	
	samtools index \\
		--threads {options['cores'] - 1} \\
		-b \\
		{outputDirectory}/{sampleName}.markdup.bam \\
		{outputDirectory}/{sampleName}.markdup.bam.bai

	mv {outputDirectory}/{sampleName}.markdupstats.txt {outputDirectory}/stats

	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec, executor=Conda(environment), group=group)

def alignment_stats(alignmentFile: str, outputDirectory: str, environment: str, group: str | None = None):
	"""
	Template: template_description
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	outputDirectory = f'{outputDirectory}/stats'
	inputs = {'alignment': alignmentFile}
	outputs = {'stats': f'{outputDirectory}/{os.path.basename(os.path.splitext(alignmentFile)[0])}.stats.txt',
			   'flagstat': f'{outputDirectory}/{os.path.basename(os.path.splitext(alignmentFile)[0])}.flagstat.tsv'}
	options = {
		'cores': 20,
		'memory': '20g',
		'walltime': '04:00:00'
	}
	spec = f"""
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	echo "Conda Environment Info:"
	conda env export --from-history
	
	[ -d {outputDirectory} ] || mkdir -p {outputDirectory}
	
	samtools stats \\
		--threads {options['cores'] - 1} \\
		{alignmentFile} \\
		> {outputDirectory}/{os.path.basename(os.path.splitext(alignmentFile)[0])}.stats.txt
	
	samtools flagstat \\
		--threads {options['cores'] - 1} \\
		--output-fmt tsv \\
		{alignmentFile} \\
		> {outputDirectory}/{os.path.basename(os.path.splitext(alignmentFile)[0])}.flagstat.prog.tsv
	
	mv {outputDirectory}/{os.path.basename(os.path.splitext(alignmentFile)[0])}.flagstat.prog.tsv {outputs['flagstat']}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec, executor=Conda(environment), group=group)

def summary_stats(markdupstatsFile: str, statsFile: str, flagstatFile: str, sampleName: str, outputDirectory: str, environment: str, group: str | None = None):
	"""
	Template: template_description
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	outputDirectory = f'{outputDirectory}/stats'
	inputs = {'markdup': markdupstatsFile,
		   	  'stats': statsFile,
			  'flagstat': flagstatFile}
	outputs = {'stats': f'{outputDirectory}/{sampleName}.summaryStats.tsv'}
	options = {
		'cores': 1,
		'memory': '5g',
		'walltime': '00:30:00'
	}
	spec = f"""
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	echo "Conda Environment Info:"
	conda env export --from-history
	
	[ -d {outputDirectory} ] || mkdir -p {outputDirectory}
	
	awk \\
		'BEGIN{{
			FS = OFS = "\\t"
			print "sample", "reads", "mapped%", "mapped#", "properlyPaired%", "properlyPaired#", "duplicates%", "duplicates#", "opticalDuplicates%", "opticalDuplicates#"
		}}
		{{
			if (FNR == NR)
			{{
				first = $0
				next
			}}
			second = $0
		}}
		END{{
			print first, second
		}}' \\
		<(awk \\
			'BEGIN{{
				FS = OFS = "\\t"
			}}
			{{
				if ($3 == "total (QC-passed reads + QC-failed reads)")
				{{
					totalReads = $1
				}}
				if ($3 == "primary mapped")
				{{
					totalMapped = $1
				}}
				if ($3 == "primary mapped %")
				{{
					percentMapped = substr($1, 1, length($1) - 1)
				}}
				if ($3 == "properly paired")
				{{
					totalPaired = $1
				}}
				if ($3 == "properly paired %")
				{{
					percentPaired = substr($1, 1, length($1) - 1)
				}}
			}}
			END{{
				print "{sampleName}", totalReads, percentMapped, totalMapped, int(totalPaired / totalReads * 10000 + 0.5) / 100, totalPaired
			}}' \\
			{flagstatFile}) \\
		<(awk \\
			'BEGIN{{
				FS = ": "
				OFS = "\\t"
			}}
			{{
				if ($1 == "READ")
				{{
					totalReads = $2
				}}
				if ($1 == "DUPLICATE PAIR OPTICAL")
				{{
					totalOptical += $2
				}}
				if ($1 == "DUPLICATE SINGLE OPTICAL")
				{{
					totalOptical += $2
				}}
				if ($1 == "DUPLICATE NON PRIMARY")
				{{
					totalOptical += $2
				}}
				if ($1 == "DUPLICATE TOTAL")
				{{
					totalDuplicate = $2
				}}
			}}
			END{{
				print int(totalDuplicate / totalReads * 10000 + 0.5) / 100, totalDuplicate, int(totalOptical / totalDuplicate * 10000 + 0.5) / 100, totalOptical
			}}' \\
			{markdupstatsFile}) \\
		> {outputDirectory}/{sampleName}.summaryStats.prog.tsv

	mv {outputDirectory}/{sampleName}.summaryStats.prog.tsv {outputs['stats']}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec, executor=Conda(environment), group=group)

def collate_stats(summaryStatsFiles: list, outputDirectory: str, environment: str, group: str | None = None):
	"""
	Template: template_description
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	outputDirectory = f'{outputDirectory}'
	inputs = {'stats': summaryStatsFiles}
	outputs = {'stats': f'{outputDirectory}/summaryStats.tsv'}
	options = {
		'cores': 1,
		'memory': '5g',
		'walltime': '00:30:00'
	}
	spec = f"""
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	echo "Conda Environment Info:"
	conda env export --from-history
	
	[ -d {outputDirectory} ] || mkdir -p {outputDirectory}
	
	awk \\
		'BEGIN{{
			FS = OFS = "\\t"
			print "sample", "reads", "mapped%", "mapped#", "properlyPaired%", "properlyPaired#", "duplicates%", "duplicates#", "opticalDuplicates%", "opticalDuplicates#"
		}}
		{{
			if (FNR == 2)
			{{
				print $0
			}}
		}}' \\
		{' '.join(summaryStatsFiles)} \\
		> {outputDirectory}/summaryStats.prog.tsv
	
	mv {outputDirectory}/summaryStats.prog.tsv {outputs['stats']}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec, executor=Conda(environment), group=group)