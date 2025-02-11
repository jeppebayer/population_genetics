#!/bin/env python3
from gwf import AnonymousTarget
import os, glob

########################## Functions ##########################

def speciesAbbreviation(speciesName: str) -> str:
	"""Creates species abbreviation from species name.

	:param str speciesName:
		Species name written as *genus* *species*"""
	genus, species = speciesName.replace(' ', '_').split('_')
	genus = genus[0].upper() + genus[1:3]
	species = species[0].upper() + species[1:3]
	return genus + species

########################## Targets ############################

def name_depth_fractions(idx: str, target: AnonymousTarget) -> str:
	return f'depth_fractions_{os.path.basename(target.outputs['tsv']).replace("-", "_")}'

def depth_fractions(bamFile: str, lowerThreshold: int, outputDirectory: str, sampleName: str, depthDistributionFractions: str = f'{os.path.dirname(os.path.realpath(__file__))}/software/depthDistributionFractions.py'):
	"""
	Template: Determine the coverage needed to capture various fractions of the data.
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'bam': bamFile}
	outputs = {'tsv': f'{outputDirectory}/{sampleName}.tsv',
			   'png': f'{outputDirectory}/{sampleName}.png'}
	protect = [outputs['png']]
	options = {
		'cores': 10,
		'memory': '50g',
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
	
	samtools depth \\
		--threads {options['cores'] - 1} \\
		{bamFile} \\
	| python {depthDistributionFractions} \\
		- \\
		{sampleName} \\
		{lowerThreshold} \\
		{outputDirectory}/{sampleName}.prog
	
	mv {outputDirectory}/{sampleName}.prog.tsv {outputs['tsv']}
	mv {outputDirectory}/{sampleName}.prog.png {outputs['png']}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, protect=protect, options=options, spec=spec)

def concat_depth_fractions_tsv(tsvFiles: list, outputDirectory: str, outputName: str):
	"""
	Template: Concatenate :format:`tsv` files generate from depthFractions.
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'tsvFiles': tsvFiles}
	outputs = {'concatTsv': f'{outputDirectory}/{outputName}.tsv'}
	protect = [outputs['concatTsv']]
	options = {
		'cores': 1,
		'memory': '10g',
		'walltime': '01:00:00'
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
	
	awk \\
		'BEGIN{{
			FS = OFS = "\\t"
		}}
		{{
			if (NR == 1 && FNR == 1)
			{{
				print $0
				next
			}}
			if (FNR >= 2)
			{{
				print $0
			}}
		}}' \
		{' '.join(tsvFiles)} \
		> {outputDirectory}/{outputName}.prog.tsv
	
	mv {outputDirectory}/{outputName}.prog.tsv {outputs['concatTsv']}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, protect=protect, options=options, spec=spec)

def concat_depth_fractions_plots(plotFiles: list, outputDirectory: str, outputName: str, nColumns: int = 4, concatenateImages: str = f'{os.path.dirname(os.path.realpath(__file__))}/software/concatenateImages.py'):
	"""
	Template: template_description
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'plotFiles': plotFiles}
	outputs = {'concatPlots': f'{outputDirectory}/{outputName}.png'}
	protect = [outputs['concatPlots']]
	options = {
		'cores': 10,
		'memory': '50g',
		'walltime': '01:00:00'
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
	
	python {concatenateImages} \\
		{outputDirectory} \\
		{outputName}.prog.png \\
		{nColumns} \\
		{' '.join(plotFiles)}
	
	mv {outputDirectory}/{outputName}.prog.png {outputs['concatPlots']}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, protect=protect, options=options, spec=spec)