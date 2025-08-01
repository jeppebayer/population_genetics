#!/bin/env python3
from gwf import AnonymousTarget
from gwf.executors import Conda
import os, yaml

########################## Functions ##########################

def speciesAbbreviation(speciesName: str) -> str:
	"""Creates species abbreviation from species name.

	:param str speciesName:
		Species name written as *genus* *species*"""
	genus, species = speciesName.replace(' ', '_').split('_')
	genus = genus[0].upper() + genus[1:3]
	species = species[0].upper() + species[1:3]
	return genus + species

############################## Templates ##############################

def sfs_count(vcfFile: str, intergenicBed: str, repeatsBed: str, outputDirectory: str, outputName: str, environment: str, group: str | None = None):
	"""
	Template: template_description
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	outputDirectory = f'{outputDirectory}'
	inputs = {'vcf': vcfFile,
		   	  'intergenic': intergenicBed,
			  'repeats': repeatsBed}
	outputs = {'sfs': f'{outputDirectory}/{outputName}.sfs'}
	options = {
		'cores': 1,
		'memory': '10g',
		'walltime': '06:00:00'
	}
	spec = f"""
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	
	[ -d {outputDirectory} ] || mkdir -p {outputDirectory}
	
	bcftools query \\
		-f '[%CHROM\\t%SAMPLE\\t%GT\\n]' \\
		--targets-file <(bedtools subtract \\
			-a {intergenicBed} \\
			-b {repeatsBed}) \\
		{vcfFile} \\
	| awk \\
		'BEGIN{{
			FS = OFS = "\\t"
		}}
		{{
			split($3, gtArray, "/")
			gtValue = 0
			for (i in gtArray)
			{{
				gtValue += gtArray[i]
			}}
			if (gtValue > 50)
			{{
				gtValue = 100 - gtValue
			}}
			sfsArray[$2, gtValue] += 1
		}}
		END{{
			for (i in sfsArray)
			{{
				split(i, sampleChromBin, "\\034")
				print sampleChromBin[1], sampleChromBin[2], sfsArray[i]
			}}
		}}' \\
		- \\
	| sort \\
	 	-k1,1 \\
		-k2,2n \\
		- \\
		> {outputDirectory}/{outputName}.prog.sfs
	
	mv {outputDirectory}/{outputName}.prog.sfs {outputs['sfs']}
	
	cat <<-VERSIONS
	{sfs_count.__name__}:
		bcftools: $(bcftools --version | head -n 1 | sed 's/bcftools //')
		bedtools: $(bedtools --version | sed 's/bedtools v//')
		GNU awk: $(awk --version | head -n 1 | sed 's/,.*//' | sed 's/GNU Awk //')
		sort: $(sort --version | head -n 1 | sed 's/sort (GNU coreutils) //')
	VERSIONS
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec, executor=Conda(environment), group=group)

def sfs_plot(sfsFile: str, outputDirectory: str, environment: str, group: str | None = None, sfsPlot: str = f'{os.path.dirname(os.path.realpath(__file__))}/software/sfsPlot.py'):
	"""
	Template: template_description
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	outputDirectory = f'{outputDirectory}'
	inputs = {'sfs': sfsFile}
	outputs = {'pdf': f'{outputDirectory}/{os.path.basename(sfsFile)}.pdf'}
	options = {
		'cores': 2,
		'memory': '30g',
		'walltime': '01:00:00'
	}
	spec = f"""
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	
	[ -d {outputDirectory} ] || mkdir -p {outputDirectory}
	
	python {sfsPlot} \\
		{sfsFile} \\
		{outputDirectory}
	
	cat <<-VERSIONS
	{sfs_plot.__name__}:
		python: $(python --version | sed 's/Python //')
	VERSIONS
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec, executor=Conda(environment), group=group)

def get_versions(versionFiles: list, environment: str, group: str | None = None):
	"""
	Template: Collect version information from all jobs
	
	Template I/O::
	
		inputs = {'versions': list}
		outputs = {'versions': str}
	
	:param list versionFiles:
		List of version files from all jobs.
	"""
	inputs = {'versions': versionFiles}
	outputs = {'versions': f'{os.getcwd()}/versions.yml'}
	options = {
		'cores': 1,
		'memory': '2g',
		'walltime': '00:10:00'
	}
	spec = f"""
	cat <<-VERSIONS > {outputs['versions']}
	workflow:
		gwf: $(gwf --version | sed 's/gwf, version //')
	VERSIONS
	
	cat \\
		{' '.join(versionFiles)} \\
		>> {os.getcwd()}/versions.yml
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec, executor=Conda(environment), group=group)