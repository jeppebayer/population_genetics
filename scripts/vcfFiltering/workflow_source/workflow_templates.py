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

########################## Templates ##########################

def filter_vcf(vcfFile: str, depthThresholdBed: str, minQuality: float, outputDirectory: str, groupStatus: str):
	"""
	Template: template_description
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'vcf': vcfFile,
		   	  'bed': depthThresholdBed}
	outputs = {'vcf': f'{outputDirectory}/{os.path.splitext(os.path.splitext(os.path.basename(vcfFile))[0])[0] if vcfFile.endswith(".gz") else os.path.splitext(os.path.basename(vcfFile))[0]}.bcftoolsfilter_AF0_SnpGap5_typesnps_biallelic_DPdynamic_AO1.{'ingroup' if groupStatus == 'i' else 'outgroup'}.vcf.gz',
			   'index': f'{outputDirectory}/{os.path.splitext(os.path.splitext(os.path.basename(vcfFile))[0])[0] if vcfFile.endswith(".gz") else os.path.splitext(os.path.basename(vcfFile))[0]}.bcftoolsfilter_AF0_SnpGap5_typesnps_biallelic_DPdynamic_AO1.{'ingroup' if groupStatus == 'i' else 'outgroup'}.vcf.gz.csi',
			   'sitetable': f'{outputDirectory}/sitetable/{os.path.basename(vcfFile).split('.')[0]}.{'ingroup' if groupStatus == 'i' else 'outgroup'}.sitetable.variable.tsv'}
	options = {
		'cores': 18,
		'memory': '30g',
		'walltime': '24:00:00'
	}
	protect = [outputs['vcf'], outputs['index'], outputs['sitetable']]
	spec = f"""
	# Sources environment
	if [ "$USER" == "jepe" ]; then
		source /home/"$USER"/.bashrc
		source activate popgen
	fi
	
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	
	[ -d {outputDirectory}/sitetable ] || mkdir -p {outputDirectory}/sitetable

	variablesitecount() {{
		awk \\
			-v first="$1" \\
			-v stage="$2" \\
			-v identifier="$3" \\
			'BEGIN{{
				FS = OFS = "\\t"
			}}
			{{
				if ($0 ~ /^##/)
				{{
					next
				}}
				if ($0 ~ /^#CHROM/)
				{{
					for (i = 10; i <= NF; i++)
					{{
						populations[i] = $i
					}}
					next
				}}
				for (i = 10; i <= NF; i++)
				{{
					split($i, genotypearray, ":")
					if (genotypearray[1] !~ /^0\\/.*0$/)
					{{
						npopulationvariants[populations[i], "whole_genome"] += 1
						npopulationvariants[populations[i], $1] += 1
					}}
				}}
				ntotalvariants["whole_genome"] += 1
				ntotalvariants[$1] += 1
			}}
			END{{
				if (first == 1)
				{{
					print "stage", "identifier", "group", "region", "type", "site_count"
				}}
				for (i in ntotalvariants)
				{{
					print stage, identifier, "all", i, "variable", ntotalvariants[i]
				}}
				for (i in npopulationvariants)
				{{
					split(i, popregion, "\\034")
					print stage, identifier, popregion[1], popregion[2], "variable", npopulationvariants[i]
				}}
			}}'
	}}
	
	bcftools view \\
		--threads {options['cores']} \\
		--output-type v \\
		{vcfFile} \\
	| tee \\
		>(variablesitecount \\
			1 \\
			0 \\
			"total" \\
			> {outputDirectory}/sitetable/{os.path.basename(vcfFile).split('.')[0]}.{'ingroup' if groupStatus == 'i' else 'outgroup'}.sitetable.variable.unsorted.tsv) \\
	| bcftools view \\
		--threads {options['cores']} \\
		--include 'INFO/AF > 0' \\
		--output-type u \\
		- \\
	| bcftools filter \\
		--threads {options['cores']} \\
		--SnpGap 5:indel \\
		--output-type u \\
		- \\
	| bcftools view \\
		--threads {options['cores']} \\
		--types snps \\
		--output-type u \\
		- \\
	| bcftools view \\
		--threads {options['cores']} \\
		--max-alleles 2 \\
		--output-type u \\
		- \\
	| bcftools view \\
		--threads {options['cores']} \\
		--include "FMT/DP>=$mindepth & FMT/DP<=$maxdepth" \\
		--output-type u \\
		- \\
	| bcftools view \\
		--threads {options['cores']} \\
		--include 'FMT/AO > 1' \\
		--output-type v \\
		- \\
	| tee \\
		>(variablesitecount \\
			0 \\
			1 \\
			"AO>1" \\
			>> {outputDirectory}/sitetable/{os.path.basename(vcfFile).split('.')[0]}.{'ingroup' if groupStatus == 'i' else 'outgroup'}.sitetable.variable.unsorted.tsv) \\
	| bcftools view \\
		--threads {options['cores']} \\
		--output-type z \\
		--output {outputDirectory}/{os.path.splitext(os.path.splitext(os.path.basename(vcfFile))[0])[0] if vcfFile.endswith(".gz") else os.path.splitext(os.path.basename(vcfFile))[0]}.bcftoolsfilter_AF0_SnpGap5_typesnps_biallelic_DPdynamic_AO1.{'ingroup' if groupStatus == 'i' else 'outgroup'}.prog.vcf.gz \\
		--write-index \\
		-
	
	awk \\
		'BEGIN{{
			FS = OFS = "\\t"
		}}
		{{
			if (NR == 1)
			{{
				print $0
				next
			}}
			print $0 | "sort -k 1,1 -k 3,3 -k 4,4"
		}}' \\
		{outputDirectory}/sitetable/{os.path.basename(vcfFile).split('.')[0]}.{'ingroup' if groupStatus == 'i' else 'outgroup'}.sitetable.variable.unsorted.tsv \\
		> {outputDirectory}/sitetable/{os.path.basename(vcfFile).split('.')[0]}.{'ingroup' if groupStatus == 'i' else 'outgroup'}.sitetable.variable.prog.tsv
		
	mv {outputDirectory}/{os.path.splitext(os.path.splitext(os.path.basename(vcfFile))[0])[0] if vcfFile.endswith(".gz") else os.path.splitext(os.path.basename(vcfFile))[0]}.bcftoolsfilter_AF0_SnpGap5_typesnps_biallelic_DPdynamic_AO1.{'ingroup' if groupStatus == 'i' else 'outgroup'}.prog.vcf.gz {outputs['vcf']}
	mv {outputDirectory}/{os.path.splitext(os.path.splitext(os.path.basename(vcfFile))[0])[0] if vcfFile.endswith(".gz") else os.path.splitext(os.path.basename(vcfFile))[0]}.bcftoolsfilter_AF0_SnpGap5_typesnps_biallelic_DPdynamic_AO1.{'ingroup' if groupStatus == 'i' else 'outgroup'}.prog.vcf.gz.csi {outputs['index']}
	mv {outputDirectory}/sitetable/{os.path.basename(vcfFile).split('.')[0]}.{'ingroup' if groupStatus == 'i' else 'outgroup'}.sitetable.variable.prog.tsv {outputs['sitetable']}
	rm {outputDirectory}/sitetable/{os.path.basename(vcfFile).split('.')[0]}.{'ingroup' if groupStatus == 'i' else 'outgroup'}.sitetable.variable.unsorted.tsv
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, protect=protect, spec=spec)