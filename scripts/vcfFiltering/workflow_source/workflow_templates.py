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

def filter_vcf(vcfFile: str, depthThresholdBed: str, minDepth: int, maxDepth: int, outputDirectory: str):
	"""
	Template: template_description
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	filename = f'{os.path.splitext(os.path.splitext(os.path.basename(vcfFile))[0])[0] if vcfFile.endswith(".gz") else os.path.splitext(os.path.basename(vcfFile))[0]}.bcftoolsfilter_AF0_SnpGap5_typeSnps_biallelic_DP{minDepth}-{maxDepth}_AO1'
	inputs = {'vcf': vcfFile,
		   	  'bed': depthThresholdBed}
	outputs = {'vcf': f'{outputDirectory}/{filename}.vcf.gz',
			   'index': f'{outputDirectory}/{filename}.vcf.gz.csi',
			   'sitetable': f'{outputDirectory}/sitetable/{filename}.sitetable.variable.tsv'}
	options = {
		'cores': 30,
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
			> {outputDirectory}/sitetable/{filename}.sitetable.variable.unsorted.tsv) \\
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
		--regions-file {depthThresholdBed} \\
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
			"filterPass0" \\
			>> {outputDirectory}/sitetable/{filename}.sitetable.variable.unsorted.tsv) \\
	| bcftools view \\
		--threads {options['cores']} \\
		--output-type z \\
		--output {outputDirectory}/{filename}.prog.vcf.gz \\
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
		{outputDirectory}/sitetable/{filename}.sitetable.variable.unsorted.tsv \\
		> {outputDirectory}/sitetable/{filename}.sitetable.variable.prog.tsv
		
	mv {outputDirectory}/{filename}.prog.vcf.gz {outputs['vcf']}
	mv {outputDirectory}/{filename}.prog.vcf.gz.csi {outputs['index']}
	mv {outputDirectory}/sitetable/{filename}.sitetable.variable.prog.tsv {outputs['sitetable']}
	rm {outputDirectory}/sitetable/{filename}.sitetable.variable.unsorted.tsv
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, protect=protect, spec=spec)

def snpeff_build_database(referenceGenome: str, gtfAnnotation: str, outputDirectory: str, speciesName: str, snpeffConfig: str = f'{os.path.dirname(os.path.realpath(__file__))}/software/snpEff.config'):
	"""
	Template: template_description
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'reference': referenceGenome,
		   	  'gtf': gtfAnnotation}
	outputs = {'cds': f'{outputDirectory}/snpeff/data/{os.path.basename(os.path.splitext(referenceGenome)[0])}/cds.fa',
			   'protein': f'{outputDirectory}/snpeff/data/{os.path.basename(os.path.splitext(referenceGenome)[0])}/protein.fa',
			   'sequences': f'{outputDirectory}/snpeff/data/{os.path.basename(os.path.splitext(referenceGenome)[0])}/sequences.fa',
			   'genes': f'{outputDirectory}/snpeff/data/{os.path.basename(os.path.splitext(referenceGenome)[0])}/genes.gtf',
			   'predictor': f'{outputDirectory}/snpeff/data/{os.path.basename(os.path.splitext(referenceGenome)[0])}/snpEffectPredictor.bin'}
	options = {
		'cores': 10,
		'memory': '80g',
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
	
	[ -d {outputDirectory}/snpeff/data/{os.path.basename(os.path.splitext(referenceGenome)[0])} ] || mkdir -p {outputDirectory}/snpeff/data/{os.path.basename(os.path.splitext(referenceGenome)[0])}
	
	echo "# {speciesName} genome, {os.path.basename(referenceGenome)}" >> {snpeffConfig}
	echo "{os.path.basename(os.path.splitext(referenceGenome)[0])}.genome : {speciesName}" >> {snpeffConfig}
	echo "{os.path.basename(os.path.splitext(referenceGenome)[0])}.file_location : {referenceGenome}" >> {snpeffConfig}
	echo "{os.path.basename(os.path.splitext(referenceGenome)[0])}.addition_date : $(date +%d'/'%m'/'%Y)" >> {snpeffConfig}
	echo -e "{os.path.basename(os.path.splitext(referenceGenome)[0])}.data_directory : {outputDirectory}/snpeff/data/{os.path.basename(os.path.splitext(referenceGenome)[0])}\\n" >> {snpeffConfig}

	cp {referenceGenome} {outputDirectory}/snpeff/data/{os.path.basename(os.path.splitext(referenceGenome)[0])}/sequences.fa
	cp {gtfAnnotation} {outputDirectory}/snpeff/data/{os.path.basename(os.path.splitext(referenceGenome)[0])}/genes.gtf

	cd {outputDirectory}/snpeff/data/{os.path.basename(os.path.splitext(referenceGenome)[0])}

	agat_sp_extract_sequences.pl \\
		--gff {gtfAnnotation} \\
		--fasta {referenceGenome} \\
		--type cds \\
		--output {outputDirectory}/snpeff/data/{os.path.basename(os.path.splitext(referenceGenome)[0])}/cds.fa
	
	agat_sp_extract_sequences.pl \\
		--gff {gtfAnnotation} \\
		--fasta {referenceGenome} \\
		--type cds \\
		--protein \\
		--output {outputDirectory}/snpeff/data/{os.path.basename(os.path.splitext(referenceGenome)[0])}/protein.fa

	export _JAVA_OPTIONS="-Xmx{options['memory']}"

	snpEff build \\
		-gtf22 \\
		-config {snpeffConfig} \\
		-dataDir {outputDirectory}/snpeff/data \\
		-nodownload \\
		-verbose \\
		{os.path.basename(os.path.splitext(referenceGenome)[0])}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)