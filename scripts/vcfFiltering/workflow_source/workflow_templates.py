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

def filter_vcf(vcfFile: str, depthThresholdBed: str, minDepth: int, maxDepth: int, outputDirectory: str):
	"""
	Template: template_description
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	filename = f'{os.path.splitext(os.path.splitext(os.path.basename(vcfFile))[0])[0] if vcfFile.endswith(".gz") else os.path.splitext(os.path.basename(vcfFile))[0]}.filter_AF0_SnpGap5_typeSnps_biallelic_AO1_DP{minDepth}-{maxDepth}'
	inputs = {'vcf': vcfFile,
		   	  'bed': depthThresholdBed}
	outputs = {'vcf': f'{outputDirectory}/{filename}.vcf.gz',
			   'index': f'{outputDirectory}/{filename}.vcf.gz.csi'}
	options = {
		'cores': 30,
		'memory': '30g',
		'walltime': '24:00:00'
	}
	protect = [outputs['vcf'], outputs['index']]
	spec = f"""
	# Sources environment
	if [ "$USER" == "jepe" ]; then
		source /home/"$USER"/.bashrc
		source activate popgen
	fi
	
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	
	[ -d {outputDirectory} ] || mkdir -p {outputDirectory}
	
	bcftools view \\
		--threads {options['cores']} \\
		--include 'INFO/AF > 0' \\
		--output-type u \\
		{vcfFile} \\
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
		--include 'INFO/AO > 1' \\
		--output-type z \\
		--output {outputDirectory}/{filename}.prog1.vcf.gz \\
		--write-index \\
		-
	
	bcftools view \\
		--threads {options['cores']} \\
		--regions-file {depthThresholdBed} \\
		--output-type z \\
		--output {outputDirectory}/{filename}.prog2.vcf.gz \\
		--write-index \\
		{outputDirectory}/{filename}.prog1.vcf.gz
	
	rm {outputDirectory}/{filename}.prog1.vcf.gz
	rm {outputDirectory}/{filename}.prog1.vcf.gz.csi
	mv {outputDirectory}/{filename}.prog2.vcf.gz {outputs['vcf']}
	mv {outputDirectory}/{filename}.prog2.vcf.gz.csi {outputs['index']}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, protect=protect, spec=spec)

def variable_site_count(vcfFileBefore: str, vcfFileAfter: str, outputDirectory: str):
	"""
	Template: template_description
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	filename = os.path.basename(os.path.splitext(os.path.splitext(vcfFileAfter)[0])[0]) if vcfFileAfter.endswith('.gz') else os.path.basename(os.path.splitext(vcfFileAfter)[0])
	inputs = {'before': vcfFileBefore,
			  'after': vcfFileAfter}
	outputs = {'sitetable': f'{outputDirectory}/sitetable/{filename}.variable.sitetable.tsv'}
	options = {
		'cores': 1,
		'memory': '50g',
		'walltime': '08:00:00'
	}
	spec = f"""
	# Sources environment
	if [ "$USER" == "jepe" ]; then
		source /home/"$USER"/.bashrc
		source activate popgen
	fi
	
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	
	[ -d {outputDirectory}/sitetable ] || mkdir -p {outputDirectory}/sitetable
	
	variableSiteCount() {{
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
					split($i, genotypeArray, ":")
					if (genotypeArray[1] !~ /^0\\/.*0$/)
					{{
						nPopulationVariants[populations[i], "whole_genome"] += 1
						nPopulationVariants[populations[i], $1] += 1
					}}
				}}
				nTotalVariants["whole_genome"] += 1
				nTotalVariants[$1] += 1
			}}
			END{{
				if (first == 1)
				{{
					print "stage", "identifier", "group", "region", "type", "site_count"
				}}
				for (i in nTotalVariants)
				{{
					print stage, identifier, "all", i, "variable", ntotalvariants[i]
				}}
				for (i in nPopulationVariants)
				{{
					split(i, popPegion, "\\034")
					print stage, identifier, popRegion[1], popRegion[2], "variable", nPopulationVariants[i]
				}}
			}}'
	}}
	
	bcftools view \\
		--threads {options['cores']} \\
		--output-type v \\
		{vcfFileBefore} \\
	| variableSiteCount \\
		1 \\
		0 \\
		"total" \\
		> {outputDirectory}/sitetable/{filename}.variable.sitetable.unsorted.tsv

	bcftools view \\
		--threads {options['cores']} \\
		--output-type v \\
		{vcfFileAfter} \\
	| variableSiteCount \\
		0 \\
		1 \\
		"filter" \\
		>> {outputDirectory}/sitetable/{filename}.variable.sitetable.unsorted.tsv

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
		{outputDirectory}/sitetable/{filename}.variable.sitetable.unsorted.tsv \\
		> {outputDirectory}/sitetable/{filename}.variable.sitetable.prog.tsv

	mv {outputDirectory}/sitetable/{filename}.variable.sitetable.prog.tsv {outputs['sitetable']}
	rm {outputDirectory}/sitetable/{filename}.variable.sitetable.unsorted.tsv
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

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

def snpeff_annotation(vcfFile: str, snpeffPredictorFile: str, outputDirectory: str, snpeffConfigFile: str = f'{os.path.dirname(os.path.realpath(__file__))}/software/snpEff.config'):
	"""
	Template: template_description
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	filename = os.path.splitext(os.path.splitext(os.path.basename(vcfFile))[0])[0] if vcfFile.endswith('.gz') else os.path.splitext(os.path.basename(vcfFile))[0]
	inputs = {'vcf': vcfFile,
		   	  'predictor': snpeffPredictorFile,
			  'config': snpeffConfigFile}
	outputs = {'vcf': f'{outputDirectory}/snpEff/{filename}.ann.vcf.gz',
			   'index': f'{outputDirectory}/snpEff/{filename}.ann.vcf.gz.csi',
			   'csv': f'{outputDirectory}/snpEff/{filename}.ann.vcf.gz.snpEffSummary.csv',
			   'html': f'{outputDirectory}/snpEff/{filename}.ann.vcf.gz.snpEffSummary.html'}
	protect = [outputs['csv'], outputs['html']]
	options = {
		'cores': 30,
		'memory': '80gg',
		'walltime': '16:00:00'
	}
	spec = f"""
	# Sources environment
	if [ "$USER" == "jepe" ]; then
		source /home/"$USER"/.bashrc
		source activate popgen
	fi
	
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	
	[ -d {outputDirectory}/snpEff ] || mkdir -p {outputDirectory}/snpEff
	
	export _JAVA_OPTIONS="-Xmx{options['memory']}"

	snpEff ann \\
		-csvStats {outputDirectory}/snpEff/{filename}.ann.vcf.gz.snpEffSummary.prog.csv \\
		-htmlStats {outputDirectory}/snpEff/{filename}.ann.vcf.gz.snpEffSummary.prog.html \\
		-nodownload \\
		-config {snpeffConfigFile} \\
		-dataDir {os.path.dirname(os.path.dirname(snpeffPredictorFile))} \\
		-verbose \\
		-i vcf \\
		-o vcf \\
		{os.path.basename(os.path.dirname(snpeffPredictorFile))} \\
		{vcfFile} \\
	| bcftools view \\
		--output-type z \\
		--output {outputDirectory}/snpEff/{filename}.ann.prog.vcf.gz \\
		--write-index \\
		-
	
	mv {outputDirectory}/snpEff/{filename}.ann.vcf.gz.snpEffSummary.prog.csv {outputs['csv']}
	mv {outputDirectory}/snpEff/{filename}.ann.vcf.gz.snpEffSummary.prog.html {outputs['html']}
	mv {outputDirectory}/snpEff/{filename}.ann.prog.vcf.gz {outputs['vcf']}
	mv {outputDirectory}/snpEff/{filename}.ann.prog.vcf.gz.csi {outputs['index']}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, protect=protect, options=options, spec=spec)

# def ancestral_allele_inferrence(outgroupVcf: str, depthThresholdBed: str, outputDirectory: str, speciesName: str):
# 	"""
# 	Template: template_description
	
# 	Template I/O::
	
# 		inputs = {}
# 		outputs = {}
	
# 	:param
# 	"""
# 	inputs = {}
# 	outputs = {}
# 	options = {
# 		'cores': 20,
# 		'memory': '30g',
# 		'walltime': '08:00:00'
# 	}
# 	spec = f"""
# 	# Sources environment
# 	if [ "$USER" == "jepe" ]; then
# 		source /home/"$USER"/.bashrc
# 		source activate popgen
# 	fi
	
# 	echo "START: $(date)"
# 	echo "JobID: $SLURM_JOBID"
	
# 	[ -d {outputDirectory} ] || mkdir -p {outputDirectory}
	
# 	bcftools view \\
# 		--threads {options['cores']} \\
# 		--include 'INFO/AF > 0' \\
# 		--output-type u \\
# 		{vcfFile} \\
# 	| bcftools filter \\
# 		--threads {options['cores']} \\
# 		--SnpGap 5:indel \\
# 		--output-type u \\
# 		- \\
# 	| bcftools view \\
# 		--threads {options['cores']} \\
# 		--types snps \\
# 		--output-type u \\
# 		- \\
# 	| bcftools view \\
# 		--threads {options['cores']} \\
# 		--max-alleles 2 \\
# 		--output-type u \\
# 		- \\
# 	| bcftools view \\
# 		--threads {options['cores']} \\
# 		--include 'INFO/AO > 1' \\
# 		--output-type z \\
# 		--output {outputDirectory}/{filename}.prog1.vcf.gz \\
# 		--write-index
	
# 	| bcftools view \\
# 		--threads {options['cores']} \\
# 		--regions-file {depthThresholdBed} \\
# 		--output-type z \\
# 		--output {outputDirectory}/{filename}.prog1.vcf.gz \\
# 		--write-index \\
# 		{outputDirectory}/{filename}.prog2.vcf.gz
	
# 	bcftools query --format '%CHROM\t%POS0\t%END\n'

# 	mv
	
# 	echo "END: $(date)"
# 	echo "$(jobinfo "$SLURM_JOBID")"
# 	"""
# 	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)