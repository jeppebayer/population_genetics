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
		--output-type u \\
		- \\
	| bcftools view \\
		--threads {options['cores']} \\
		--targets-file {depthThresholdBed} \\
		--output-type z \\
		--output {outputDirectory}/{filename}.prog.vcf.gz \\
		--write-index
	
	mv {outputDirectory}/{filename}.prog.vcf.gz {outputs['vcf']}
	mv {outputDirectory}/{filename}.prog.vcf.gz.csi {outputs['index']}
	
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

def snpeff_build_database(referenceGenome: str, gtfAnnotation: str, outputDirectory: str, speciesName: str, snpeffConfig: str = f'{os.path.dirname(os.path.realpath(__file__))}/software/snpeff/snpEff.config'):
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
	protect = [file for file in outputs]
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
	return AnonymousTarget(inputs=inputs, outputs=outputs, protect=protect, options=options, spec=spec)

def snpeff_annotation(vcfFile: str, snpeffPredictorFile: str, vcfOutputDirectory: str, snpeffOutputDirectory: str, snpeffConfigFile: str = f'{os.path.dirname(os.path.realpath(__file__))}/software/snpeff/snpEff.config'):
	"""
	Template: template_description
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	filename = os.path.splitext(os.path.splitext(os.path.basename(vcfFile))[0])[0] if vcfFile.endswith('.gz') else os.path.splitext(os.path.basename(vcfFile))[0]
	inputs = {'vcf': vcfFile,
		   	  'predictor': snpeffPredictorFile}
	outputs = {'vcf': f'{vcfOutputDirectory}/{filename}.ann.vcf.gz',
			   'index': f'{vcfOutputDirectory}/{filename}.ann.vcf.gz.csi',
			   'csv': f'{snpeffOutputDirectory}/snpEff/{filename}.ann.vcf.gz.snpEffSummary.csv',
			   'html': f'{snpeffOutputDirectory}/snpEff/{filename}.ann.vcf.gz.snpEffSummary.html',
			   'genes': f'{snpeffOutputDirectory}/snpEff/{filename}.ann.vcf.gz.snpEffSummary.genes.txt'}
	protect = [outputs['csv'], outputs['html']]
	options = {
		'cores': 30,
		'memory': '80g',
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
	
	[ -d {snpeffOutputDirectory}/snpEff ] || mkdir -p {snpeffOutputDirectory}/snpEff
	[ -d {vcfOutputDirectory} ] || mkdri -p {vcfOutputDirectory}
	
	export _JAVA_OPTIONS="-Xmx{options['memory']}"

	snpEff ann \\
		-csvStats {snpeffOutputDirectory}/snpEff/{filename}.ann.vcf.gz.snpEffSummary.prog.csv \\
		-htmlStats {snpeffOutputDirectory}/snpEff/{filename}.ann.vcf.gz.snpEffSummary.prog.html \\
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
		--output {vcfOutputDirectory}/{filename}.ann.prog.vcf.gz \\
		--write-index \\
		-
	
	mv {snpeffOutputDirectory}/snpEff/{filename}.ann.vcf.gz.snpEffSummary.prog.csv {outputs['csv']}
	mv {snpeffOutputDirectory}/snpEff/{filename}.ann.vcf.gz.snpEffSummary.prog.html {outputs['html']}
	mv {snpeffOutputDirectory}/snpEff/{filename}.ann.vcf.gz.snpEffSummary.prog.genes.txt {outputs['genes']}
	mv {vcfOutputDirectory}/{filename}.ann.prog.vcf.gz {outputs['vcf']}
	mv {vcfOutputDirectory}/{filename}.ann.prog.vcf.gz.csi {outputs['index']}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, protect=protect, options=options, spec=spec)

def ancestral_allele_inference_reference(referenceGenome: str, outgroupSitesWithinCoverageBed: str, outputDirectory: str, outputName: str):
	"""
	Template: template_description
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'reference': referenceGenome,
		   	  'bed': outgroupSitesWithinCoverageBed}
	outputs = {'annotation': f'{outputDirectory}/ancestral_allele/{outputName}.annotation.reference.tsv.gz'}
	options = {
		'cores': 20,
		'memory': '30g',
		'walltime': '20:00:00'
	}
	spec = f"""
	# Sources environment
	if [ "$USER" == "jepe" ]; then
		source /home/"$USER"/.bashrc
		source activate popgen
	fi
	
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	
	[ -d {outputDirectory}/ancestral_allele ] || mkdir -p {outputDirectory}/ancestral_allele

	bedtools getfasta \\
		-fi {referenceGenome} \\
		-bed <(awk \\
			'BEGIN{{
				FS = OFS = "\\t"
			}}
			{{
				for (pos = $2; pos < $3; pos++)
				{{
					print $1, pos, pos + 1
				}}
			}}' \\
			{outgroupSitesWithinCoverageBed}) \\
		-bedOut \\
	| awk \\
		'BEGIN{{
			FS = OFS ="\\t"
		}}
		{{
			print $1, $3, $4
		}}' \\
		- \\
	| bgzip \\
		--threads {options['cores']} \\
		--output {outputDirectory}/ancestral_allele/{outputName}.annotation.reference.prog.tsv.gz \\
		-

	mv {outputDirectory}/ancestral_allele/{outputName}.annotation.reference.prog.tsv.gz {outputs['annotation']}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def ancestral_allele_inference_variant(outgroupVcf: str, outgroupSitesWithinCoverageBed: str, outputDirectory: str, outputName: str):
	"""
	Template: template_description
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'vcf': outgroupVcf,
		   	  'bed': outgroupSitesWithinCoverageBed}
	outputs = {'annotation': f'{outputDirectory}/ancestral_allele/{outputName}.annotation.variant.tsv.gz'}
	options = {
		'cores': 20,
		'memory': '30g',
		'walltime': '15:00:00'
	}
	spec = f"""
	# Sources environment
	if [ "$USER" == "jepe" ]; then
		source /home/"$USER"/.bashrc
		source activate popgen
	fi
	
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	
	[ -d {outputDirectory}/ancestral_allele ] || mkdir -p {outputDirectory}/ancestral_allele
	
	cat \\
		<(bcftools filter \\
			--threads {options['cores']} \\
			--SnpGap 5:indel \\
			--output-type u \\
			{outgroupVcf} \\
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
			--output-type u \\
			- \\
		| bcftools view \\
			--threads {options['cores']} \\
			--min-af 1 \\
			--max-af 1 \\
			--output-type u \\
			- \\
		| bcftools view \\
			--threads {options['cores']} \\
			--targets-file {outgroupSitesWithinCoverageBed} \\
			--output-type u \\
			- \\
		| bcftools query \\
			--format '%CHROM\\t%POS\\t%ALT\\n' \\
			-) \\
		<(bcftools view \\
			--threads {options['cores']} \\
			--min-af 0.00000001 \\
			--max-af 0.99999999 \\
			--output-type u \\
			{outgroupVcf} \\
		| bcftools view \\
			--threads {options['cores']} \\
			--targets-file {outgroupSitesWithinCoverageBed} \\
			--output-type u \\
			- \\
		| bcftools query \\
			--format '%CHROM\\t%POS\t.\\n' \\
			-) \\
	| sort \\
		--field-separator=$'\\t' \\
		--key=1,1 \\
		--key=2,2g \\
		- \\
	| bgzip \\
		--threads {options['cores']} \\
		--output {outputDirectory}/ancestral_allele/{outputName}.annotation.variant.prog.tsv.gz \\
		-

	mv {outputDirectory}/ancestral_allele/{outputName}.annotation.variant.prog.tsv.gz {outputs['annotation']}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def ancestral_allele_inference_merge(variantAnnotation: str, referenceAnnotation: str, outputDirectory: str, outputName: str):
	"""
	Template: template_description
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'variant': variantAnnotation,
		   	  'reference': referenceAnnotation}
	outputs = {'annotation': f'{outputDirectory}/ancestral_allele/{outputName}.annotation.tsv.gz',
			   'index': f'{outputDirectory}/ancestral_allele/{outputName}.annotation.tsv.gz.tbi'}
	options = {
		'cores': 1,
		'memory': '60g',
		'walltime': '12:00:00'
	}
	spec = f"""
	# Sources environment
	if [ "$USER" == "jepe" ]; then
		source /home/"$USER"/.bashrc
		source activate popgen
	fi
	
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	
	[ -d {outputDirectory}/ancestral_allele ] || mkdir -p {outputDirectory}/ancestral_allele
	
	awk \\
		'BEGIN{{
			FS = OFS = "\\t"
		}}
		{{
			if (FNR == NR)
			{{
				originalArray[$1, $2] = $3
				next
			}}
			if (($1, $2) in originalArray)
			{{
				print $1, $2, originalArray[$1, $2]
				next
			}}
			print $1, $2, $3
		}}' \\
		{'<(zcat ' + variantAnnotation + ')' if variantAnnotation.endswith('.gz') else variantAnnotation} \\
		{'<(zcat ' + referenceAnnotation + ')' if referenceAnnotation.endswith('.gz') else referenceAnnotation} \\
	| bgzip \\
		--threads {options['cores']} \\
		--output {outputDirectory}/ancestral_allele/{outputName}.annotation.prog.tsv.gz \\
		-

	tabix \\
		--threads {options['cores']} \\
		--sequence 1 \\
		--begin 2 \\
		--end 2 \\
		{outputDirectory}/ancestral_allele/{outputName}.annotation.prog.tsv.gz
	
	mv {outputDirectory}/ancestral_allele/{outputName}.annotation.prog.tsv.gz {outputs['annotation']}
	mv {outputDirectory}/ancestral_allele/{outputName}.annotation.prog.tsv.gz.tbi {outputs['index']}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def vcf_to_bed(vcfFile: str, outputDirectory: str):
	"""
	Template: template_description
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	filename = os.path.basename(os.path.splitext(os.path.splitext(vcfFile)[0])[0]) if vcfFile.endswith('.gz') else os.path.basename(os.path.splitext(vcfFile)[0])
	inputs = {'vcf': {vcfFile}}
	outputs = {'bed': f'{outputDirectory}/bed/{filename}.bed'}
	options = {
		'cores': 1,
		'memory': '10g',
		'walltime': '10:00:00'
	}
	spec = f"""
	# Sources environment
	if [ "$USER" == "jepe" ]; then
		source /home/"$USER"/.bashrc
		source activate popgen
	fi
	
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	
	[ -d {outputDirectory}/bed ] || mkdir -p {outputDirectory}/bed
	
	bedtools merge \\
		-i <(bcftools query \\
			--format '%CHROM\t%POS0\t%END\n' \\
			{vcfFile}) \\
		> {outputDirectory}/bed/{filename}.prog.bed
	
	mv {outputDirectory}/bed/{filename}.prog.bed {outputs['bed']}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def update_ancetral_allele_information(vcfFile: str, ancestralAnnotationFile: str, outputDirectory: str):
	"""
	Template: template_description
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	filename = os.path.basename(os.path.splitext(os.path.splitext(vcfFile)[0])[0]) if vcfFile.endswith('.gz') else os.path.basename(os.path.splitext(vcfFile)[0])
	inputs = {'vcf': vcfFile,
		   	  'annotation': ancestralAnnotationFile}
	outputs = {'vcf': f'{outputDirectory}/{filename}.aa.vcf.gz',
			   'index': f'{outputDirectory}/{filename}.aa.vcf.gz.csi'}
	options = {
		'cores': 30,
		'memory': '20g',
		'walltime': '12:00:00'
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
	
	bcftools annotate \\
		--threads {options['cores']} \\
		--header-line '##INFO=<ID=AA,Number=1,Type=String,Description="Inferred ancestral allele">' \\
		--annotations {ancestralAnnotationFile} \\
		--columns CHROM,POS,.INFO/AA \\
		--mark-sites -AA=. \\
		--output-type v \\
		{vcfFile} \\
	| awk \\
		'BEGIN{{
			FS = OFS = "\\t"
		}}
		{{
			if ($0 == "##INFO=<ID=AA=.,Number=0,Type=Flag,Description=\\"Sites not listed in AA=.\\">")
			{{
				next
			}}
			print $0
		}}' \\
		- \\
	| bcftools view \\
		--threads {options['cores']} \\
		--output-type z \\
		--output {outputDirectory}/{filename}.aa.prog.vcf.gz \\
		--write-index \\
		-
	
	mv {outputDirectory}/{filename}.aa.prog.vcf.gz {outputs['vcf']}
	mv {outputDirectory}/{filename}.aa.prog.vcf.gz.csi {outputs['index']}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)