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

def snpeff_results(vcfFile: str, gtfAnnotationFile: str, outputDirectory: str):
	"""
	Template: template_description
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	filename = os.path.basename(os.path.splitext(os.path.splitext(vcfFile)[0])[0]) if vcfFile.endswith('.gz') else os.path.basename(os.path.splitext(vcfFile)[0])
	inputs = {'vcf': vcfFile,
		   	  'gtf': gtfAnnotationFile}
	outputs = {'impactFreqs': f'{outputDirectory}/{filename}.impactFreqs.tsv'}
	options = {
		'cores': 10,
		'memory': '20g',
		'walltime': '02:00:00'
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
	
	bcftools query \
		--format "%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/AA\\t%INFO/ANN\\t%INFO/LOF\\t%INFO/NMD[\\t%SAMPLE:%GT]\\n" \
		--regions-file <(awk \
			'BEGIN{{
				FS = OFS = "\\t"
			}}
			{{
				if ($3 == "CDS")
				{{
					print $1, $4 - 1, $5
				}}
			}}' \
			{gtfAnnotationFile}) \
		{vcfFile} \
	| awk \
		'BEGIN{{
			FS = OFS = "\\t"
			annField = 6
			lofField = 7
			nmdField = 8
			genotypeFieldStart = 9
			print "sample", "chromosome", "position", "reference", "alternative", "ancestral", "gene", "exonNum/exonTotal", "effect", "impact", "frequency", "nAnnotations", "lof" > "{outputDirectory}/{filename}.impactFreqs.prog.tsv"
		}}
		{{
			nAnnotations = split($annField, annotationsArray, ",")
			split(annotationsArray[1], fieldsArray, "|")
			allele = fieldsArray[1]
			effect = fieldsArray[2]
			impact = fieldsArray[3]
			geneName = fieldsArray[4]
			geneId = fieldsArray[5]
			featureType = fieldsArray[6]
			featureId = fieldsArray[7]
			transcriptBiotype = fieldsArray[8]
			rankTotal = fieldsArray[9]
			hgvsC = fieldsArray[10]
			hgvsP = fieldsArray[11]
			cdnaPositionLength = fieldsArray[12]
			cdsPositionLength = fieldsArray[13]
			proteinPositionLength = fieldsArray[14]
			distanceToFeature = fieldsArray[15]
			warnings = fieldsArray[16]
			if (impact != "LOW" && impact != "MODERATE" && impact != "HIGH")
			{{
				next
			}}
			if (length(rankTotal) == 0)
			{{
				rankTotal = "NA"
			}}
			gsub(/\)/, "", $lofField)
			lenLofArray = split($lofField, lofArray, "|")
			if (lenLofArray == 1)
			{{
				lof = "0.00"
			}}
			else
			{{
				lof = lofArray[4]
			}}
			for (i = genotypeFieldStart; i <= NF; i++)
			{{
				split($i, sampleGenotype, ":")
				split(sampleGenotype[2], zeroesAndOnes, "/")
				frequency = 0
				for (j in zeroesAndOnes)
				{{
					frequency += zeroesAndOnes[j]
				}}
				if (frequency > 0)
				{{
					frequency = frequency / 100
				}}
				print sampleGenotype[1], $1, $2, touppper($3), toupper($4), toupper($5), geneId, rankTotal, effect, impact, frequency, nAnnotations, lof | "sort -k1,1 -k2,2 -k3,3n >> {outputDirectory}/{filename}.impactFreqs.prog.tsv"
			}}
		}}' \
		-
	
	mv {outputDirectory}/{filename}.impactFreqs.prog.tsv {outputs['impactFreqs']}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def snpgenie(referenceGenome: str, gtfFile: str, vcfFile: str, sampleName: str, region: str, outputDirectory: str):
	"""
	Template: template_description
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'reference': referenceGenome,
		   	  'gtf': gtfFile,
			  'vcf': vcfFile}
	outputs = {'plusCodon': f'{outputDirectory}/{sampleName}/{region}/plus/{sampleName}.{region}.codon_results.txt',
			   'plusSummary': f'{outputDirectory}/{sampleName}/{region}/plus/{sampleName}.{region}.population_summary.txt',
			   'plusProduct': f'{outputDirectory}/{sampleName}/{region}/plus/{sampleName}.{region}.product_results.txt',
			   'plusLog': f'{outputDirectory}/{sampleName}/{region}/plus/{sampleName}.{region}.SNPGenie_LOG.txt',
			   'plusParameters': f'{outputDirectory}/{sampleName}/{region}/plus/{sampleName}.{region}.SNPGenie_parameters.txt',
			   'minusCodon': f'{outputDirectory}/{sampleName}/{region}/minus/{sampleName}.{region}.codon_results.txt',
			   'minusSummary': f'{outputDirectory}/{sampleName}/{region}/minus/{sampleName}.{region}.population_summary.txt',
			   'minusProduct': f'{outputDirectory}/{sampleName}/{region}/minus/{sampleName}.{region}.product_results.txt',
			   'minusLog': f'{outputDirectory}/{sampleName}/{region}/minus/{sampleName}.{region}.SNPGenie_LOG.txt',
			   'minusParameters': f'{outputDirectory}/{sampleName}/{region}/minus/{sampleName}.{region}.SNPGenie_parameters.txt'}
	options = {
		'cores': 10,
		'memory': '50g',
		'walltime': '24:00:00'
	}
	spec = f"""
	# Sources environment
	if [ "$USER" == "jepe" ]; then
		source /home/"$USER"/.bashrc
		source activate snpgenie
	fi
	
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	
	[ -d {outputDirectory}/tmp/{sampleName}/{region} ] && rm -rf {outputDirectory}/tmp/{sampleName}/{region}
	[ -d {outputDirectory}/tmp/{sampleName}/{region} ] || mkdir -p {outputDirectory}/tmp/{sampleName}/{region}
	[ -d {outputDirectory}/{sampleName}/{region} ] && rm -rf {outputDirectory}/{sampleName}/{region}
	[ -d {outputDirectory}/{sampleName}/{region} ] || mkdir -p {outputDirectory}/{sampleName}/{region}
	
	bcftools view \\
		--threads {options['cores']} \\
		--samples {sampleName} \\
		--targets {region} \\
		--output-type v \\
		--output {outputDirectory}/tmp/{sampleName}/{region}/{sampleName}.{region}.vcf \\
		{vcfFile}

	awk \\
		-v region={region} \\
		'BEGIN{{
			RS = ">"
			ORS = ""
			FS = OFS = "\\n"
		}}
		{{
			if (NR > 1 && $1 == region)
			{{
				print ">" $0
				exit
			}}
		}}' \\
		{referenceGenome} \\
		> {outputDirectory}/tmp/{sampleName}/{region}/{os.path.splitext(os.path.basename(referenceGenome))[0]}.{region}.fasta

	awk \\
		-v region={region} \\
		'BEGIN{{
			FS = OFS = "\\t"
		}}
		{{
			if ($1 == region)
			{{
				print $0
			}}
		}}' \\
		{gtfFile} \\
		> {outputDirectory}/tmp/{sampleName}/{region}/{os.path.splitext(os.path.basename(gtfFile))[0]}.{region}.gtf

	vcf2revcom.pl \\
		{outputDirectory}/tmp/{sampleName}/{region}/{sampleName}.{region}.vcf \\
		<(awk \\
			-v region={region} \\
			'BEGIN{{
				RS = ">"
				ORS = ""
				FS = OFS = "\\n"
			}}
			{{
				if (NR > 1 && $1 == region)
				{{
					for (i = 2; i <= NF; i++)
					{{
						sumLength += length($i)
					}}
					print sumLength
					exit
				}}
			}}' \\
			{referenceGenome}) \\
		> {outputDirectory}/tmp/{sampleName}/{region}/{sampleName}.{region}.revcom.vcf


	fasta2revcom.pl \\
		{outputDirectory}/tmp/{sampleName}/{region}/{os.path.splitext(os.path.basename(referenceGenome))[0]}.{region}.fasta \\
		> {outputDirectory}/tmp/{sampleName}/{region}/{os.path.splitext(os.path.basename(referenceGenome))[0]}.{region}.revcom.fasta

	gtf2revcom.pl \\
		{outputDirectory}/tmp/{sampleName}/{region}/{os.path.splitext(os.path.basename(gtfFile))[0]}.{region}.gtf \\
		<(awk \\
			-v region={region} \\
			'BEGIN{{
				RS = ">"
				ORS = ""
				FS = OFS = "\\n"
			}}
			{{
				if (NR > 1 && $1 == region)
				{{
					for (i = 2; i <= NF; i++)
					{{
						sumLength += length($i)
					}}
					print sumLength
					exit
				}}
			}}' \
			{referenceGenome}) \
		> {outputDirectory}/tmp/{sampleName}/{region}/{os.path.splitext(os.path.basename(gtfFile))[0]}.{region}.revcom.gtf

	echo -e "#########################\\n# Processing '+' strand #\\n#########################"

	snpgenie.pl \\
		--vcfformat=4 \\
		--snpreport={outputDirectory}/tmp/{sampleName}/{region}/{sampleName}.{region}.vcf \\
		--fastafile={outputDirectory}/tmp/{sampleName}/{region}/{os.path.splitext(os.path.basename(referenceGenome))[0]}.{region}.fasta \\
		--gtffile={outputDirectory}/tmp/{sampleName}/{region}/{os.path.splitext(os.path.basename(gtfFile))[0]}.{region}.gtf \\
		--minfreq=0 \\
		--slidingwindow=9 \\
		--workdir={outputDirectory}/tmp/{sampleName}/{region} \\
		--outdir={outputDirectory}/{sampleName}/{region}/plus

	echo -e "#########################\\n# Processing '-' strand #\\n#########################"

	snpgenie.pl \\
		--vcfformat=4 \\
		--snpreport={outputDirectory}/tmp/{sampleName}/{region}/{sampleName}.{region}.revcom.vcf \\
		--fastafile={outputDirectory}/tmp/{sampleName}/{region}/{os.path.splitext(os.path.basename(referenceGenome))[0]}.{region}.revcom.fasta \\
		--gtffile={outputDirectory}/tmp/{sampleName}/{region}/{os.path.splitext(os.path.basename(gtfFile))[0]}.{region}.revcom.gtf \\
		--minfreq=0 \\
		--slidingwindow=9 \\
		--workdir={outputDirectory}/tmp/{sampleName}/{region} \\
		--outdir={outputDirectory}/{sampleName}/{region}/minus
	
	rm -rf {outputDirectory}/tmp/{sampleName}/{region}

	mv {outputDirectory}/{sampleName}/{region}/plus/codon_results.txt {outputs['plusCodon']}
	mv {outputDirectory}/{sampleName}/{region}/plus/population_summary.txt {outputs['plusSummary']}
	mv {outputDirectory}/{sampleName}/{region}/plus/product_results.txt {outputs['plusProduct']}
	mv {outputDirectory}/{sampleName}/{region}/plus/SNPGenie_LOG.txt {outputs['plusLog']}
	mv {outputDirectory}/{sampleName}/{region}/plus/SNPGenie_parameters.txt {outputs['plusParameters']}

	mv {outputDirectory}/{sampleName}/{region}/minus/codon_results.txt {outputs['minusCodon']}
	mv {outputDirectory}/{sampleName}/{region}/minus/population_summary.txt {outputs['minusSummary']}
	mv {outputDirectory}/{sampleName}/{region}/minus/product_results.txt {outputs['minusProduct']}
	mv {outputDirectory}/{sampleName}/{region}/minus/SNPGenie_LOG.txt {outputs['minusLog']}
	mv {outputDirectory}/{sampleName}/{region}/minus/SNPGenie_parameters.txt {outputs['minusParameters']}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)