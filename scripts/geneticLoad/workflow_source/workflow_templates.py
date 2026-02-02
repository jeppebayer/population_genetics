#!/bin/env python3
from gwf import AnonymousTarget
from gwf.executors import Conda
import os, yaml, gzip

########################## Functions ##########################

def speciesAbbreviation(speciesName: str) -> str:
	"""Creates species abbreviation from species name.

	:param str speciesName:
		Species name written as *genus* *species*"""
	genus, species = speciesName.replace(' ', '_').split('_')
	genus = genus[0].upper() + genus[1:3]
	species = species[0].upper() + species[1:3]
	return genus + species

def sequenceNamesFasta(fastaFile: str):
	"""
	Parses :format:`FASTA` file returning all sequence names in a list.
	
	:param str fastaFile:
		Sequence file in :format:`FASTA` format.
	"""
	with open(fastaFile, 'r') as fasta:
		sequenceNameList = [entry.strip().split(" ", 1)[0][1:] for entry in fasta if entry.startswith(">")]
	return sequenceNameList

def sampleNamesVCF(vcfFile: str):
	"""
	Return list of sample names in :format:`VCF` file.

	:param str vcfFile:
		Variant file in :format:`VCF` format.
	"""
	if vcfFile.endswith(".gz"):
		with gzip.open(vcfFile, 'rt') as vcf:
			for entry in vcf:
				if entry.startswith('#CHROM'):
					sampleList = [sampleName.strip() for sampleName in entry.split('\t')[10:]]
					break
	else:
		with open(vcfFile, 'r') as vcf:
			for entry in vcf:
				if entry.startswith('#CHROM'):
					sampleList = [sampleName.strip() for sampleName in entry.split('\t')[10:]]
					break
	return sampleList

########################## Templates ##########################

def snpeff_results(vcfFile: str, gtfAnnotationFile: str, speciesName: str, outputDirectory: str, environment: str, group: str | None = None):
	"""
	Template: template_description
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	filename = os.path.basename(os.path.splitext(os.path.splitext(vcfFile)[0])[0]) if vcfFile.endswith('.gz') else os.path.basename(os.path.splitext(vcfFile)[0])
	outputDirectory = f'{outputDirectory}/functional_effect_classes'
	inputs = {'vcf': vcfFile,
		   	  'gtf': gtfAnnotationFile}
	outputs = {'impact': f'{outputDirectory}/{filename}.impact.tsv'}
	options = {
		'cores': 5,
		'memory': '20g',
		'walltime': '04:00:00'
	}
	spec = f"""
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	
	[ -d {outputDirectory} ] || mkdir -p {outputDirectory}
	
	if [[ $(bcftools view -h {vcfFile} | grep "##INFO=<ID=AA,") ]]; then
		formatString="%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/AA\\t%INFO/ANN\\t%INFO/LOF\\t%INFO/NMD[\\t%SAMPLE:%GT]\\n"
	else
		formatString="%CHROM\\t%POS\\t%REF\\t%ALT\\t.\\t%INFO/ANN\\t%INFO/LOF\\t%INFO/NMD[\\t%SAMPLE:%GT]\\n"
	fi

	bcftools query \\
		--format "$formatString" \\
		--regions-file <(awk \\
			'BEGIN{{
				FS = OFS = "\\t"
			}}
			{{
				if ($3 == "CDS")
				{{
					print $1, $4 - 1, $5
				}}
			}}' \\
			{gtfAnnotationFile}) \\
		{vcfFile} \\
	| awk \\
	    -v speciesName="{speciesName}" \\
		'BEGIN{{
			FS = OFS = "\\t"
	        chromField = 1
			posField =  2
			refField = 3
			altField = 4
			aaField = 5
			annField = 6
			lofField = 7
			nmdField = 8
			genotypeFieldStart = 9
			print "sample", "species", "chromosome", "position", "reference", "alternative", "ancestral", "geneName", "geneId", "effect", "impact", "frequency", "lofPoT", "lofNoT", "nmdPoT", "nmdNoT", "exonNum_exonTotal", "nAnnotations", "annotationRank", "geneRank" > '{outputDirectory}/{filename}.impact.prog.tsv'
		}}
		{{
			if ($4 == ".")
			{{
				next
			}}

	        gsub(/\\)/, "", $lofField)
			lenLofArray = split($lofField, lofArray, "|")
			if (lenLofArray == 1)
			{{
				lofGeneCheckName = "NA"
			}}
			else
			{{
				lofGeneCheckName = lofArray[1]
			}}

			gsub(/\\)/, "", $nmdField)
			lenNmdArray = split($nmdField, nmdArray, "|")
			if (lenNmdArray == 1)
			{{
				nmdGeneCheckName = "NA"
			}}
			else
			{{
				nmdGeneCheckName = nmdArray[1]
			}}

			nAnnotations = split($annField, annotationsArray, ",")
	        lastGene = ""

			for (annNum = 1; annNum <= nAnnotations; annNum++)
			{{
				split(annotationsArray[annNum], fieldsArray, "|")
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

				if (geneName != lastGene)
				{{
					lastGene = geneName
					inGeneCount = 1
				}}
				else
				{{
					inGeneCount++
				}}

				if (geneName != lofGeneCheckName)
				{{
					lofPoT = 0.00
					lofNoT = 0
					lofGeneName = "NA"
					lofGeneId = "NA"
				}}
				else
				{{
					lofPoT = lofArray[4]
					lofNoT = lofArray[3]
					lofGeneName = lofArray[1]
					lofGeneId = lofArray[2]
				}}

				if (geneName != nmdGeneCheckName)
				{{
					nmdPoT = 0.00
					nmdNoT = 0
					nmdGeneName = "NA"
					nmdGeneId = "NA"
				}}
				else
				{{
					nmdPoT = nmdArray[4]
					nmdNoT = nmdArray[3]
					nmdGeneName = nmdArray[1]
					nmdGeneId = nmdArray[2]
				}}

				for (currentGenotypeField = genotypeFieldStart; currentGenotypeField <= NF; currentGenotypeField++)
				{{
					split($currentGenotypeField, sampleGenotype, ":")
					alleleSum = split(sampleGenotype[2], zeroesAndOnes, "/")
					frequency = 0
					for (haplotype in zeroesAndOnes)
					{{
						frequency += zeroesAndOnes[haplotype]
					}}
					if (frequency > 0)
					{{
						frequency = frequency / alleleSum
					}}
					else
					{{
						next
					}}
					print sampleGenotype[1], speciesName, $chromField, $posField, toupper($refField), toupper($altField), toupper($aaField), geneName, geneId, effect, impact, frequency, lofPoT, lofNoT, nmdPoT, nmdNoT, rankTotal, nAnnotations, annNum, inGeneCount > '{outputDirectory}/{filename}.impact.prog.tsv'
				}}
			}}
		}}' \\
		-
	
	mv {outputDirectory}/{filename}.impact.prog.tsv {outputs['impact']}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec, executor=Conda(environment), group=group)

def vcf_reformat(vcfFile: str, outputDirectory: str, environment: str, group: str | None = None):
	"""
	Template: template_description
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	filename = os.path.basename(os.path.splitext(os.path.splitext(vcfFile)[0])[0]) if vcfFile.endswith('.gz') else os.path.basename(os.path.splitext(vcfFile)[0])
	compression = f'<(zcat {vcfFile})' if vcfFile.endswith('.gz') else vcfFile
	outputDirectory = f'{outputDirectory}/reformat_vcf'
	inputs = {'vcf': vcfFile}
	outputs = {'vcf': f'{outputDirectory}/{filename}.reformat.gz',
			   'index': f'{outputDirectory}/{filename}.reformat.gz.csi'}
	options = {
		'cores': 5,
		'memory': '10g',
		'walltime': '20:00:00'
	}
	spec = f"""
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	
	[ -d {outputDirectory} ] || mkdir -p {outputDirectory}
	
	awk \\
		'BEGIN{{
			FS = OFS = "\\t"
		}}
		{{
			if ($0 !~ /^#/)
			{{
				CHROM = $1
				POS = $2
				ID = $3
				REF = $4
				ALT = $5
				QUAL = $6
				FILTER = $7
				INFO = $8
				FORMAT = $9
				sampleFieldStart = 10
				formatLength = split(FORMAT, formatField, ":")
				sampleFields = ""
				for (i = sampleFieldStart; i <= NF; i++)
				{{
					split($i, sampleCurrent, ":")
					alleleSum = split(sampleCurrent[1], zeroesAndOnes, "/")
					alternateCount = 0
					for (j in zeroesAndOnes)
					{{
						alternateCount += zeroesAndOnes[j]
					}}
					if (formatLength == 7)
					{{
						if (sampleFields)
						{{
							sampleFields = sampleFields "\\t" sampleCurrent[1] ":" alleleSum ":" (alleleSum - alternateCount) "," alternateCount ":" (alleleSum - alternateCount) ":" sampleCurrent[5] ":" alternateCount ":" sampleCurrent[7]
						}}
						else
						{{
							sampleFields = sampleCurrent[1] ":" alleleSum ":" (alleleSum - alternateCount) "," alternateCount ":" (alleleSum - alternateCount) ":" sampleCurrent[5] ":" alternateCount ":" sampleCurrent[7]
						}}
					}}
					if (formatLength == 5)
					{{
						if (sampleFields)
						{{
							sampleFields = sampleFields "\\t" sampleCurrent[1] ":" alleleSum ":" (alleleSum - alternateCount) ":" (alleleSum - alternateCount) ":" sampleCurrent[5]
						}}
						else
						{{
							sampleFields = sampleCurrent[1] ":" alleleSum ":" (alleleSum - alternateCount) ":" (alleleSum - alternateCount) ":" sampleCurrent[5]
						}}
					}}
				}}
				print CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, sampleFields
				next
			}}
			print $0
		}}' \\
		{compression} \\
	| bcftools view \\
		--threads {options['cores']} \\
		--output-type z \\
		--write-index \\
		--output {outputDirectory}/{filename}.reformat.prog.gz \\
		- 
	
	mv {outputDirectory}/{filename}.reformat.prog.gz {outputs['vcf']}
	mv {outputDirectory}/{filename}.reformat.prog.gz.csi {outputs['index']}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec, executor=Conda(environment), group=group)

def snpgenie_sequence(referenceGenome: str, gtfFile: str, vcfFile: str, sampleName: str, region: str, outputDirectory: str, environment: str, group: str | None = None):
	"""
	Template: template_description
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	tmpDirectory = f'{outputDirectory}/snpgenie/tmp/{sampleName}/{region}'
	outputDirectory = f'{outputDirectory}/snpgenie/{sampleName}/{region}'
	inputs = {'reference': referenceGenome,
		   	  'gtf': gtfFile,
			  'vcf': vcfFile}
	outputs = {'plusCodon': f'{outputDirectory}/plus/{sampleName}.{region}.codon_results.txt',
			   'plusSummary': f'{outputDirectory}/plus/{sampleName}.{region}.population_summary.txt',
			   'plusProduct': f'{outputDirectory}/plus/{sampleName}.{region}.product_results.txt',
			   'plusLog': f'{outputDirectory}/plus/{sampleName}.{region}.SNPGenie_LOG.txt',
			   'plusParameters': f'{outputDirectory}/plus/{sampleName}.{region}.SNPGenie_parameters.txt',
			   'minusCodon': f'{outputDirectory}/minus/{sampleName}.{region}.codon_results.txt',
			   'minusSummary': f'{outputDirectory}/minus/{sampleName}.{region}.population_summary.txt',
			   'minusProduct': f'{outputDirectory}/minus/{sampleName}.{region}.product_results.txt',
			   'minusLog': f'{outputDirectory}/minus/{sampleName}.{region}.SNPGenie_LOG.txt',
			   'minusParameters': f'{outputDirectory}/minus/{sampleName}.{region}.SNPGenie_parameters.txt'}
	options = {
		'cores': 10,
		'memory': '30g',
		'walltime': '06:00:00'
	}
	spec = f"""
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	
	[ -d {outputDirectory} ] && rm -rf {outputDirectory}
	[ -d {outputDirectory} ] || mkdir -p {outputDirectory}
	[ -d {tmpDirectory} ] || mkdir -p {tmpDirectory}
	
	cd {tmpDirectory}

	bcftools view \\
		--threads {options['cores']} \\
		--samples {sampleName} \\
		--targets {region} \\
		--output-type v \\
		--output {tmpDirectory}/{sampleName}.{region}.vcf \\
		--trim-alt-alleles \\
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
		> {tmpDirectory}/{os.path.splitext(os.path.basename(referenceGenome))[0]}.{region}.fasta

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
		> {tmpDirectory}/{os.path.splitext(os.path.basename(gtfFile))[0]}.{region}.gtf

	vcf2revcom.pl \\
		{tmpDirectory}/{sampleName}.{region}.vcf \\
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
			{referenceGenome})


	fasta2revcom.pl \\
		{tmpDirectory}/{os.path.splitext(os.path.basename(referenceGenome))[0]}.{region}.fasta

	gtf2revcom.pl \\
		{tmpDirectory}/{os.path.splitext(os.path.basename(gtfFile))[0]}.{region}.gtf \\
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
			{referenceGenome})

	echo -e "#########################\\n# Processing '+' strand #\\n#########################"

	snpgenie.pl \\
		--vcfformat=4 \\
		--snpreport={tmpDirectory}/{sampleName}.{region}.vcf \\
		--fastafile={tmpDirectory}/{os.path.splitext(os.path.basename(referenceGenome))[0]}.{region}.fasta \\
		--gtffile={tmpDirectory}/{os.path.splitext(os.path.basename(gtfFile))[0]}.{region}.gtf \\
		--minfreq=0 \\
		--slidingwindow=9 \\
		--workdir={tmpDirectory} \\
		--outdir={outputDirectory}/plus

	echo -e "#########################\\n# Processing '-' strand #\\n#########################"

	snpgenie.pl \\
		--vcfformat=4 \\
		--snpreport={tmpDirectory}/{sampleName}.{region}_revcom.vcf \\
		--fastafile={tmpDirectory}/{os.path.splitext(os.path.basename(referenceGenome))[0]}.{region}_revcom.fasta \\
		--gtffile={tmpDirectory}/{os.path.splitext(os.path.basename(gtfFile))[0]}.{region}_revcom.gtf \\
		--minfreq=0 \\
		--slidingwindow=9 \\
		--workdir={tmpDirectory} \\
		--outdir={outputDirectory}/minus
	
	rm -rf {tmpDirectory}

	mv {outputDirectory}/plus/codon_results.txt {outputs['plusCodon']}
	mv {outputDirectory}/plus/population_summary.txt {outputs['plusSummary']}
	mv {outputDirectory}/plus/product_results.txt {outputs['plusProduct']}
	mv {outputDirectory}/plus/SNPGenie_LOG.txt {outputs['plusLog']}
	mv {outputDirectory}/plus/SNPGenie_parameters.txt {outputs['plusParameters']}

	mv {outputDirectory}/minus/codon_results.txt {outputs['minusCodon']}
	mv {outputDirectory}/minus/population_summary.txt {outputs['minusSummary']}
	mv {outputDirectory}/minus/product_results.txt {outputs['minusProduct']}
	mv {outputDirectory}/minus/SNPGenie_LOG.txt {outputs['minusLog']}
	mv {outputDirectory}/minus/SNPGenie_parameters.txt {outputs['minusParameters']}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec, executor=Conda(environment), group=group)

def snpgenie_summarize_results(outputDirectory: str, environment: str, group: str | None = None):
	"""
	Template: template_description
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	outputDirectory = f'{outputDirectory}'
	inputs = {}
	outputs = {}
	options = {
		'cores': 1,
		'memory': '20g',
		'walltime': '04:00:00'
	}
	spec = f"""
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	
	[ -d {outputDirectory} ] || mkdir -p {outputDirectory}
	
	
	
	mv
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec, executor=Conda(environment), group=group)