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