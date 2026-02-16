#!/bin/awk -f

# Takes bcftools query input using format string "%CHROM\t%POS\t%REF\t%ALT\t%INFO/AA\t%INFO/ANN\t%INFO/LOF\t%INFO/NMD[\t%SAMPLE:%GT]\n"

BEGIN{
	FS = OFS = "\t"
	# Associate each field with a variable name
	chromField = 1
	posField =  2
	refField = 3
	altField = 4
	aaField = 5
	annField = 6
	lofField = 7
	nmdField = 8
	genotypeFieldStart = 9
	# Print output file header
	print "sample", "species", "chromosome", "position", "reference", "alternative", "ancestral", "geneName", "geneId", "effect", "impact", "frequency", "lofPoT", "lofNoT", "nmdPoT", "nmdNoT", "exonNum_exonTotal", "nAnnotations", "annotationRank", "geneRank"
}
{
	# If there is no variant at site, skip
	if ($altField == ".")
	{
		next
	}

	# Remove parenthesis from LoF field
	gsub(/\)/, "", $lofField)
	gsub(/\(/, "", $lofField)
	# Split LoF field, separating LoF annotations
	lenLofArray = split($lofField, lofArray, ",")
	# Make sure lofGeneCheckArray is empty
	split("", lofGeneCheckArray)
	# Do for each LoF annotation
	for (lofGene in lofArray)
	{
		# Split LoF annotation into each of its components: gene, geneID, number of transcripts affected, proportion of transcripts affected
		lenLofAnn = split(lofArray[lofGene], lofAnn, "|")
		# If length of split is 1 there is no LoF annotation
		if (lenLofAnn == 1)
		{
			lofGeneCheckArray["NA"] = "NA|NA|0|0.00"
		}
		# If LoF annotation is present store it in an array based on gene name
		else
		{
			lofGeneCheckArray[lofAnn[1]] = lofArray[lofGene]
		}
	}

	# Remove parenthesis from NMD field
	gsub(/\)/, "", $nmdField)
	gsub(/\(/, "", $nmdField)
	# Split NMD field, separating NMD annotations
	lenNmdArray = split($nmdField, nmdArray, ",")
	# Make sure nmdGeneCheckArray is empty
	split("", nmdGeneCheckArray)
	# Do for each NMD annotation
	for (nmdGene in nmdArray)
	{
		# Split NMD annotation into each of its components: gene, geneID, number of transcripts affected, proportion of transcripts affected
		lenNmdAnn = split(nmdArray[nmdGene], nmdAnn, "|")
		if (lenNmdAnn == 1)
		{
			nmdGeneCheckArray["NA"] = "NA|NA|0|0.00"
		}
		# If NMD annotation is present store it in an array based on gene name
		else
		{
			nmdGeneCheckArray[nmdAnn[1]] = nmdArray[nmdGene]
		}
	}

	# Split SnpEff annotation field, separating each SnpEff annotation
	nAnnotations = split($annField, annotationsArray, ",")
	# Create empty lastGene variable
	lastGene = ""
	# Do for each SnpEff annotation
	for (annNum = 1; annNum <= nAnnotations; annNum++)
	{
		# Split annotation into each of its components
		split(annotationsArray[annNum], fieldsArray, "|")
		# Associate each component with a variable for ease
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

		# If the impact of the annotation is not LOW, MODERATE or HIGH disregard it
		if (impact != "LOW" && impact != "MODERATE" && impact != "HIGH")
		{
			continue
		}
		# Should the component denoting exon or intron rank out of the total number of exons or introns be empty assign it a value of NA
		if (length(rankTotal) == 0)
		{
			rankTotal = "NA"
		}
		# If the current gene name is not the same as the last gene name, update name of last gene and reset within gene annotation counter
		if (geneName != lastGene)
		{
			lastGene = geneName
			inGeneCount = 1
		}
		# If current gene name is the same as the last gene name, increment within gene annotation counter
		else
		{
			inGeneCount++
		}
		# Do for each registered LoF annotated gene
	    for (gene in lofGeneCheckArray)
	    {
			# If a LoF gene is found to match the current gene, retrieve LoF annotation components from array
	    	if (gene == geneName)
	    	{
	    		split(lofGeneCheckArray[gene], currentLofGeneArray, "|")
	    		lofPoT = currentLofGeneArray[4]
	    		lofNoT = currentLofGeneArray[3]
	    		lofGeneName = currentLofGeneArray[1]
	    		lofGeneId = currentLofGeneArray[2]
	    		break	
	    	}
			# If there is no LoF gene on current position, retrieve NA values from array
	    	if (gene == "NA")
	    	{
	    		split(lofGeneCheckArray[gene], currentLofGeneArray, "|")
	    		lofPoT = currentLofGeneArray[4]
	    		lofNoT = currentLofGeneArray[3]
	    		lofGeneName = currentLofGeneArray[1]
	    		lofGeneId = currentLofGeneArray[2]
	    		break
	    	}
	    }
		# Do for each registered NMD annotated gene
		for (gene in nmdGeneCheckArray)
		{
			# If a NMD gene is found to match the current gene, retrieve NMD annotation components from array
			if (gene == geneName)
			{
				split(nmdGeneCheckArray[gene], currentNmdGeneArray, "|")
				nmdPoT = currentNmdGeneArray[4]
				nmdNoT = currentNmdGeneArray[3]
				nmdGeneName = currentNmdGeneArray[1]
				nmdGeneId = currentNmdGeneArray[2]
				break	
			}
			# If there is no NMD gene on current position, retrieve NA values from array
			if (gene == "NA")
			{
				split(nmdGeneCheckArray[gene], currentNmdGeneArray, "|")
				nmdPoT = currentNmdGeneArray[4]
				nmdNoT = currentNmdGeneArray[3]
				nmdGeneName = currentNmdGeneArray[1]
				nmdGeneId = currentNmdGeneArray[2]
				break
			}
		}
		# Do for each sample-genotype field
		for (currentGenotypeField = genotypeFieldStart; currentGenotypeField <= NF; currentGenotypeField++)
		{
			# Split sample-genotype field into its components
			split($currentGenotypeField, sampleGenotype, ":")
			# Split genotype notation into haplotypes (0's and 1's). alleleSum is equivalent to number of chromosomes
			alleleSum = split(sampleGenotype[2], zeroesAndOnes, "/")
			# Set initial frequency to 0
			frequency = 0
			# Do for each haplotype
			for (haplotype in zeroesAndOnes)
			{
				# Add haplotype value to frequency
				frequency += zeroesAndOnes[haplotype]
			}
			# If frequency is greater than 0, convert from total to proportion
			if (frequency > 0)
			{
				frequency = frequency / alleleSum
			}
			# If frequency is 0 disregard sample
			else
			{
				continue
			}
			# Print values to output file
			print sampleGenotype[1], speciesName, $chromField, $posField, toupper($refField), toupper($altField), toupper($aaField), geneName, geneId, effect, impact, frequency, lofPoT, lofNoT, nmdPoT, nmdNoT, rankTotal, nAnnotations, annNum, inGeneCount
		}
	}
}