#!/bin/env python3
from gwf import AnonymousTarget
import os, glob

########################## Functions ##########################

def species_abbreviation(species_name: str) -> str:
	"""Creates species abbreviation from species name.

	:param str species_name:
		Species name written as *genus* *species*"""
	genus, species = species_name.replace(' ', '_').split('_')
	genus = genus[0].upper() + genus[1:3]
	species = species[0].upper() + species[1:3]
	return genus + species

def sequence_names_fasta(fasta_file: str):
    """
    Parses :format:`FASTA` file returning all sequence names in a list.
    
    ::
    
        return [str, ...]
    
    :param str fasta_file:
        Sequence file in :format:`FASTA` format.
    """
    fasta_list = []
    seq_name = None
    with open(fasta_file, 'r') as fasta:
        for entry in fasta:
            entry = entry.strip()
            if entry.startswith(">"):
                if seq_name:
                    fasta_list.append(seq_name)
                entry = entry.split(" ", 1)
                seq_name = entry[0][1:]
        fasta_list.append(seq_name)
    return fasta_list

def parse_fasta(fasta_file: str):
    """
    Parses :format:`FASTA` file returning all sequence names and lengths paired in a list of dictionaries.
    
    ::
    
        return [{'sequence_name': str, 'sequence_length': int}, ...]
    
    :param str fasta_file:
        Sequence file in :format:`FASTA` format.
    """
    fasta_list = []
    seq_name = None
    length = 0
    with open(fasta_file, 'r') as fasta:
        for entry in fasta:
            entry = entry.strip()
            if entry.startswith(">"):
                if seq_name:
                    fasta_list.append({'sequence_name': seq_name, 'sequence_length': length})
                    length = 0
                entry = entry.split(" ", 1)
                seq_name = entry[0][1:]
            else:
                length += len(entry)
        fasta_list.append({'sequence_name': seq_name, 'sequence_length': length})
    return fasta_list

########################## SnpEff ##########################

def snpeff_database_build(gtf_annotation_file: str, reference_genome_file: str, species_name: str, snpeff_directory: str = f'{os.path.dirname(os.path.realpath(__file__))}/software/snpeff'):
	"""
	Template: Constructs custom SnpEff database entry.
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'gtf': gtf_annotation_file,
		   	  'reference': reference_genome_file}
	outputs = {'sequences': f'{snpeff_directory}/data/{os.path.splitext(os.path.basename(reference_genome_file))[0].split(sep="_genomic")[0]}/sequences.fa',
			   'genes': f'{snpeff_directory}/data/{os.path.splitext(os.path.basename(reference_genome_file))[0].split(sep="_genomic")[0]}/genes.gtf',
			   'proteins': f'{snpeff_directory}/data/{os.path.splitext(os.path.basename(reference_genome_file))[0].split(sep="_genomic")[0]}/protein.fa',
			   'cds': f'{snpeff_directory}/data/{os.path.splitext(os.path.basename(reference_genome_file))[0].split(sep="_genomic")[0]}/cds.fa',
			   'predictor': f'{snpeff_directory}/data/{os.path.splitext(os.path.basename(reference_genome_file))[0].split(sep="_genomic")[0]}/snpEffectPredictor.bin'}
	options = {
		'cores': 1,
		'memory': '80g',
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
	
	[ -d {snpeff_directory}/data ] || mkdir -p {snpeff_directory}/data

	snpeffconfig={snpeff_directory}/snpEff.config

	if [ ! -e "$snpeffconfig" ]; then
		echo -e "#-------------------------------------------------------------------------------\\n#\\n# SnpEff configuration file\\n#\\n#\\n#-------------------------------------------------------------------------------\\n" > "$snpeffconfig"
		echo -e "#---\\n# Databases are stored here\\n# E.g.: Information for 'hg19' is stored in data.dir/hg19/\\n#\\n# You can use tilde ('~') as first character to refer to your home directory.\\n# Also, a non-absolute path will be relative to config's file dir\\n#\\n#---\\ndata.dir = ./data/\\n" >> "$snpeffconfig"
		echo -e "#-------------------------------------------------------------------------------\\n# Loss of function (LOF)\\n#-------------------------------------------------------------------------------\\n\\n# It is assumed that even with a protein coding change at the\\n# last 5% of the protein, the protein could still be functional.\\nlof.ignoreProteinCodingAfter  : 0.95\\n\\n# It is assumed that even with a protein coding change at the\\n# first 5% of the protein:\\n#\\t\\t\\t\\t"..suggesting some disrupted transcripts are\\n#\\t\\t\\t\\trescued by transcriptional reinitiation at an\\n#\\t\\t\\t\\talternative start codon."\\nlof.ignoreProteinCodingBefore : 0.05\\n\\n# Larger deletions removing either the first exon or more than\\n# 50% of the protein-coding sequence of the affected transcript\\nlof.deleteProteinCodingBases : 0.50\\n" >> "$snpeffconfig"
		echo -e "#-------------------------------------------------------------------------------\\n# Codon tables\\n#\\n# Format:\\tIt's a comma separated "codon/aminoAcid[+*]" list\\n#\\t\\t\\tWhere 'codon' is in uppper case, aminoAcid is a one letter\\n#\\t\\t\\tcode, '+' denotes start codon and '*' denotes stop codon.\\n#\\n# References:\\thttp://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi\\n#\\t\\t\\t\\tftp://ftp.ncbi.nih.gov/entrez/misc/data/gc.prt\\n#-------------------------------------------------------------------------------\\n\\ncodon.Standard\\t\\t\\t\\t: TTT/F, TTC/F, TTA/L, TTG/L+, TCT/S, TCC/S, TCA/S, TCG/S, TAT/Y, TAC/Y, TAA/*, TAG/*, TGT/C, TGC/C, TGA/*, TGG/W, CTT/L, CTC/L, CTA/L, CTG/L+, CCT/P, CCC/P, CCA/P, CCG/P, CAT/H, CAC/H, CAA/Q, CAG/Q, CGT/R, CGC/R, CGA/R, CGG/R, ATT/I, ATC/I, ATA/I, ATG/M+, ACT/T, ACC/T, ACA/T, ACG/T, AAT/N, AAC/N, AAA/K, AAG/K, AGT/S, AGC/S, AGA/R, AGG/R, GTT/V, GTC/V, GTA/V, GTG/V, GCT/A, GCC/A, GCA/A, GCG/A, GAT/D, GAC/D, GAA/E, GAG/E, GGT/G, GGC/G, GGA/G, GGG/G\\n" >> "$snpeffconfig"
		echo -e "#-------------------------------------------------------------------------------\\n# Databases & Genomes\\n#\\n# One entry per genome version.\\n#\\n# For genome version 'ZZZ' the entries look like\\n#\\t\\tZZZ.genome\\t\\t\\t\\t: Real name for ZZZ (e.g. 'Human')\\n#\\t\\tZZZ.reference\\t\\t\\t: [Optional] Comma separated list of URL to site/s where information for building ZZZ database was extracted.\\n#\\t\\tZZZ.chrName.codonTable\\t: [Optional] Define codon table used for chromosome 'chrName' (Default: 'codon.Standard')\\n#\\n#-------------------------------------------------------------------------------\\n" >> "$snpeffconfig"
		echo -e "#---\\n# Centre for Ecological Genetics (AU) - Custom Databases\\n#---\\n" >> "$snpeffconfig"
	fi
	
	line=$(awk '{{if ($0 ~ /^{os.path.splitext(os.path.basename(reference_genome_file))[0].split(sep="_genomic")[0]}.genome/) {{print NR; exit}}}}' {snpeff_directory}/snpEff.config)

	if [ ! -z "$line" ]; then
		sed -i "$((line - 1)),$((line + 3))d" {snpeff_directory}/snpEff.config
	fi

	echo -e "# {species_name} genome, version {os.path.splitext(os.path.basename(reference_genome_file))[0].split(sep="_genomic")[0]}" >> "$snpeffconfig"
	echo -e "{os.path.splitext(os.path.basename(reference_genome_file))[0].split(sep="_genomic")[0]}.genome : {species_name}" >> "$snpeffconfig"
	echo -e "{os.path.splitext(os.path.basename(reference_genome_file))[0].split(sep="_genomic")[0]}.file_location : {reference_genome_file}" >> "$snpeffconfig"
	echo -e "{os.path.splitext(os.path.basename(reference_genome_file))[0].split(sep="_genomic")[0]}.addition_date : $(date +%d'/'%m'/'%Y)\\n" >> "$snpeffconfig"

	[ -d {snpeff_directory}/data/{os.path.splitext(os.path.basename(reference_genome_file))[0].split(sep="_genomic")[0]} ] || mkdir -p {snpeff_directory}/data/{os.path.splitext(os.path.basename(reference_genome_file))[0].split(sep="_genomic")[0]}

	if [ ! -e {snpeff_directory}/data/{os.path.splitext(os.path.basename(reference_genome_file))[0].split(sep="_genomic")[0]}/genes.gtf ]; then
		ln \
			-s \
			{gtf_annotation_file} \
			{snpeff_directory}/data/{os.path.splitext(os.path.basename(reference_genome_file))[0].split(sep="_genomic")[0]}/genes.gtf
	fi

	if [ ! -e {snpeff_directory}/data/{os.path.splitext(os.path.basename(reference_genome_file))[0].split(sep="_genomic")[0]}/sequences.fa ]; then
		ln \
			-s \
			{reference_genome_file} \
			{snpeff_directory}/data/{os.path.splitext(os.path.basename(reference_genome_file))[0].split(sep="_genomic")[0]}/sequences.fa
	fi

	agat_sp_extract_sequences.pl \
		--gff {gtf_annotation_file} \
		--fasta {reference_genome_file} \
		--type cds \
		--output {snpeff_directory}/data/{os.path.splitext(os.path.basename(reference_genome_file))[0].split(sep="_genomic")[0]}/cds.prog.fa

	agat_sp_extract_sequences.pl \
		--gff {gtf_annotation_file} \
		--fasta {reference_genome_file} \
		--type cds \
		--protein \
		--output {snpeff_directory}/data/{os.path.splitext(os.path.basename(reference_genome_file))[0].split(sep="_genomic")[0]}/protein.prog.fa

	mv {snpeff_directory}/data/{os.path.splitext(os.path.basename(reference_genome_file))[0].split(sep="_genomic")[0]}/cds.prog.fa {snpeff_directory}/data/{os.path.splitext(os.path.basename(reference_genome_file))[0].split(sep="_genomic")[0]}/cds.fa
	mv {snpeff_directory}/data/{os.path.splitext(os.path.basename(reference_genome_file))[0].split(sep="_genomic")[0]}/protein.prog.fa {snpeff_directory}/data/{os.path.splitext(os.path.basename(reference_genome_file))[0].split(sep="_genomic")[0]}/protein.fa

	export _JAVA_OPTIONS="-Xmx{options['memory']}"

	snpEff build \
		-gtf22 \
		-config {snpeff_directory}/snpEff.config \
		-nodownload \
		-verbose \
		{os.path.splitext(os.path.basename(reference_genome_file))[0].split(sep="_genomic")[0]}
	
	rm *.agat.log

	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def name_snpeff_annotation(idx: str, target: AnonymousTarget) -> str:
	return f'snpeff_{os.path.basename(os.path.dirname(target.outputs['ann'])).replace("-", "_")}'

def snpeff_annotation(vcf_file: str, snpeff_predictor_file: str, snpeff_config_file: str, output_directory: str, sample_group: str, sample_name: str):
	"""
	Template: Annotates :format:`VCF` file with variant function.
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'vcf': vcf_file,
		   	  'predictor': snpeff_predictor_file,
			  'config': snpeff_config_file}
	outputs = {'ann': f'{output_directory}/snpEff/{sample_group}/{sample_name}/{os.path.splitext(os.path.basename(vcf_file))[0] if vcf_file.endswith(".vcf") else os.path.splitext(os.path.splitext(os.path.basename(vcf_file))[0])[0]}.ann.vcf.gz',
			   'csv': f'{output_directory}/snpEff/{sample_group}/{sample_name}/{sample_name}.snpEff_summary.csv',
			   'txt': f'{output_directory}/snpEff/{sample_group}/{sample_name}/{sample_name}.snpEff_summary.genes.txt',
			   'html': f'{output_directory}/snpEff/{sample_group}/{sample_name}/{sample_name}.snpEff_summary.html'}
	options = {
		'cores': 18,
		'memory': '80g',
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
	
	[ -d {output_directory}/snpEff/{sample_group}/{sample_name} ] || mkdir -p {output_directory}/snpEff/{sample_group}/{sample_name}
	
	export _JAVA_OPTIONS="-Xmx{options['memory']}"

	snpEff ann \
		-csvStats {output_directory}/snpEff/{sample_group}/{sample_name}/{sample_name}.snpEff_summary.prog.csv \
		-htmlStats {output_directory}/snpEff/{sample_group}/{sample_name}/{sample_name}.snpEff_summary.prog.html \
		-nodownload \
		-config {snpeff_config_file} \
		-verbose \
		-i vcf \
		-o vcf \
		{os.path.basename(os.path.dirname(snpeff_predictor_file))} \
		{vcf_file} \
	| gzip \
		--stdout \
		- \
		> {output_directory}/snpEff/{sample_group}/{sample_name}/{os.path.splitext(os.path.basename(vcf_file))[0] if vcf_file.endswith('.vcf') else os.path.splitext(os.path.splitext(os.path.basename(vcf_file))[0])[0]}.ann.prog.vcf.gz
	
	mv {output_directory}/snpEff/{sample_group}/{sample_name}/{os.path.splitext(os.path.basename(vcf_file))[0] if vcf_file.endswith('.vcf') else os.path.splitext(os.path.splitext(os.path.basename(vcf_file))[0])[0]}.ann.prog.vcf.gz {outputs['ann']}
	mv {output_directory}/snpEff/{sample_group}/{sample_name}/{sample_name}.snpEff_summary.prog.csv {outputs['csv']}
	mv {output_directory}/snpEff/{sample_group}/{sample_name}/{sample_name}.snpEff_summary.prog.genes.txt {outputs['txt']}
	mv {output_directory}/snpEff/{sample_group}/{sample_name}/{sample_name}.snpEff_summary.prog.html {outputs['html']}

	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def name_snpeff_freqs(idx: str, target: AnonymousTarget) -> str:
	return f'snpeff_freqs_{os.path.basename(os.path.dirname(target.inputs['vcf'])).replace("-", "_")}'

def snpeff_freqs(ann: str):
	"""
	Template: Summarizes annotated variant functions.
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	annotated_vcf = f'<(zcat {ann})' if ann.endswith('.gz') else ann
	inputs = {'vcf': ann}
	outputs = {'csv': f'{os.path.dirname(ann)}/effectsummary.tsv'}
	options = {
		'cores': 1,
		'memory': '12g',
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
	
	awk \
        'BEGIN {{
			FS = OFS = "\\t"
		}}
		{{
			if ($0 !~ /^#/)
			{{
				chromosome = $1
				split($8, formatfield, ";")
				split(formatfield[4], affield, "=")
				split(formatfield[42], annsection, "|")
				af = affield[2]
				effect = annsection[3]
				if (af == 0)
				{{
					next
				}}
				if (effect == "LOW" || effect == "MODERATE" || effect == "HIGH")
				{{
					if (length(frequencies[chromosome, effect]) == 0)
					{{
						frequencies[chromosome, effect] = af
					}}
					else
					{{
						frequencies[chromosome, effect] = frequencies[chromosome, effect] "," af
					}}
				}}
			}}
		}}
		END {{
			print "chromosome", "effect", "frequencies", "sum_frequencies", "num_frequencies", "mean_frequency"
			for (i in frequencies)
			{{
				split(i, chreffect, "\\034")
				split(frequencies[i], sumarray, ",")
				for (j in sumarray)
				{{
					n += 1
					sum += sumarray[j]
				}}
				print chreffect[1], chreffect[2] , frequencies[i], sum, n, sum / length(sumarray)
				sum = 0
				n = 0
			}}
		}}' \
		{annotated_vcf} \
        > {os.path.dirname(ann)}/effectsummary.prog.csv
	
	mv {os.path.dirname(ann)}/effectsummary.prog.csv {outputs['csv']}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def name_cds_site_count(idx: str, target: AnonymousTarget) -> str:
	return f'cds_site_count_{os.path.basename(target.outputs['sites']).replace("-", "_")}'

def cds_site_count(bam_file: str, gtf_annotation_file: str, output_directory: str, sample_group: str, sample_name: str, min_coverage: int = 300, max_coverage: int = 600):
	"""
	Template: Count number of potential sites in CDS region within coverage threshold.
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'bam': bam_file,
           	  'gtf': gtf_annotation_file}
	outputs = {'bed': f'{output_directory}/snpEff/{sample_group}/{sample_name}/calculated_genetic_load/tmp/{os.path.basename(os.path.splitext(gtf_annotation_file)[0])}.cds.bed',
               'sites': f'{output_directory}/snpEff/{sample_group}/{sample_name}/calculated_genetic_load/tmp/{sample_name}.sitecount.tsv'}
	options = {
		'cores': 18,
		'memory': '30g',
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
	
	[ -d {output_directory}/snpEff/{sample_group}/{sample_name}/calculated_genetic_load/tmp ] || mkdir -p {output_directory}/snpEff/{sample_group}/{sample_name}/calculated_genetic_load/tmp
	
	awk \
        'BEGIN{{
			FS = OFS = "\\t"
		}}
        {{
			if ($3 == "CDS") {{
				if ($7 == "+") {{
					npos += 1
					print $1, $4 - 1, $5, "pos" npos, ".", $7
				}}
				if ($7 == "-") {{
					nneg += 1
					print $1, $4 - 1, $5, "neg" nneg, ".", $7
				}}
			}}
		}}' \
		{gtf_annotation_file} \
		> {output_directory}/snpEff/{sample_group}/{sample_name}/calculated_genetic_load/tmp/{os.path.basename(os.path.splitext(gtf_annotation_file)[0])}.cds.bed

    samtools depth \
		-@ {options['cores'] - 1} \
        -b {output_directory}/snpEff/{sample_group}/{sample_name}/calculated_genetic_load/tmp/{os.path.basename(os.path.splitext(gtf_annotation_file)[0])}.cds.bed \
        {bam_file} \
    | awk \
		'BEGIN{{
			FS = OFS = "\\t"
		}}
		{{
            if ($3 > {min_coverage} && $3 < {max_coverage}) {{
				chromosomearray[$1] += 1
			}}
		}}
        END{{
			print "chromosome", "n_positions"
			for (i in chromosomearray) {{
				print i, chromosomearray[i]
			}}
		}}' \
		> {output_directory}/snpEff/{sample_group}/{sample_name}/calculated_genetic_load/tmp/{sample_name}.sitecount.prog.tsv

	mv {output_directory}/snpEff/{sample_group}/{sample_name}/calculated_genetic_load/tmp/{sample_name}.sitecount.prog.tsv {outputs['sites']}
    
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def snpeff_result(effectsummary_file: str, sitecount_file: str, output_directory: str, sample_group: str, sample_name: str, snpeff_calcresult: str = f'{os.path.dirname(os.path.realpath(__file__))}/software/snpeff_calcresult.py'):
	"""
	Template: template_description
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'effectsummary': effectsummary_file,
           	  'sitecount': sitecount_file}
	outputs = {'tsv': f'{output_directory}/snpEff/{sample_group}/{sample_name}/calculated_genetic_load/snpeff_results.tsv'}
	options = {
		'cores': 1,
		'memory': '18g',
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
	
	[ -d {output_directory}/snpEff/{sample_group}/{sample_name}/calculated_genetic_load ] || mkdir -p {output_directory}/snpEff/{sample_group}/{sample_name}/calculated_genetic_load
	
	python {snpeff_calcresult} \
		{sample_name} \
		{sample_group} \
		{effectsummary_file} \
		{sitecount_file} \
        {output_directory}/snpEff/{sample_group}/{sample_name}/calculated_genetic_load/snpeff_results.prog.tsv
	
	mv {output_directory}/snpEff/{sample_group}/{sample_name}/calculated_genetic_load/snpeff_results.prog.tsv {outputs['tsv']}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def concatenate_snpeff_results(files: list, output_name: str, output_directory: str):
	"""
	Template: template_description
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'files': files}
	outputs = {'concat_file': f'{output_directory}/{output_name}.tsv'}
	options = {
		'cores': 1,
		'memory': '10g',
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
	
	[ -d {output_directory} ] || mkdir -p {output_directory}
	
	awk \
		'BEGIN{{FS=OFS="\\t"}}
        {{
        	if (FNR == 1 && NR != 1)
				{{next}}
            else
				{{print}}
        }}' \
        {' '.join(files)} \
        > {output_directory}/{output_name}.prog.tsv
        
    mv {output_directory}/{output_name}.prog.tsv {outputs['concat_file']}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def concat(files: list, output_name: str, output_directory: str = None, compress: bool = False):
    """
    Template: Name-sorts and concatenates files. Optionally compresses output using :script:`gzip`.
    
    Template I/O::
    
        inputs = {'files': files}
        outputs = {'concat_file': output_name.ext | output_name.ext.gzip}
    
    :param list files:
        List containing files to concatenate.
    :param str output_name:
        Desired name of output file, no extension.
    :param str output_directory:
        Path to output directory. Default is directory of 'files[0]'.
    :param bool compress:
        Bool indicating whether the output file should be compressed or not.
    """
    if output_directory is None:
        output_directory = os.path.dirname(files[0])
    inputs = {'files': files}
    if compress:
        outputs = {'concat_file': f'{output_directory}/{output_name}{os.path.splitext(files[0])[1]}.gz'}
    else:
        outputs = {'concat_file': f'{output_directory}/{output_name}{os.path.splitext(files[0])[1]}'}
    options = {
        'cores': 2,
        'memory': '16g',
        'walltime': '24:00:00'
    }
    protect = outputs['concat_file']
    spec = f"""
    # Sources environment
    if [ "$USER" == "jepe" ]; then
        source /home/"$USER"/.bashrc
        source activate popgen
    fi
    
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    
    [ -d {output_directory}] || mkdir -p {output_directory}

    if [ {compress} == 'False' ]; then
        cat \
            {' '.join(files)} \
            > {output_directory}/{output_name}.prog{os.path.splitext(files[0])[1]}
        
        mv {output_directory}/{output_name}.prog{os.path.splitext(files[0])[1]} {outputs['concat_file']}
    else
        cat \
            {' '.join(files)} \
        | gzip \
            -c \
            - \
            > {output_directory}/{output_name}.prog{os.path.splitext(files[0])[1]}.gz
        
        mv {output_directory}/{output_name}.prog{os.path.splitext(files[0])[1]}.gz {outputs['concat_file']}
    fi

    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, protect=protect, spec=spec)

def name_snpgenie(idx: str, target: AnonymousTarget) -> str:
	return f'snpgenie_{os.path.basename(os.path.dirname(os.path.dirname(os.path.dirname(target.outputs['plus']['param'])))).replace("-", "_")}_{os.path.basename(os.path.dirname(os.path.dirname(target.outputs['plus']['param']))).replace("-", "_")}'

def snpgenie_withinpool(reference_genome_file: str, gtf_annotation_file: str, vcf_file: str, sample_group: str, sample_name: str, region: str, region_length: int,  output_directory: str, min_allele_frequency: int | float = 0, sliding_window_size: int = 9):
	"""
	Template: Estimate pi_N/pi_S for each chromosome, separately for 'x' and '-' strand
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	vcf = f'<(zcat {vcf_file})' if vcf_file.endswith('.gz') else vcf_file
	inputs = {'reference': reference_genome_file,
		   	  'gtf': gtf_annotation_file,
			  'vcf': vcf_file}
	outputs = {'plus': {'vcf': f'{output_directory}/snpgenie/{sample_group}/tmp/{sample_name}/{region}/{sample_group}.{sample_name}.{region}.vcf',
                        'reference': f'{output_directory}/snpgenie/{sample_group}/tmp/{sample_name}/{region}/{os.path.splitext(os.path.basename(reference_genome_file))[0]}.{region}.fasta',
						'gtf': f'{output_directory}/snpgenie/{sample_group}/tmp/{sample_name}/{region}/{os.path.splitext(os.path.basename(gtf_annotation_file))[0]}.{region}.gtf',
						'param': f'{output_directory}/snpgenie/{sample_group}/{sample_name}/{region}/plus/SNPGenie_parameters.txt',
						'log': f'{output_directory}/snpgenie/{sample_group}/{sample_name}/{region}/plus/SNPGenie_LOG.txt',
						'site': f'{output_directory}/snpgenie/{sample_group}/{sample_name}/{region}/plus/site_results.txt',
						'codon': f'{output_directory}/snpgenie/{sample_group}/{sample_name}/{region}/plus/codon_results.txt',
						'product': f'{output_directory}/snpgenie/{sample_group}/{sample_name}/{region}/plus/product_results.txt',
						'summary': f'{output_directory}/snpgenie/{sample_group}/{sample_name}/{region}/plus/population_summary.txt',
						'window': f'{output_directory}/snpgenie/{sample_group}/{sample_name}/{region}/plus/sliding_window_length{sliding_window_size}_results.txt'},
               'minus': {'vcf': f'{output_directory}/snpgenie/{sample_group}/tmp/{sample_name}/{region}/{sample_group}.{sample_name}.{region}_revcom.vcf',
                         'reference': f'{output_directory}/snpgenie/{sample_group}/tmp/{sample_name}/{region}/{os.path.splitext(os.path.basename(reference_genome_file))[0]}.{region}_revcom.fasta',
						 'gtf': f'{output_directory}/snpgenie/{sample_group}/tmp/{sample_name}/{region}/{os.path.splitext(os.path.basename(gtf_annotation_file))[0]}.{region}_revcom.gtf',
						 'param': f'{output_directory}/snpgenie/{sample_group}/{sample_name}/{region}/minus/SNPGenie_parameters.txt',
						 'log': f'{output_directory}/snpgenie/{sample_group}/{sample_name}/{region}/minus/SNPGenie_LOG.txt',
						 'site': f'{output_directory}/snpgenie/{sample_group}/{sample_name}/{region}/minus/site_results.txt',
						 'codon': f'{output_directory}/snpgenie/{sample_group}/{sample_name}/{region}/minus/codon_results.txt',
						 'product': f'{output_directory}/snpgenie/{sample_group}/{sample_name}/{region}/minus/product_results.txt',
						 'summary': f'{output_directory}/snpgenie/{sample_group}/{sample_name}/{region}/minus/population_summary.txt',
						 'window': f'{output_directory}/snpgenie/{sample_group}/{sample_name}/{region}/minus/sliding_window_length{sliding_window_size}_results.txt'}}
	options = {
		'cores': 18,
		'memory': '30g',
		'walltime': '48:00:00'
	}
	spec = f"""
	# Sources environment
	if [ "$USER" == "jepe" ]; then
		source /home/"$USER"/.bashrc
		source activate snpgenie
	fi
	
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	
	[ -d {output_directory}/snpgenie/{sample_group}/{sample_name}/{region}/plus ] && rm -rf {output_directory}/snpgenie/{sample_group}/{sample_name}/{region}/plus
    [ -d {output_directory}/snpgenie/{sample_group}/{sample_name}/{region}/minus ] && rm -rf {output_directory}/snpgenie/{sample_group}/{sample_name}/{region}/minus
    [ -d {output_directory}/snpgenie/{sample_group}/{sample_name}/{region} ] || mkdir -p {output_directory}/snpgenie/{sample_group}/{sample_name}/{region}
	[ -d {output_directory}/snpgenie/{sample_group}/tmp/{sample_name}/{region} ] || mkdir -p {output_directory}/snpgenie/{sample_group}/tmp/{sample_name}/{region}
    
    cd {output_directory}/snpgenie/{sample_group}/tmp/{sample_name}/{region}
      
    awk -v region={region} \
        'BEGIN{{OFS=FS="\\t"}}
        {{if ($0 ~ /^#/ || $1 == region)
            {{print}}
        }}' \
        {vcf} \
        > {output_directory}/snpgenie/{sample_group}/tmp/{sample_name}/{region}/{sample_group}.{sample_name}.{region}.vcf
    
    [ -e {output_directory}/snpgenie/{sample_group}/tmp/{sample_name}/{region}/{sample_group}.{sample_name}.{region}_revcom.vcf ] && rm -f {output_directory}/snpgenie/{sample_group}/tmp/{sample_name}/{region}/{sample_group}.{sample_name}.{region}_revcom.vcf
    
    vcf2revcom.pl \
		{output_directory}/snpgenie/{sample_group}/tmp/{sample_name}/{region}/{sample_group}.{sample_name}.{region}.vcf \
        {region_length}
      
    awk -v region={region} \
		'BEGIN{{RS=">"; ORS=""; FS=OFS="\\n"}}
        {{if (NR > 1 && $1 == region)
            {{print ">"$0}}
        }}' \
		{reference_genome_file} \
		> {output_directory}/snpgenie/{sample_group}/tmp/{sample_name}/{region}/{os.path.splitext(os.path.basename(reference_genome_file))[0]}.{region}.fasta
    
	[ -e {output_directory}/snpgenie/{sample_group}/tmp/{sample_name}/{region}/{os.path.splitext(os.path.basename(reference_genome_file))[0]}.{region}_revcom.fasta ] && rm -f {output_directory}/snpgenie/{sample_group}/tmp/{sample_name}/{region}/{os.path.splitext(os.path.basename(reference_genome_file))[0]}.{region}_revcom.fasta

    fasta2revcom.pl \
		{output_directory}/snpgenie/{sample_group}/tmp/{sample_name}/{region}/{os.path.splitext(os.path.basename(reference_genome_file))[0]}.{region}.fasta

    awk -v region={region} \
        'BEGIN{{OFS=FS="\\t"}}
        {{if ($1 == region)
            {{print}}
        }}' \
        {gtf_annotation_file} \
        > {output_directory}/snpgenie/{sample_group}/tmp/{sample_name}/{region}/{os.path.splitext(os.path.basename(gtf_annotation_file))[0]}.{region}.gtf
    
	[ -e {output_directory}/snpgenie/{sample_group}/tmp/{sample_name}/{region}/{os.path.splitext(os.path.basename(gtf_annotation_file))[0]}.{region}_revcom.gtf ] && rm -f {output_directory}/snpgenie/{sample_group}/tmp/{sample_name}/{region}/{os.path.splitext(os.path.basename(gtf_annotation_file))[0]}.{region}_revcom.gtf
    
    gtf2revcom.pl \
		{output_directory}/snpgenie/{sample_group}/tmp/{sample_name}/{region}/{os.path.splitext(os.path.basename(gtf_annotation_file))[0]}.{region}.gtf \
        {region_length}

    echo -e "#########################\\n# Processing '+' strand #\\n#########################"
    
	snpgenie.pl \
		--vcfformat=2 \
		--snpreport={output_directory}/snpgenie/{sample_group}/tmp/{sample_name}/{region}/{sample_group}.{sample_name}.{region}.vcf \
		--fastafile={output_directory}/snpgenie/{sample_group}/tmp/{sample_name}/{region}/{os.path.splitext(os.path.basename(reference_genome_file))[0]}.{region}.fasta \
		--gtffile={output_directory}/snpgenie/{sample_group}/tmp/{sample_name}/{region}/{os.path.splitext(os.path.basename(gtf_annotation_file))[0]}.{region}.gtf \
		--workdir={output_directory}/snpgenie/{sample_group}/tmp/{sample_name}/{region} \
		--outdir={output_directory}/snpgenie/{sample_group}/{sample_name}/{region}/plus \
		--minfreq={min_allele_frequency} \
		--slidingwindow={sliding_window_size}
	
    echo -e "#########################\\n# Processing '-' strand #\\n#########################"

    snpgenie.pl \
		--vcfformat=2 \
		--snpreport={output_directory}/snpgenie/{sample_group}/tmp/{sample_name}/{region}/{sample_group}.{sample_name}.{region}_revcom.vcf \
		--fastafile={output_directory}/snpgenie/{sample_group}/tmp/{sample_name}/{region}/{os.path.splitext(os.path.basename(reference_genome_file))[0]}.{region}_revcom.fasta \
		--gtffile={output_directory}/snpgenie/{sample_group}/tmp/{sample_name}/{region}/{os.path.splitext(os.path.basename(gtf_annotation_file))[0]}.{region}_revcom.gtf \
		--workdir={output_directory}/snpgenie/{sample_group}/tmp/{sample_name}/{region}/ \
		--outdir={output_directory}/snpgenie/{sample_group}/{sample_name}/{region}/minus \
		--minfreq={min_allele_frequency} \
		--slidingwindow={sliding_window_size}

	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def snpgenie_summarize_results_population(population_summary_files: list, output_directory: str, species_name: str):
	"""
	Template: template_description
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'summaries': population_summary_files}
	outputs = {'tsv': f'{output_directory}/snpgenie/{species_abbreviation(species_name)}.snpgenie_results.tsv'}
	options = {
		'cores': 1,
		'memory': '18g',
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
	
	[ -d {output_directory}/snpgenie ] || mkdir -p {output_directory}/snpgenie
	
	awk \
		'BEGIN{{
			FS=OFS="\\t"
            print "sample", "group", "chromosome", "strand", "piN", "piS", "piN/piS"
		}}
        {{
			if (FNR == 1) {{
				next
			}}
            else {{
				split(FILENAME, filenamearray, "/")
                group=filenamearray[13]
                sample=filenamearray[14]
                chromosome=filenamearray[15]
                strand=filenamearray[16]
                piN=$10
                piS=$11
                print sample, group, chromosome, strand, piN, piS, piN_piS
			}}
		}}' \
        {' '.join(population_summary_files)} \
        > {output_directory}/snpgenie/{species_abbreviation(species_name)}.snpgenie_results.prog.tsv
	
	mv {output_directory}/snpgenie/{species_abbreviation(species_name)}.snpgenie_results.prog.tsv {outputs['tsv']}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)