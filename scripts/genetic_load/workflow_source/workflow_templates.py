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

def snpeff_freqs(ann: str, csv: str, txt: str, html: str):
	"""
	Template: Summarizes annotated variant functions.
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'vcf': ann}
	outputs = {'csv': f'{os.path.dirname(ann)}/effectsummary.tsv'}
	options = {
		'cores': 1,
		'memory': '12g',
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
	
	awk \
        'BEGIN{{FS=OFS="\\t"}}
        {{if ($0 !~ /^#/) 
            {{print $8}}
        }}' \
        {ann} \
    | awk \
        'BEGIN{{FS=";"}}
        {{
        id = 0;
        for (i = 1; i <= NF; i++)
            {{
            if ($i ~ /AF=/)
				{{
                split($i, af, "=");
                frq = af[2]
                }}
            if ($i ~ /\\|LOW\\|/)
                {{
                id = 1
				}}
            if ($i ~ /\\|MODERATE\\|/)
                {{
                id = 2
				}}
            if ($i ~ /\\|HIGH\\|/)
                {{
                id = 3
				}}
            if (length(effectarray[id]) < 0)
				{{
                effectarray[id] = frq
				}}
            else
            	{{
				effectarray[id] = effectarray[id] "," frq
				}}
            id = 0;
            }}
        }}
        END{{
        for (j = 1; j <= 3; j++)
			{{
            print effectarray[j]
			}}
		}}' \
        > {os.path.dirname(ann)}/effectsummary.prog.csv
	
	mv {os.path.dirname(ann)}/effectsummary.prog.csv {outputs['csv']}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def name_snpgenie(idx: str, target: AnonymousTarget) -> str:
	return f'snpgenie_{os.path.basename(os.path.dirname(os.path.dirname(target.outputs['param']))).replace("-", "_")}_{os.path.basename(os.path.dirname(target.outputs['param'])).replace("-", "_")}'

def snpgenie_withinpool(reference_genome_file: str, gtf_annotation_file: str, vcf_file: str, sample_group: str, sample_name: str, region: str,  output_directory: str, min_allele_frequency: int | float = 0, sliding_window_size: int = 9):
	"""
	Template: Estimate pi_N/pi_S
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'reference': reference_genome_file,
		   	  'gtf': gtf_annotation_file,
			  'vcf': vcf_file}
	outputs = {'vcf': f'{output_directory}/snpgenie/{sample_group}/tmp/{sample_name}/{region}/{sample_group}.{sample_name}.{region}.vcf',
               'reference': f'{output_directory}/snpgenie/{sample_group}/tmp/{sample_name}/{region}/{os.path.splitext(os.path.basename(reference_genome_file))[0]}.{region}.fasta',
               'gtf': f'{output_directory}/snpgenie/{sample_group}/tmp/{sample_name}/{region}/{sample_group}.{sample_name}.{region}.gtf',
               'param': f'{output_directory}/snpgenie/{sample_group}/{sample_name}/{region}/SNPGenie_parameters.txt',
			   'log': f'{output_directory}/snpgenie/{sample_group}/{sample_name}/{region}/SNPGenie_LOG.txt',
			   'site': f'{output_directory}/snpgenie/{sample_group}/{sample_name}/{region}/site_results.txt',
			   'codon': f'{output_directory}/snpgenie/{sample_group}/{sample_name}/{region}/codon_results.txt',
			   'product': f'{output_directory}/snpgenie/{sample_group}/{sample_name}/{region}/product_results.txt',
			   'summary': f'{output_directory}/snpgenie/{sample_group}/{sample_name}/{region}/population_summary.txt',
			   'window': f'{output_directory}/snpgenie/{sample_group}/{sample_name}/{region}/sliding_window_length_{sliding_window_size}_results.txt'}
	options = {
		'cores': 18,
		'memory': '30g',
		'walltime': '36:00:00'
	}
	spec = f"""
	# Sources environment
	if [ "$USER" == "jepe" ]; then
		source /home/"$USER"/.bashrc
		source activate snpgenie
	fi
	
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	
	[ -d {output_directory}/snpgenie/{sample_group}/{sample_name}/{region} ] && rm -rf {output_directory}/snpgenie/{sample_group}/{sample_name}/{region}
    [ -d {output_directory}/snpgenie/{sample_group}/{sample_name} ] || mkdir -p {output_directory}/snpgenie/{sample_group}/{sample_name}
	[ -d {output_directory}/snpgenie/{sample_group}/tmp/{sample_name}/{region} ] || mkdir -p {output_directory}/snpgenie/{sample_group}/tmp/{sample_name}/{region}
    
    awk -v region={region} \
        'BEGIN{{OFS=FS="\\t"}}
        {{if ($0 ~ /^#/ || $1 == region)
            {{print}}
        }}' \
        {vcf_file} \
        > {output_directory}/snpgenie/{sample_group}/tmp/{sample_name}/{region}/{sample_group}.{sample_name}.{region}.vcf
      
    awk -v region={region} \
		'BEGIN{{RS=">"; ORS=""; FS=OFS="\\t"}}
        {{if (NR > 1 && $1 == region)
            {{print ">"$0}}
        }}' \
		{reference_genome_file} \
		> {output_directory}/snpgenie/{sample_group}/tmp/{sample_name}/{region}/{os.path.splitext(os.path.basename(reference_genome_file))[0]}.{region}.fasta
    
    awk -v region={region} \
        'BEGIN{{OFS=FS="\\t"}}
        {{if ($1 == region)
            {{print}}
        }}' \
        {gtf_annotation_file} \
        > {output_directory}/snpgenie/{sample_group}/tmp/{sample_name}/{region}/{sample_group}.{sample_name}.{region}.gtf
            
	snpgenie.pl \
		--vcfformat=4 \
		--snpreport={output_directory}/snpgenie/{sample_group}/tmp/{sample_name}/{region}/{sample_group}.{sample_name}.{region}.vcf \
		--fastafile={output_directory}/snpgenie/{sample_group}/tmp/{sample_name}/{region}/{os.path.splitext(os.path.basename(reference_genome_file))[0]}.{region}.fasta \
		--gtffile={output_directory}/snpgenie/{sample_group}/tmp/{sample_name}/{region}/{sample_group}.{sample_name}.{region}.gtf \
		--workdir={output_directory}/snpgenie/{sample_group}/tmp/{sample_name}/{region} \
		--outdir={output_directory}/snpgenie/{sample_group}/{sample_name}/{region} \
		--minfreq={min_allele_frequency} \
		--slidingwindow={sliding_window_size}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)