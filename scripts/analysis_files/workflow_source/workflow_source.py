#!/bin/env python3
from gwf import Workflow
from gwf.workflow import collect
import os, yaml, glob, sys
from workflow_templates import *

def poolsnp_vcf_workflow(config_file: str = glob.glob('*config.y*ml')[0]):
    """
    Workflow: Create :format:`mpileup`, :format:`sync` and :foramt:`VCF` files using :script:`samtools`, :script:`popoolation2` and :script:`PoolSNP`.
    
    :param str config_file:
        Configuration file containing pre-defined set of variables
    """
    # --------------------------------------------------
    #                  Configuration
    # --------------------------------------------------
    
    config = yaml.safe_load(open(config_file))
    ACCOUNT: str = config['account']
    SPECIES_NAME: str = config['species_name']
    SAMPLE_LIST: list = config['sample_list']
    sample_names: list = [os.path.basename(os.path.dirname(path)) for path in SAMPLE_LIST]
    REFERENCE_GENOME: str = config['reference_genome_path']
    WORKING_DIR: str = config['working_directory_path']
    OUTPUT_DIR: str = config['output_directory_path']
    MAXCOV: float = config['max_cov']
    MINCOV: int = config['min_cov']
    MINCOUNT: int = config['min_count']
    MINFREQ: float = config['min_freq']
    MISSFRAC: float = config['miss_frac']
    BASEQUAL: int = config['bq']
    ALLSITES: bool = config['all_sites']
    if ALLSITES == 'True':
        ALLSITES = 1,
    elif ALLSITES == 'False':
        ALLSITES = 0
    partition_size = 200000

    # --------------------------------------------------
    #                  Workflow
    # --------------------------------------------------
    
    gwf = Workflow(
        defaults={'account': ACCOUNT}
    )
    
    if os.path.exists(f'reference_partitions.{partition_size}bp.txt') and os.path.exists('reference_sequences.txt'):
        # If files exists reads data directly from files
        # Loads reference genome partitioning
        with open(f'reference_partitions.{partition_size}bp.txt', 'r') as infile:
            partitions = [{'num': entry.split(sep='\t')[0].strip(), 'region': entry.split(sep='\t')[1].strip(), 'start':entry.split(sep='\t')[2].strip(), 'end': entry.split(sep='\t')[3].strip()} for entry in infile]
        # Loads list of contigs in reference genome
        with open('reference_sequences.txt', 'r') as infile:
            contigs = [{'contig': entry.split(sep='\t')[0].strip()} for entry in infile]
    else:
        # If files don't exist, generate data and write files
        sequences = parse_fasta(REFERENCE_GENOME)
        with open(f'reference_sequences.txt', 'w') as outfile:
            outfile.write('\n'.join('\t'.join(str(i) for i in entry.values()) for entry in sequences))
        # Partitions reference genome
        partitions = partition_chrom(parse_fasta=sequences, size=partition_size)
        with open(f'reference_partitions.{partition_size}bp.txt', 'w') as outfile:
            outfile.write('\n'.join('\t'.join(str(i) for i in entry.values()) for entry in partitions))
        # Creates list of contigs in reference genome
        contigs = [{'contig': contig['sequence_name']} for contig in sequences]
    
    top_dir = f'{WORKING_DIR}/analysis_files/{SPECIES_NAME.replace(" ", "_")}'
    output_dir = f'{OUTPUT_DIR}/{species_abbreviation(SPECIES_NAME)}'

    mpileup = gwf.map(
        name=name_mpileup,
        template_func=mpileup_parts,
        inputs=partitions,
        extra={'bam_files': SAMPLE_LIST,
               'reference_genome': REFERENCE_GENOME,
               'species_name': SPECIES_NAME,
               'output_directory': top_dir}
    )

    return gwf