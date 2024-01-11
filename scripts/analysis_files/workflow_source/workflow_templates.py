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

def parse_fasta(fasta_file: str):
    """Parses `FASTA` file returning all sequence names and lengths paired in a list of dictionaries.
    
    ::
    
        return [{'sequence_name': str, 'sequence_length': int}, ...]
    
    :param str fasta_file:
        Sequence file in `FASTA` format.
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

def partition_chrom(parse_fasta: list, size: int = 500000):
    """
    Partitions `FASTA` file parsed with **parse_fasta**.
    
    Uses the list of dictionaries from **parse_fasta** to creates a list of dictionaries
    containing with partition number, sequence name, start and end postion (0 based).
    By default the partition size is 500kbs.
    
    ::
    
        return [{'num': int, 'region': str, 'start': int, 'end': int}, ...]
    
    :param list parse_fasta:
        List of dictionaries produced by **parse_fasta**.
    :param int size:
        Size of partitions. Default 500kb.
    """
    chrom_partition = []
    num = 1
    for chrom in parse_fasta:
        whole_chunks = chrom['sequence_length'] // size
        partial_chunk = chrom['sequence_length'] - whole_chunks * size
        start = 0
        for chunk in range(whole_chunks):
            end = start + size
            chrom_partition.append({'num': num, 'region': chrom['sequence_name'], 'start': start, 'end': end})
            start = end
            num += 1
        chrom_partition.append({'num': num, 'region': chrom['sequence_name'], 'start': start, 'end': start + partial_chunk})
        num += 1
    return chrom_partition

########################## PoolSNP ##########################

def mpileup_parts(bam_files: list, reference_genome: str, species_name: str, region: str, num: int, start: int, end: int, output_directory: str):
    """
    Template: Create :format:`mpileup` files for each partition of reference genome from multiple :format:`BAM` files using :script:`samtools mpileup`.
    
    Template I/O::
    
        inputs = {'bam_files': bam_files,
                  'reference': reference_genome}
        outputs = {'mpileup': *.mpileup}
    
    :param list bam_files:
        List of all :format:`BAM` files to be included in :format:`mpileup` file.
    :param str reference_genome:
        Path to genome reference file in `FASTA`format.
    :param str species_name:
        Name of species being worked on.
    :param str region:
        Name of chromosome from **partition_chrom**.
    :param int num:
        Partition number from **partition_chrom**.
    :param int start:
        Start position from **partition_chrom**.
    :param int end:
        End position from **partition_chrom**.
    :param str output_directory:
        Path to desired output directory. Creates directories 'tmp/bed' and 'tmp/mpileup' at location.
    """
    inputs = {'bam_files': bam_files,
              'reference': reference_genome}
    outputs = {'mpileup': f'{output_directory}/tmp/mpileup/{species_abbreviation(species_name)}_{num}_{region}.mpileup'}
    options = {
        'cores': 1,
        'memory': '16g',
        'walltime': '10:00:00'
    }
    spec = f"""
    # Sources environment
    if [ "$USER" == "jepe" ]; then
        source /home/"$USER"/.bashrc
        source activate vcf
    fi
    
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    
    sleep 2m

    [ -d {output_directory}/tmp/bed ] || mkdir -p {output_directory}/tmp/bed
    [ -d {output_directory}/tmp/mpileup ] || mkdir -p {output_directory}/tmp/mpileup

    echo -e '{region}\t{start}\t{end}' > {output_directory}/tmp/bed/{num}.bed
    
    samtools mpileup \
        --max-depth 0 \
        --fasta-ref {reference_genome} \
        --positions {output_directory}/tmp/bed/{num}.bed \
        --min-BQ 0 \
        --region {region} \
        --output {output_directory}/tmp/mpileup/{species_abbreviation(species_name)}_{num}_{region}.prog.mpileup \
        {' '.join(bam_files)}
    
    mv {output_directory}/tmp/mpileup/{species_abbreviation(species_name)}_{num}_{region}.prog.mpileup {outputs['mpileup']}
    
    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def mpileup2sync(mpileup_file: str, output_directory: str, mpileup2sync: str = glob.glob(f'{os.path.dirname(os.path.realpath(__file__))}/scripts/popoolation2*/mpileup2sync.pl')[0]):
    """
    Template: Makes a :format:`sync` file for a corresponding :format:`mpileup` file using :script:`popoolation2`'s :script:`mpileup2sync.pl`
    
    Template I/O::
    
        inputs = {'mpileup': mpileup_file}
        outputs = {'sync': *.sync}
    
    :param str mpileup_file:
        Input :format:`mpileup` file.
    :param str output_directory:
        Desired output directory for :format:`sync` file. Creates directories 'tmp/sync/' at location.
    :param str mpileup2sync:
        Path to :script:`mpile2sync.pl`.
    """
    inputs = {'mpileup': mpileup_file}
    outputs = {'sync': f'{output_directory}/{os.path.splitext(os.path.basename(mpileup_file))[0]}.sync',
               'params': f'{output_directory}/{os.path.splitext(os.path.basename(mpileup_file))[0]}.sync.params'}
    options = {
        'cores': 2,
        'memory': '16g',
        'walltime': '06:00:00'
    }
    spec = f"""
    # Sources environment
    if [ "$USER" == "jepe" ]; then
        source /home/"$USER"/.bashrc
        source activate vcf
    fi
    
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    
    [ -d {output_directory}/tmp/sync ] || mkdir -p {output_directory}/tmp/sync

    perl {mpileup2sync} \
        --input {mpileup_file} \
        --output {output_directory}/{os.path.splitext(os.path.basename(mpileup_file))[0]}.prog.sync \
        --fastq-type sanger
    
    mv {output_directory}/{os.path.splitext(os.path.basename(mpileup_file))[0]}.prog.sync {outputs['sync']}
    mv {output_directory}/{os.path.splitext(os.path.basename(mpileup_file))[0]}.prog.sync.params {outputs['params']}

    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def max_cov(mpileup: str, contig: str, cutoff: float, output_directory: str, script: str = f'{os.path.dirname(os.path.realpath(__file__))}/PoolSNP/scripts/max-cov.py'):
    """
    Template: Calculates coverage thresholds using :script:`max-cov.py`.
    
    Template I/O::
    
        inputs = {}
        outputs = {}
    
    :param
    """
    file_name = '{output_directory}/tmp/cov/cutoffs/{contig}'.format(output_directory=output_directory, contig=contig)
    inputs = {'mpileup': mpileup}
    outputs = {'cutoff': '{}.txt'.format(file_name)}
    options = {
        'cores': 1,
        'memory': '10g',
        'walltime': '96:00:00'
    }
    spec = """
    # Sources environment
    if [ "$USER" == "jepe" ]; then
        source /home/"$USER"/.bashrc
        source activate vcf
    fi
    
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    
    [ -d {output_directory}/tmp/cov/cutoffs ] || mkdir -p {output_directory}/tmp/cov/cutoffs

    presence=$(awk -v contig={contig} 'BEGIN{{presence = "no"}} {{if ($1 == contig) {{presence = "yes"; exit}} }} END{{print presence}}' {mpileup})

    if [ "$presence" == "yes" ]; then
        awk \
            -v contig={contig} \
            '{{if ($1 == contig) {{print $0}}}}' \
            {mpileup} \
        | python {script} \
            --mpileup - \
            --cutoff {cutoff} \
            --contig {contig} \
            --out {file_name}.prog.txt
        
        mv {file_name}.prog.txt {cutoff_file}
    else
        echo -n "" > {cutoff_file}
    fi
    
    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """.format(output_directory=output_directory, contig=contig, mpileup=mpileup, script=script, cutoff=cutoff, file_name=file_name, cutoff_file=outputs['cutoff'])
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
    files.sort()
    if output_directory is None:
        output_directory = os.path.dirname(files[0])
    inputs = {'files': files}
    if compress:
        outputs = {'concat_file': f'{output_directory}/{output_name}{os.path.splitext(files[0])[1]}.gzip'}
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
        source activate vcf
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
            > {output_directory}/{output_name}.prog{os.path.splitext(files[0])[1]}.gzip
        
        mv {output_directory}/{output_name}.prog{os.path.splitext(files[0])[1]}.gzip {outputs['concat_file']}
    fi

    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, protect=protect, spec=spec)