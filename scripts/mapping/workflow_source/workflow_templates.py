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

def get_sample_data(path: str) -> dict:
	"""Create dictionary of lists each containing all resequencing for for a sample in one read direction.
	
	:param str path:
		Path to sample directory containing resequencing data."""
	r1 = ('_1', '_R1')
	r2 = ('_2', '_R2')
	extformat = ('.fq', '.fa', '.fasta', 'fastq', '.fn')
	compress = ('', '.gz', '.gzip')
	reseq = {'sample_name': os.path.basename(path),
		  	 'r1': [os.path.join(path, file) for file in os.listdir(path) if file.endswith(tuple([i+j+k for i in r1 for j in extformat for k in compress]))],
		  	 'r2': [os.path.join(path, file) for file in os.listdir(path) if file.endswith(tuple([i+j+k for i in r2 for j in extformat for k in compress]))]}
	return reseq

########################## Mapping ##########################

def adapterremoval(sample_name: str, read1_files: list, read2_files: list, output_directory: str, min_qulaity: int = 25, min_length: int = 20, adapter1: str = 'AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA', adapter2: str = 'AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG'):
	"""
	Template: Remove remnant adapter sequences from read data using :script:`AdapterRemoval`.
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'read1': read1_files,
		   	  'read2': read2_files}
	outputs = {'pair1': f'{output_directory}/adapterremoval/{sample_name}/{sample_name}.pair1.truncated',
			   'pair2':f'{output_directory}/adapterremoval/{sample_name}/{sample_name}.pair2.truncated',
			   'collapsed': [f'{output_directory}/adapterremoval/{sample_name}/{sample_name}.collapsed',
			   				 f'{output_directory}/adapterremoval/{sample_name}/{sample_name}.collapsed.truncated'],
			   'misc': [f'{output_directory}/adapterremoval/{sample_name}/{sample_name}.settings',
			   			f'{output_directory}/adapterremoval/{sample_name}/{sample_name}.singleton.truncated',
			   			f'{output_directory}/adapterremoval/{sample_name}/{sample_name}.discarded']}
	options = {
		'cores': 16,
		'memory': '60g',
		'walltime': '06:00:00'
	}
	spec = f"""
	# Sources environment
	if [ "$USER" == "jepe" ]; then
		source /home/"$USER"/.bashrc
		source activate mapping
	fi
	
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	
	[ -d {output_directory}/adapterremoval/{sample_name} ] || mkdir -p {output_directory}/adapterremoval/{sample_name}
	
	AdapterRemoval \
		--threads {options['cores']} \
		--file1 {' '.join(read1_files)} \
		--file2 {' '.join(read2_files)} \
		--adapter1 {adapter1} \
		--adapter2 {adapter2} \
		--minquality {min_qulaity} \
		--minlength {min_length} \
		--basename {output_directory}/adapterremoval/{sample_name}.prog \
		--trimns \
		--trimqualities \
		--collapse
	
	mv {output_directory}/adapterremoval/{sample_name}/{sample_name}.prog.settings {outputs['misc'][0]}
	mv {output_directory}/adapterremoval/{sample_name}/{sample_name}.prog.pair1.truncated {outputs['pair1']}
	mv {output_directory}/adapterremoval/{sample_name}/{sample_name}.prog.pair2.truncated {outputs['pair2']}
	mv {output_directory}/adapterremoval/{sample_name}/{sample_name}.prog.singleton.truncated {outputs['misc'][1]}
	mv {output_directory}/adapterremoval/{sample_name}/{sample_name}.prog.collapsed {outputs['collapsed'][0]}
	mv {output_directory}/adapterremoval/{sample_name}/{sample_name}.prog.collapsed.truncated {outputs['collapsed'][1]}
	mv {output_directory}/adapterremoval/{sample_name}/{sample_name}.prog.discarded {outputs['misc'][2]}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)