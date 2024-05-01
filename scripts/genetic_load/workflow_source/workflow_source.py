#!/bin/env python3
from gwf import Workflow
from gwf.workflow import collect
import os, yaml, glob, sys, gzip
from workflow_templates import *

def genetic_load_workflow(config_file: str = glob.glob('*config.y*ml')[0]):
	"""
	Workflow: description
	
	:param str config_file:
		Configuration file containing pre-defined set of variables
	"""
	# --------------------------------------------------
	#                  Configuration
	# --------------------------------------------------
	
	CONFIG = yaml.safe_load(open(config_file))
	ACCOUNT: str = CONFIG['account']
	SPECIES_NAME: str = CONFIG['species_name']
	REFERENCE: str = CONFIG['reference_genome_file']
	GTF: str = CONFIG['gtf_annotation_file']
	WORK_DIR: str = CONFIG['working_directory_path']
	OUTPUT_DIR: str = CONFIG['output_directory_path']
	SNPGENIE_SETTINGS: dict = CONFIG['snpgenie_settings']
	SNPGENIE_MINFREQ: float | int = SNPGENIE_SETTINGS['minfreq'] if SNPGENIE_SETTINGS['minfreq'] else 0
	SNPGENIE_SLIDINGWINDOW: int = SNPGENIE_SETTINGS['slidingwindow'] if SNPGENIE_SETTINGS['slidingwindow'] else 9
	SNPGGENIE_VCF_FORMAT: int = SNPGENIE_SETTINGS['vcf_format'] if SNPGENIE_SETTINGS['vcf_format'] else 2
	DOS_SETTINGS: dict = CONFIG['dos_settings']
	DOS_OUTGROUP_NAME: str = DOS_SETTINGS['outgroup_name'].replace(' ', '_') if DOS_SETTINGS['outgroup_name'] else None
	DOS_OUTGROUP_BAM: str = DOS_SETTINGS['outgroup_bam']
	DOS_OUTGROUP_VCF: str = DOS_SETTINGS['outgroup_vcf']
	SAMPLES: list = CONFIG['sample_list']
	
	SNPEFF_DIR = f'{os.path.dirname(os.path.realpath(__file__))}/software/snpeff'

	# --------------------------------------------------
	#                  Workflow
	# --------------------------------------------------
	
	gwf = Workflow(
		defaults={'account': ACCOUNT}
	)

	top_dir = f'{WORK_DIR}/{SPECIES_NAME.replace(" ", "_")}/genetic_load'

	if os.path.exists(f'reference_sequences_{os.path.splitext(os.path.basename(REFERENCE))[0]}.txt'):
		with open(f'reference_sequences_{os.path.splitext(os.path.basename(REFERENCE))[0]}.txt', 'r') as infile:
			sequences = [{'sequence_name': entry.split(sep='\t')[0].strip(), 'sequence_length': entry.split(sep='\t')[1].strip()} for entry in infile]
	else:
		sequences = parse_fasta(REFERENCE)
		with open(f'reference_sequences_{os.path.splitext(os.path.basename(REFERENCE))[0]}.txt', 'w') as outfile:
			outfile.write('\n'.join('\t'.join(str(i) for i in entry.values()) for entry in sequences))

	# Construct SnpEff database entry for reference genome.
	database_entry = gwf.target_from_template(
		name=f'{species_abbreviation(SPECIES_NAME)}_snpeff_database_entry',
		template=snpeff_database_build(
			gtf_annotation_file=GTF,
			reference_genome_file=REFERENCE,
			species_name=SPECIES_NAME,
			snpeff_directory=SNPEFF_DIR
		)
	)
	
	# Create bed file of gene regions.
	gene_bed = gwf.target_from_template(
		name=f'gene_bed_{species_abbreviation(SPECIES_NAME)}',
		template=gtf2gene_bed(
			gtf_annotation_file=GTF,
			output_directory=top_dir
		)
	)

	# Create a consensus sequence from outgroup alignment file to approximate ancestral state.
	# Job only available if outgroup alignment file has been supplied.
	if DOS_OUTGROUP_BAM:
		outgroup_consensus = gwf.target_from_template(
			name=f'consensus_sequence_outgroup',
			template=outgroup_consensus_sequence(
				outgroup_bam_file=DOS_OUTGROUP_BAM,
				output_directory=top_dir,
				outgroup_name=DOS_OUTGROUP_NAME
			)
		)

		outgroup_pileup = gwf.target_from_template(
			name=f'consensus_pileup_outgroup',
			template=outgroup_consensus_pileup(
				outgroup_bam_file=DOS_OUTGROUP_BAM,
				output_directory=top_dir,
				outgroup_name=DOS_OUTGROUP_NAME
			)
		)

	# Annotate outgroup VCF with SnpEff.
	# Job only available if outgroup VCF has been supplied.
	if DOS_OUTGROUP_VCF:
		outgroup_variant_annotation = gwf.target_from_template(
			name=f'snpeff_annotation_outgroup',
			template=snpeff_annotation_outgroup(
				vcf_file=DOS_OUTGROUP_VCF,
				snpeff_predictor_file=database_entry.outputs['predictor'],
				snpeff_config_file=f'{SNPEFF_DIR}/snpEff.config',
				output_directory=top_dir,
				outgroup_name=DOS_OUTGROUP_NAME
			)
		)

		ancestral_allele = gwf.target_from_template(
			name=f'ancestral_allele_information',
			template=ancestral_allele_information(
				outgroup_vcf_file=DOS_OUTGROUP_VCF,
				output_directory=top_dir,
				outgroup_name=DOS_OUTGROUP_NAME
			)
		)

	# Initialize lists to hold output files
	snpeff_results_list = []
	snpgenie_results_list = []
	dos_results_list = []

	# Iterates through all supplied samples.
	for sample in SAMPLES:
		SAMPLE_NAME = sample['sample_name']
		SAMPLE_GROUP = sample['sample_group']
		VCF = sample['vcf_file']
		BAM = sample['bam_file']

		# Annotate VCF file with SnpEff
		variant_annotation = gwf.target_from_template(
			name=f'snpeff_annotation_{SAMPLE_GROUP}_{SAMPLE_NAME.replace("-", "_")}',
			template=snpeff_annotation(
				vcf_file=VCF,
				snpeff_predictor_file=database_entry.outputs['predictor'],
				snpeff_config_file=f'{SNPEFF_DIR}/snpEff.config',
				output_directory=top_dir,
				sample_group=SAMPLE_GROUP,
				sample_name=SAMPLE_NAME
			)
		)

		# Calculate frequencies of each impact class from SnpEff annotation.
		annotation_frequencies = gwf.target_from_template(
			name=f'snpeff_frequencies_{SAMPLE_GROUP}_{SAMPLE_NAME.replace("-", "_")}',
			template=snpeff_freqs(
				ann=variant_annotation.outputs['ann']
			)
		)

		# Jobs only available if BAM alignment files have been supplied.
		if BAM:
			# Count number of potential sites in CDS regions
			site_count = gwf.target_from_template(
				name=f'cds_site_count_{SAMPLE_GROUP}_{SAMPLE_NAME.replace("-", "_")}',
				template=cds_site_count(
					bam_file=BAM,
					gtf_annotation_file=GTF,
					output_directory=top_dir,
					sample_group=SAMPLE_GROUP,
					sample_name=SAMPLE_NAME,
					min_coverage=300,
					max_coverage=600
				)
			)

			# Calculate genetic load based on SnpEff annotation.
			snpeff_results = gwf.target_from_template(
				name=f'snpeff_results_{SAMPLE_GROUP}_{SAMPLE_NAME.replace("-", "_")}',
				template=snpeff_result(
					effectsummary_file=annotation_frequencies.outputs['csv'],
					sitecount_file=site_count.outputs['sites'],
					output_directory=top_dir,
					sample_group=SAMPLE_GROUP,
					sample_name=SAMPLE_NAME)
			)

			snpeff_results_list.append(snpeff_results.outputs['tsv'])

		# Calculate pi values for each chromosome in each direction.
		for sequence in sequences:
			REGION = sequence['sequence_name']
			REGION_LENGTH = sequence['sequence_length']

			snpgenie_pi = gwf.target_from_template(
				name=f'snpgenie_{SAMPLE_GROUP}_{SAMPLE_NAME.replace("-", "_")}_{REGION}',
				template=snpgenie_withinpool(
					reference_genome_file=REFERENCE,
					gtf_annotation_file=GTF,
					vcf_file=VCF,
					sample_group=SAMPLE_GROUP,
					sample_name=SAMPLE_NAME,
					region=REGION,
					region_length=REGION_LENGTH,
					output_directory=top_dir,
					min_allele_frequency=SNPGENIE_MINFREQ,
					sliding_window_size=SNPGENIE_SLIDINGWINDOW,
					vcf_format=SNPGGENIE_VCF_FORMAT
				)
			)

			snpgenie_results_list.append(snpgenie_pi.outputs['plus']['summary'])
			snpgenie_results_list.append(snpgenie_pi.outputs['minus']['summary'])

		# Count number of synonymous and non-synonymous polymorphisms in gene regions.
		dos_polymorphic = gwf.target_from_template(
			name=f'dos_polymorphic_site_count_{SAMPLE_GROUP}_{SAMPLE_NAME.replace("-", "_")}',
			template=dos_count_polymorphic_sites(
				snpeff_annotated_vcf_file=variant_annotation.outputs['ann'],
				gene_bed_file=gene_bed.outputs['bed'],
				output_directory=top_dir,
				sample_group=SAMPLE_GROUP,
				sample_name=SAMPLE_NAME
			)
		)

		# Jobs only available if outgroup VCF has been supplied.
		if DOS_OUTGROUP_VCF:
			# Count number of synonymous and non-synonymous substitutions in gene regions.
			dos_substitution = gwf.target_from_template(
				name=f'dos_substitution_site_count_{SAMPLE_GROUP}_{SAMPLE_NAME.replace("-", "_")}',
				template=dos_count_substitution_sites(
					snpeff_annotated_vcf_file=outgroup_variant_annotation.outputs['ann'],
					gene_bed_file=gene_bed.outputs['bed'],
					polymorphism_variants_file=dos_polymorphic.outputs['var'],
					output_directory=top_dir,
					sample_group=SAMPLE_GROUP,
					sample_name=SAMPLE_NAME
				)
			)

			# Combine counts of polymorphisms and substitutions.
			dos_combine = gwf.target_from_template(
				name=f'dos_combine_{SAMPLE_GROUP}_{SAMPLE_NAME.replace("-", "_")}',
				template=dos_combine_d_p(
					polymorphisms_file=dos_polymorphic.outputs['poly'],
					substitutions_file=dos_substitution.outputs['sub'],
					output_directory=top_dir,
					sample_group=SAMPLE_GROUP,
					sample_name=SAMPLE_NAME
				)
			)

			# Calculate DoS.
			dos_result = gwf.target_from_template(
				name=f'dos_result_{SAMPLE_GROUP}_{SAMPLE_NAME.replace("-", "_")}',
				template=dos_results(
					combined_file=dos_combine.outputs['combined'],
					output_directory=top_dir,
					sample_group=SAMPLE_GROUP,
					sample_name=SAMPLE_NAME
				)
			)

			dos_results_list.append(dos_result.outputs['dos'])

			update_aa = gwf.target_from_template(
				name=f'update_ancestral_allele_{SAMPLE_GROUP}_{SAMPLE_NAME.replace("-", "_")}',
				template=update_ancestral_allele(
					ancestral_allele_file=ancestral_allele.outputs['aa'],
					vcf_file=VCF,
					output_directory=top_dir,
					sample_group=SAMPLE_GROUP,
					sample_name=SAMPLE_NAME
				)
			)

	if len(snpeff_results_list) >= 1:
		# Concatenates results from different populations into one file.
		concat_snpeff_results = gwf.target_from_template(
			name=f'snpeff_concatenate_{SPECIES_NAME.replace(" ", "_")}',
			template=snpeff_concatenate_results(
				files=snpeff_results_list,
				output_name=f'{species_abbreviation(SPECIES_NAME)}.snpeff_results',
				output_directory=f'{top_dir}/snpEff'
			)
		)

	# Concatenates results from different populations into one file.
	summarize_snpgenie = gwf.target_from_template(
		name=f'snpgenie_summarize_result_{SPECIES_NAME.replace(" ", "_")}',
		template=snpgenie_summarize_results_population(
			population_summary_files=snpgenie_results_list,
			output_directory=top_dir,
			species_name=SPECIES_NAME
		)
	)

	# Concatenates results from different populations into one file.
	if len(dos_results_list) >= 1:
		concat_dos_results = gwf.target_from_template(
			name=f'dos_concatenate_{SPECIES_NAME.replace(" ", "_")}',
			template=dos_concatenate_results(
				dos_results_files=dos_results_list,
				output_directory=top_dir,
				species_name=SPECIES_NAME
			)
		)

	return gwf