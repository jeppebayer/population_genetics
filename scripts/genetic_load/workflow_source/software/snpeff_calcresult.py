#!/bin/env python3
import os, sys

sample_name = sys.argv[1]
sample_group = sys.argv[2]

with open(sys.argv[3], 'r') as effectsummary_file:
	next(effectsummary_file)
	effectsummary_table = [{'chr': row.split('\t')[0].rstrip(), 'effect': row.split('\t')[1].rstrip(), 'sum_frq': float(row.split('\t')[3].rstrip())} for row in effectsummary_file]

with open(sys.argv[4], 'r') as sitecount_file:
	next(sitecount_file)
	sitecount_table = [{'chr': row.split('\t')[0].rstrip(), 'n_positions': float(row.split('\t')[1].rstrip())} for row in sitecount_file]
	
result_table = [{'chr': i['chr'], 'effect': i['effect'], 'sum_frq': i['sum_frq'], 'n_positions': j['n_positions'], 'frequency': i['sum_frq'] / j['n_positions']} for i in effectsummary_table for j in sitecount_table if i['chr'] == j['chr']]

sorted_result_table = sorted(result_table, key=lambda d: (d['chr'], d['effect']))

with open(sys.argv[5], 'w') as outfile:
	outfile.write('sample\tgroup\tchromosome\teffect\tsum_frequencies\tnum_positions\teffect_frequency\n')
	for entry in sorted_result_table:
		outfile.write(f'{sample_name}\t{sample_group}\t{'\t'.join(str(k) for k in entry.values())}\n')