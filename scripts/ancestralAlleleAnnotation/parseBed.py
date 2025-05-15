#!/bin/env python3
import sys

def parse_bed(bedFile: str):
	"""
	Parses :format:`BED` file returning a list of dictionaries
	containing chromsome name, number of sites within chromosome, and
	the start and end lines of the chromsome within the :format:`BED` file.

	::

		return [{'chrom': str, 'nSites': int, 'start': int, 'end': int}, ...]
	
	:param str bedFile:
		:format:`BED` file.
	"""
	bedList = []
	chromName = None
	siteCount = 0
	with open(bedFile, 'r') as bed:
		for lineNo, entry in enumerate(bed, start = 1):
			entry = entry.strip()
			entry = entry.split('\t')
			if chromName and chromName != entry[0]:
				end = lineNo - 1
				bedList.append({'chromName': chromName, 'nSites': siteCount, 'start': start, 'end': end})
				start = lineNo
				siteCount = 0
				chromName = entry[0]
			if not chromName:
				chromName = entry[0]
				start = lineNo
			siteCount += int(entry[2]) - int(entry[1])
		bedList.append({'chromName': chromName, 'nSites': siteCount, 'start': start, 'end': lineNo})
	return bedList

def partition_bed(parseBed: list, nLines: int = 250000):
	"""
	Takes the output from **parse_bed()** and splits each chromosome into chunks of nLines.

	::

		return [{'num': int, 'chromName': str, 'start' int, 'end': int}]

	:param list parseBed:
		Output from **parse_bed()** function. A list of dictionaries.
	:param int nLines:
		Number of lines to split each chromosome into.
	"""
	padding = 1
	for chrom in parseBed:
		wholeChunks = (chrom['end'] - chrom['start'] + 1) // nLines
		padding += (wholeChunks + 1)
	nPad = len(str(padding))
	bedPartition = []
	num = 1
	for chrom in parseBed:
		wholeChunks = (chrom['end'] - chrom['start'] + 1) // nLines
		partialChunk = (chrom['end'] - chrom['start'] + 1) - wholeChunks * nLines
		start = chrom['start']
		for chunk in range(wholeChunks):
			end = start + nLines - 1
			bedPartition.append({'num': f'{num:0{nPad}}', 'chromName': chrom['chromName'], 'start': start, 'end': end})
			start = end + 1
			num += 1
		if partialChunk:
			bedPartition.append({'num': f'{num:0{nPad}}', 'chromName': chrom['chromName'], 'start': start, 'end': start + partialChunk - 1})
			num += 1
	return bedPartition

print(parse_bed(sys.argv[1]))

print(partition_bed(parse_bed(sys.argv[1])))