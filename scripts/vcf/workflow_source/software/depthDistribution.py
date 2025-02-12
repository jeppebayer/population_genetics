#!/bin/env python3
import sys, os.path, matplotlib.pyplot as plt, pandas as pd, numpy as np
from statsmodels.stats.weightstats import DescrStatsW

# Takes in positional arguments
# For mode:
# '0' uses static lower threshold and dynamic upper threshold = 2 * coverage > min. coverage threshold.
# '1' uses dynamic lower and upper threshold = mode +- 2 * std which is lower bound.
if sys.argv[1] == '-':
	depthFile = sys.stdin
	minThres = int(sys.argv[2])
	mode = int(sys.argv[3])
	entryName = sys.argv[4]
	outputDirectory = os.path.abspath(sys.argv[5])
	outputName = sys.argv[6]
else:
	depthFile = os.path.abspath(sys.argv[1])
	minThres = int(sys.argv[2])
	mode = int(sys.argv[3])
	entryName = sys.argv[4]
	outputDirectory = os.path.abspath(sys.argv[5])
	outputName = sys.argv[6]

# Get number of columns in depth file
nColumns = pd.read_table(depthFile, header=None, nrows=1).shape[1]

# Sets up to read depth file in chunks
depthFileIterator = pd.read_table(depthFile, header=None, iterator=True, chunksize=10000, usecols=range(2, nColumns))
depthDistribution = pd.DataFrame()

# Reads depth file in chunks turning multi-column observations into single-column, summing observations
for chunk in depthFileIterator:
	depthDistribution = pd.concat([depthDistribution, chunk.stack().value_counts()], axis=1, join='outer').sum(axis=1)

# Creates column of depth values, sets type of columns to int, names columns, sorts observations by count and resets index
depthDistribution = depthDistribution.reset_index().convert_dtypes()
depthDistribution.columns = ['Depth', 'Count']
depthDistribution = depthDistribution.sort_values(by='Count', ascending=False).reset_index(drop=True)

# Remove row of 0 depth observations
depthDistribution = depthDistribution[depthDistribution['Depth'] > 0]

# Calculate weighted unbiased summary statistics
depthDistributionWeightedStatsBessel = DescrStatsW(depthDistribution['Depth'], weights=depthDistribution['Count'], ddof=1)
mean = depthDistributionWeightedStatsBessel.mean
stdDev = depthDistributionWeightedStatsBessel.std

# Find highest observed coverage value
maxCov = depthDistribution.max(axis=0)['Depth']

if mode == 0:
	# Locate most observed coverage value fulfilling the criteria: 2 * coverage > min. coverage threshold.
	for row in depthDistribution.itertuples():
		if row.Depth * 2 > minThres:
			maxThres = round(row.Depth * 2)
			peakCoverage = row.Depth
			rank = row.Index
			break
if mode == 1:
	# Locate most observed coverage value fulfilling the criteria: coverage > min. coverage threshold.
	for row in depthDistribution.itertuples():
		if row.Depth > minThres:
			maxThres = round(row.Depth + 2 * stdDev)
			if row.Depth - 2 * stdDev > minThres:
				minThres = round(row.Depth - 2 * stdDev)
			peakCoverage = row.Depth
			rank = row.Index
			break

# Set upper threshold for depth plot
depthDistributionMaxThres = depthDistribution[depthDistribution['Depth'] <= maxThres]

# Create tsv file of summary values
with open(f'{outputDirectory}/{outputName}.tsv', 'w') as outfile:
	outfile.write(f'name\tmean\tstd\tpeak\tpeak_rank\tmin_coverage_threshold\tmax_coverage_threshold\tmax_coverage_value\n')
	outfile.write(f'{entryName}\t{mean}\t{stdDev}\t{peakCoverage}\t{rank}\t{minThres}\t{maxThres}\t{maxCov}\n')

# Create depth plot
plt.figure(figsize = (12, 6))
plt.bar(x=depthDistributionMaxThres['Depth'], height=depthDistributionMaxThres['Count'], color='r', width=1)
plt.ylabel("Count", fontweight = 'bold')
plt.xlabel("Coverage", fontweight = 'bold')
plt.xticks(np.arange(0, maxThres, 50.0), rotation='vertical')
plt.axvline(mean, color = 'black', linestyle = 'dashed', label = 'Mean')
plt.axvline(peakCoverage, color = 'blue', linestyle = 'dashed', label = 'Peak')
plt.axvline(minThres, color = 'lightgreen', linestyle = 'solid', label='Lower threshold')
plt.axvline(maxThres, color = 'green', linestyle = 'solid', label='Upper threshold')
plt.text(round(peakCoverage + stdDev), plt.ylim()[1]/2, f'Mean: {round(mean, ndigits=2)}x\nStd: {round(stdDev, ndigits=2)}\nPeak: {peakCoverage}x\nUpper threshold: {round(maxThres, ndigits=2)}x\nLower threshold: {round(minThres, ndigits=2)}x\nMax value: {maxCov}x', va = 'center', color = 'black', bbox=dict(facecolor='white'))
plt.legend()
plt.grid(True)
plt.title(f'{entryName} - Coverage distribution within upper threshold', fontweight = 'bold')
plt.tight_layout()
plt.savefig(f'{outputDirectory}/{outputName}.png')