#!/bin/env python3
import sys, os.path, matplotlib.pyplot as plt, pandas as pd, numpy as np
from statsmodels.stats.weightstats import DescrStatsW

# First positional argument is used as a static lower coverage threshold
minThres = int(sys.argv[1])
# Second positional argument indicates which of two threshold modes to use.
# '0' uses static lower threshold and dynamic upper threshold = 2 * coverage > min. coverage threshold.
# '1' uses dynamic lower and upper threshold = mode +- 2 * std which is lower bound.
mode = int(sys.argv[2])

fullPath = os.path.abspath(sys.argv[3])
plotDir = os.path.abspath(sys.argv[4])
plotName = os.path.basename(fullPath)

# Transform samtools multi-column depth file into single-column dataframe of depth values.
nColumns = pd.read_table(fullPath, header=None, nrows=1).shape[1]
depthFileIterator = pd.read_table(fullPath, header=None, iterator=True, chunksize=10000, usecols=range(2, nColumns))
depthDistribution = pd.DataFrame()

for chunk in depthFileIterator:
	depthDistribution = pd.concat([depthDistribution, chunk.stack().value_counts()], axis=1, join='outer').sum(axis=1)

# Create observation column from index
depthDistribution = depthDistribution.reset_index().convert_dtypes()
depthDistribution.columns = ['Depth', 'Count']
depthDistribution = depthDistribution.sort_values(by='Count', ascending=False).reset_index(drop=True)

# Remove depth 0 observations
depthDistribution = depthDistribution[depthDistribution['Depth'] > 0]

# Calculate summary statistics
depthDistributionweightedStatsBessel = DescrStatsW(depthDistribution['Depth'], weights=depthDistribution['Count'], ddof=1)
mean = depthDistributionweightedStatsBessel.mean
stdDev = depthDistributionweightedStatsBessel.std

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

# Create dataframe with no values above the maximum threshold for more meaningful summary statistics
depthDistributionMaxThres = depthDistribution[depthDistribution['Depth'] <= maxThres]

# Create tsv file of summary values
with open(f'{plotDir}/{plotName}.tsv', 'w') as outfile:
	outfile.write(f'file\tmean\tstd\tpeak\tpeak_rank\tmin_coverage_threshold\tmax_coverage_threshold\tmax_coverage_value\n{plotName}\t{mean}\t{stdDev}\t{peakCoverage}\t{rank}\t{minThres}\t{maxThres}\t{maxCov}')

# Create plot
plt.figure(figsize = (12, 6))

plt.bar(x=depthDistributionMaxThres['Depth'], height=depthDistributionMaxThres['Count'], color='r', width=1)
plt.ylabel("Count", fontweight = 'bold')
plt.xlabel("Coverage", fontweight = 'bold')
plt.xticks(np.arange(0, maxThres, 100.0), rotation='vertical')
plt.axvline(mean, color = 'black', linestyle = 'dashed', label = 'Mean')
plt.axvline(peakCoverage, color = 'blue', linestyle = 'dashed', label = 'Peak')
plt.axvline(minThres, color = 'lightgreen', linestyle = 'solid', label='Lower threshold')
plt.axvline(maxThres, color = 'green', linestyle = 'solid', label='Upper threshold')
plt.text(round(peakCoverage + stdDev), plt.ylim()[1]/2, f'Mean: {round(mean, ndigits=2)}x\nStd: {round(stdDev, ndigits=2)}\nPeak: {peakCoverage}x\nUpper threshold: {maxThres}x\nLower threshold: {minThres}x\nMax value: {maxCov}x', va = 'center', color = 'black', bbox=dict(facecolor='white'))
plt.legend()
plt.grid(True)
plt.title('Coverage distribution within maximum threshold', fontweight = 'bold')

plt.tight_layout()
plt.savefig(f'{plotDir}/{plotName}.png')