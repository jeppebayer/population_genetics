#!/bin/env python3
import sys, os.path, matplotlib.pyplot as plt, pandas as pd, numpy as np

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
# Find highest observed coverage value
maxCov = depthDistribution.max(axis=0)['Depth']
if mode == 0:
	# Locate most observed coverage value fulfilling the criteria: 2 * coverage > min. coverage threshold.
	for row in depthDistribution.itertuples():
		if row.Depth * 2 > minThres:
			maxThres = row.Depth * 2
			coverageBin = row.Depth
			rank = row.Index
			break
if mode == 1:
	stdDev = depthDistribution['Depth'].std()
	# Locate most observed coverage value fulfilling the criteria: coverage > min. coverage threshold.
	for row in depthDistribution.itertuples():
		if row.Depth > minThres:
			maxThres = round(row.Depth + 2 * stdDev)
			if row.Depth - 2 * stdDev > minThres:
				minThres = round(row.Depth - 2 * stdDev)
			coverageBin = row.Depth
			rank = row.Index
			break

# Create dataframe with no values above the maximum threshold for more meaningful summary statistics
depthDistributionMaxThres = depthDistribution[depthDistribution['Depth'] <= maxThres]

# Calculate summary values
mean = depthDistributionMaxThres['Depth'].mean()
median = depthDistributionMaxThres['Depth'].median()
stdDev = depthDistributionMaxThres['Depth'].std()

# Create tsv file of summary values
with open(f'{plotDir}/{plotName}.tsv', 'w') as outfile:
	outfile.write(f'file\tmean\tmedian\tstd\tcoverage_bin\tbin_rank\tmin_coverage_threshold\tmax_coverage_threshold\tmax_coverage_value\n{plotName}\t{mean}\t{median}\t{stdDev}\t{coverageBin}\t{rank}\t{minThres}\t{maxThres}\t{maxCov}')

# Two subplots, one of the total coverage distribution and one within the maximum threshold
plt.figure(figsize = (12, 10))

tickspace = round(maxCov / 10, -2)
plt.subplot(2, 1, 1)
plt.bar(x=depthDistribution['Depth'], height=depthDistribution['Count'], color='r', width=1)
plt.ylabel("Count", fontweight = 'bold')
plt.xlabel("Coverage", fontweight = 'bold')
plt.xticks(np.arange(0, maxCov, tickspace), rotation='vertical')
plt.axvline(mean, color = 'black', linestyle = 'dashed', label = 'Mean')
plt.axvline(median, color = 'blue', linestyle = 'dashed', label = 'Median')
plt.axvline(minThres, color = 'green', linestyle = 'solid')
plt.text(minThres + 0.2, plt.ylim()[1]/2, f'Minimum threshold: {minThres}x', rotation = 90, va = 'center', color = 'green')
plt.axvline(maxThres, color = 'green', linestyle = 'solid')
plt.text(maxThres + 0.2, plt.ylim()[1]/2, f'Maximum threshold: {maxThres}x', rotation = 90, va = 'center', color = 'green')
plt.legend()
plt.grid(True)
plt.title('Coverage distribution', fontweight = 'bold')

plt.subplot(2, 1, 2)
plt.bar(x=depthDistributionMaxThres['Depth'], height=depthDistributionMaxThres['Count'], color='r', width=1)
plt.ylabel("Count", fontweight = 'bold')
plt.xlabel("Coverage", fontweight = 'bold')
plt.xticks(np.arange(0, maxThres, 100.0), rotation='vertical')
plt.axvline(mean, color = 'black', linestyle = 'dashed', label = 'Mean')
plt.axvline(median, color = 'blue', linestyle = 'dashed', label = 'Median')
plt.axvline(minThres, color = 'green', linestyle = 'solid')
plt.text(minThres + 0.2, plt.ylim()[1]/2, f'Minimum threshold: {minThres}x', rotation = 90, va = 'center', color = 'green')
plt.axvline(maxThres, color = 'green', linestyle = 'solid')
plt.text(maxThres + 0.2, plt.ylim()[1]/2, f'Maximum threshold: {maxThres}x', rotation = 90, va = 'center', color = 'green')
plt.legend()
plt.grid(True)
plt.title('Coverage distribution within maximum threshold', fontweight = 'bold')

plt.tight_layout()
plt.savefig(f'{plotDir}/{plotName}.png')