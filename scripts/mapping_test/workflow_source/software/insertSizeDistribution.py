#!/bin/env python3
import sys, os.path, matplotlib.pyplot as plt, pandas as pd, numpy as np
from statsmodels.stats.weightstats import DescrStatsW

# Takes in positional arguments
if sys.argv[1] == '-':
	insertFile = sys.stdin
	entryName = sys.argv[2]
	outputName = sys.argv[3]
else:
	insertFile = os.path.abspath(sys.argv[1])
	entryName = sys.argv[2]
	outputName = sys.argv[3]

insertDistribution = pd.read_table(insertFile, header=None, usecols=[1]).value_counts()

# # Sets up to read depthFile in chunks
# insertFileIterator = pd.read_table(insertFile, header=None, iterator=True, chunksize=10000, usecols=[1])
# insertDistribution = pd.DataFrame()

# # Reads depthFile in chunks summing observations
# for chunk in insertFileIterator:
# 	insertDistribution = pd.concat([insertDistribution, chunk.value_counts()], axis=1, join='outer').sum(axis=1)

print(insertDistribution)

# Creates column of depth values, sets type of columns to int, names columns, sorts observations by count and resets index
insertDistribution = insertDistribution.reset_index().convert_dtypes()
print(insertDistribution)
insertDistribution.columns = ['Size', 'Count']
print(insertDistribution)
insertDistribution = insertDistribution.sort_values(by='Count', ascending=False).reset_index(drop=True)
print(insertDistribution)

# # Removes row of 0 depth observations
# insertDistribution = insertDistribution[insertDistribution['Size'] > 0]

# Creates weighted unbiased summary statistics 
insertDistributionWeightedStatsBessel = DescrStatsW(insertDistribution['Size'], weights=insertDistribution['Count'], ddof=1)
mean = insertDistributionWeightedStatsBessel.mean
stdDev = insertDistributionWeightedStatsBessel.std

print(mean)
print(stdDev)

# Find highest observed coverage value
maxSize = insertDistribution.max(axis=0)['Size']

# Calculate total amount of observations
totalCount = insertDistribution.sum(axis=0)['Count']

# Sort observations by depth values
insertDistribution.sort_values(by='Size', ascending=False, inplace=True)

# Create depth distribution plot with fraction indicators
plt.figure(figsize = (12, 6))
plt.bar(x=insertDistribution['Size'], height=insertDistribution['Count'], color='r', width=1)
plt.ylabel("Count", fontweight = 'bold')
plt.xlabel("Size", fontweight = 'bold')
plt.xticks(np.arange(0, maxSize, 10.0), rotation='vertical')
plt.axvline(mean, color = 'black', linestyle = 'solid', label = 'Mean')
plt.text(round(mean + stdDev), plt.ylim()[1]/2, f'Mean: {round(mean, ndigits=2)}x\nStd: {round(stdDev, ndigits=2)}\nMax value: {maxSize}x', va = 'center', color = 'black', bbox=dict(facecolor='white'))
plt.legend()
plt.grid(True)
plt.title(f'{entryName} - Insert size distribution', fontweight = 'bold')
plt.tight_layout()
plt.savefig(f'{outputName}.png')