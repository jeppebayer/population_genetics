#!/bin/env python3
import sys, os.path, matplotlib.pyplot as plt, pandas as pd, numpy as np
from statsmodels.stats.weightstats import DescrStatsW

# Takes in positional arguments
if sys.argv[1] == '-':
	depthFile = sys.stdin
	entryName = sys.argv[2]
	outputName = sys.argv[3]
else:
	depthFile = os.path.abspath(sys.argv[1])
	entryName = sys.argv[2]
	outputName = sys.argv[3]

# Sets up to read depthFile in chunks
depthFile = pd.read_table(depthFile, header=None, names=('Chromosome', 'Position', 'Depth'))

for index, chromosome in enumerate(pd.unique(depthFile['Chromosome']), 1):
	if index == 1:
		lastPos = depthFile[depthFile['Chromosome'] == chromosome].max(axis=0)['Position']
		continue
	depthFile.loc[depthFile['Chromosome'] == chromosome, 'Position'] = depthFile.loc[depthFile['Chromosome'] == chromosome, 'Position'] + lastPos
	lastPos = depthFile[depthFile['Chromosome'] == chromosome].max(axis=0)['Position']

# Create depth distribution plot with fraction indicators
plt.figure(figsize = (12, 6))
for index, chromosome in enumerate(pd.unique(depthFile['Chromosome'])):
	if index % 2 == 0:
		color = 'lightblue'
	else:
		color = 'steelblue'
	plt.plot(depthFile.loc[depthFile['Chromosome'] == chromosome, 'Position'], depthFile.loc[depthFile['Chromosome'] == chromosome, 'Depth'], label = chromosome, color = color)
plt.ylabel("Coverage", fontweight = 'bold')
plt.xlabel("Position", fontweight = 'bold')
plt.legend(loc = 'upper right')
plt.grid(True)
plt.title(f'{entryName} - Coverage along reference genome', fontweight = 'bold')
plt.tight_layout()
plt.savefig(f'{outputName}.png')