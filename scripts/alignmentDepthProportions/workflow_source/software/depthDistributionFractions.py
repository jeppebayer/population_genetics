#!/bin/env python3
import sys, os.path, matplotlib.pyplot as plt, pandas as pd, numpy as np
from statsmodels.stats.weightstats import DescrStatsW

colorList = ['darkgreen', 'green', 'lightgreen', 'paleturquoise', 'turquoise', 'darkturquoise', 'deepskyblue', 'steelblue', 'dodgerblue', 'royalblue', 'blue', 'mediumpurple', 'darkorchid', 'indigo']

if sys.argv[1] == '-':
	depthFile = sys.stdin
	entryName = sys.argv[2]
	lowerThres = int(sys.argv[3])
	outputName = sys.argv[4]
else:
	depthFile = os.path.abspath(sys.argv[1])
	entryName = sys.argv[2]
	lowerThres = int(sys.argv[3])
	outputName = sys.argv[4]

depthFileIterator = pd.read_table(depthFile, header=None, iterator=True, chunksize=10000, usecols=[2])
depthDistribution = pd.DataFrame()

for chunk in depthFileIterator:
	depthDistribution = pd.concat([depthDistribution, chunk.value_counts()], axis=1, join='outer').sum(axis=1)

depthDistribution = depthDistribution.reset_index().convert_dtypes()
depthDistribution.columns = ['Depth', 'Count']
depthDistribution = depthDistribution.sort_values(by='Count', ascending=False).reset_index(drop=True)

depthDistribution = depthDistribution[depthDistribution['Depth'] > 0]

depthDistributionWeightedStatsBessel = DescrStatsW(depthDistribution['Depth'], weights=depthDistribution['Count'], ddof=1)
mean = depthDistributionWeightedStatsBessel.mean
stdDev = depthDistributionWeightedStatsBessel.std

maxCov = depthDistribution.max(axis=0)['Depth']

totalCount = depthDistribution.sum(axis=0)['Count']

for row in depthDistribution.itertuples():
	if row.Depth > lowerThres:
		peakCoverage = row.Depth
		break

depthDistribution.sort_values(by='Depth', ascending=False, inplace=True)

fractions = {0.98: None, 0.95: None, 0.90: None,0.85: None, 0.80: None,  0.75: None, 0.70: None, 0.65: None, 0.60: None, 0.55: None, 0.50: None, 0.45: None, 0.40: None}
sumCount = 0
for fraction in fractions:
	for row in depthDistribution.itertuples():
		sumCount += row.Count
		if sumCount >= totalCount * fraction:
			fractions[fraction] = row.Depth
			sumCount = 0
			break

upperThres = peakCoverage + 2 * stdDev
depthDistributionUpperThres = depthDistribution[depthDistribution['Depth'] <= upperThres]

resultString = f'{entryName}\t{mean}\t{stdDev}\t{maxCov}\t{peakCoverage}\t{lowerThres}\t' + '\t'.join([str(fractions[i]) for i in fractions]) + '\n'

with open(f'{outputName}.tsv', 'w') as outfile:
	outfile.write('name\tmean\tstd\tmax_value\tpeak\tlower_threshold\t98%\t95%\t90%\t85%\t80%\t75%\t70%\t65%\t60%\t55%\t50%\t45%\t40%\n')
	outfile.write(resultString)

plt.figure(figsize = (12, 6))

plt.bar(x=depthDistributionUpperThres['Depth'], height=depthDistributionUpperThres['Count'], color='r', width=1)
plt.ylabel("Count", fontweight = 'bold')
plt.xlabel("Coverage", fontweight = 'bold')
plt.xticks(np.arange(0, upperThres, 100.0), rotation='vertical')
plt.axvline(mean, color = 'black', linestyle = 'solid', label = 'Mean')
plt.axvline(peakCoverage, color = 'gold', linestyle = 'solid', label = 'Peak')
for i, j in enumerate(fractions):
	if fractions[j] > upperThres:
		continue
	plt.axvline(fractions[j], color=colorList[i], linestyle='dashed', label=f'{f'{j}0' if len(str(j)) < 4 else j}%')
plt.text(round(peakCoverage + stdDev), plt.ylim()[1]/2, f'Mean: {round(mean, ndigits=2)}x\nStd: {round(stdDev, ndigits=2)}\nPeak: {peakCoverage}x\nMax value: {maxCov}x', va = 'center', color = 'black', bbox=dict(facecolor='white'))
plt.legend()
plt.grid(True)
plt.title(f'{entryName} - Coverage distribution within upper threshold', fontweight = 'bold')

plt.tight_layout()
plt.savefig(f'{outputName}.png')