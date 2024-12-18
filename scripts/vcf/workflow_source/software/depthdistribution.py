#!/bin/env python3
import sys, os.path, matplotlib.pyplot as plt, pandas as pd, numpy as np

# First positional argument is used as a static lower coverage threshold
minthres = int(sys.argv[1])
# Second positional argument indicates which of two threshold modes to use.
# '0' uses static lower threshold and dynamic upper threshold = 2 * coverage > min. coverage threshold.
# '1' uses dynamic lower and upper threshold = mode +- 2 * std which is lower bound.
mode = int(sys.argv[2])

fullpath = os.path.abspath(sys.argv[3])
plotdir = os.path.abspath(sys.argv[4])
plotname = os.path.basename(fullpath)
depthfile = pd.read_table(fullpath, header=None)

# Transform samtools multi-column depth file into single-column dataframe of depth values.
ncolumns = depthfile.shape[1]
depthfilereformat = pd.DataFrame()
for i in range(2, ncolumns):
	depthfilereformat = depthfilereformat._append(depthfile[[i]].rename(columns={i:0}), ignore_index=True)

# Remove all depth observations with value 0
depthfilereformat = depthfilereformat[depthfilereformat[0] > 0]
# Find highest observed coverage value
maxcov = np.max(depthfilereformat[0])

# Create coverage distribution.
bins = np.bincount(depthfilereformat[0])
# Explicitly attach coverage value to number of observation.
namedbins = [[i, bins[i]] for i in range(len(bins))]
# Sort in decending number of observations.
namedbins.sort(key = lambda count: count[1], reverse = True)

if mode == 0:
	# Locate most observed coverage value fulfilling the criteria: 2 * coverage > min. coverage threshold.
	for i, j in enumerate(namedbins):
		if int(j[0]) * 2 > minthres:
			maxthres = int(j[0]) * 2
			coveragebin = j[0]
			rank = i + 1
			break

if mode == 1:
	stddev = np.std(depthfilereformat[0])
	# Locate most observed coverage value fulfilling the criteria: coverage > min. coverage threshold.
	for i, j in enumerate(namedbins):
		if int(j[0]) > minthres:
			maxthres = round(int(j[0]) + 2 * stddev)
			if int(j[0]) - 2 * stddev > minthres:
				minthres = round(int(j[0]) - 2 * stddev)
			coveragebin = j[0]
			rank = i + 1
			break

# Create dataframe with no values above the maximum threshold for more meaningful summary statistics
depthfilereformatmaxthres = depthfilereformat[depthfilereformat[0] <= maxthres]
# Calculate summary values
mean = np.mean(depthfilereformatmaxthres[0])
median = np.median(depthfilereformatmaxthres[0])
stddev = np.std(depthfilereformatmaxthres[0])

# Create tsv file of summary values
with open(f'{plotdir}/{plotname}.tsv', 'w') as outfile:
	outfile.write(f'file\tmean\tmedian\tstd\tcoverage_bin\tbin_rank\tmin_coverage_threshold\tmax_coverage_threshold\tmax_coverage_value\n{plotname}\t{mean}\t{median}\t{stddev}\t{coveragebin}\t{rank}\t{minthres}\t{maxthres}\t{maxcov}')

# Two subplots, one of the total coverage distribution and one within the maximum threshold
plt.figure(figsize = (12, 10))

tickspace = round(maxcov / 10, -2)
plt.subplot(2, 1, 1)
plt.hist(depthfilereformat[0], bins = maxcov, color = 'r')
plt.ylabel("Count", fontweight = 'bold')
plt.xlabel("Coverage", fontweight = 'bold')
plt.xticks(np.arange(0, maxcov, tickspace), rotation='vertical')
plt.axvline(mean, color = 'black', linestyle = 'dashed', label = 'Mean')
plt.axvline(median, color = 'blue', linestyle = 'dashed', label = 'Median')
plt.axvline(minthres, color = 'green', linestyle = 'solid')
plt.text(minthres + 0.2, plt.ylim()[1]/2, f'Minimum threshold: {minthres}x', rotation = 90, va = 'center', color = 'green')
plt.axvline(maxthres, color = 'green', linestyle = 'solid')
plt.text(maxthres + 0.2, plt.ylim()[1]/2, f'Maximum threshold: {maxthres}x', rotation = 90, va = 'center', color = 'green')
plt.legend()
plt.grid(True)
plt.title('Coverage distribution', fontweight = 'bold')

plt.subplot(2, 1, 2)
plt.hist(depthfilereformatmaxthres[0], bins = maxthres, color = 'r')
plt.ylabel("Count", fontweight = 'bold')
plt.xlabel("Coverage", fontweight = 'bold')
plt.xticks(np.arange(0, maxthres, 100.0), rotation='vertical')
plt.axvline(mean, color = 'black', linestyle = 'dashed', label = 'Mean')
plt.axvline(median, color = 'blue', linestyle = 'dashed', label = 'Median')
plt.axvline(minthres, color = 'green', linestyle = 'solid')
plt.text(minthres + 0.2, plt.ylim()[1]/2, f'Minimum threshold: {minthres}x', rotation = 90, va = 'center', color = 'green')
plt.axvline(maxthres, color = 'green', linestyle = 'solid')
plt.text(maxthres + 0.2, plt.ylim()[1]/2, f'Maximum threshold: {maxthres}x', rotation = 90, va = 'center', color = 'green')
plt.legend()
plt.grid(True)
plt.title('Coverage distribution within maximum threshold', fontweight = 'bold')

plt.tight_layout()
plt.savefig(f'{plotdir}/{plotname}.png')

# A single plot within the maximum threshold
# plt.figure(figsize = (15, 8.5))
# plt.hist(depthfilereformat[0], bins = maxthres, color = 'r')
# plt.ylabel("Count")
# plt.xlabel("Coverage")
# plt.xticks(np.arange(0, maxthres, 100.0), rotation='vertical')
# plt.axvline(mean, color = 'black', linestyle = 'dashed', label = 'Mean')
# plt.axvline(median, color = 'blue', linestyle = 'dashed', label = 'Median')
# plt.axvline(minthres, color = 'green', linestyle = 'solid')
# plt.text(minthres + 0.1, plt.ylim()[1]/2, f'Minimum threshold: {minthres}x', rotation = 90, va = 'center', color = 'green')
# plt.axvline(maxthres, color = 'green', linestyle = 'solid')
# plt.text(maxthres + 0.1, plt.ylim()[1]/2, f'Maximum threshold: {maxthres}x', rotation = 90, va = 'center', color = 'green')
# plt.legend()
# plt.grid(True)
# plt.title('Coverage distribution within maximum threshold')

# plt.savefig(f'{plotdir}/{plotname}.png')