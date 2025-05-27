#!/bin/env python3
import sys, os, pandas as pd, matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

usage = f'{os.path.basename(__file__)} <sfsFile> <outputDirectory>'

def sfsPlot():
	sfsFile = pd.read_table(os.path.abspath(sys.argv[1]), sep='\t', header=None, names=['Sample', 'Bin', 'Count'])
	outputDirectory = sys.argv[2]

	with PdfPages(f'{outputDirectory}/{os.path.basename(sys.argv[1])}.pdf') as pdf:
		for sample in sfsFile['Sample'].unique():
			subsetSfsFile = sfsFile[(sfsFile['Sample'] == sample) & (sfsFile['Bin'] > 0)]
			plt.figure(figsize = (6, 6))
			plt.bar(x=subsetSfsFile['Bin'], height=subsetSfsFile['Count'], color='k', width=1)
			plt.ylabel('Count', fontweight='bold')
			plt.xlabel('Bin', fontweight='bold')
			plt.xticks(range(0, len(subsetSfsFile) + 1, 5))
			plt.title(f'SFS {sample}')
			plt.tight_layout()
			pdf.savefig()
			plt.close()

if __name__ == "__main__":
	if len(sys.argv) != 3:
		print(usage)
		sys.exit(1)
	sfsPlot()