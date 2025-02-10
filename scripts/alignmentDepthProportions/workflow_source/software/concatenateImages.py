#!/bin/env python3
from PIL import Image
import sys, os

outputDirectory = os.path.abspath(sys.argv[1])
outputName = sys.argv[2]
nCols = int(sys.argv[3])
images = sys.argv[4:]

imageList = [Image.open(image) for image in images]
nImages = len(imageList)
nRows = 1 + ((nImages - 1) // nCols)
baseWidth, baseHeight = imageList[0].size
imageListResize = [image.resize((baseWidth, baseHeight)) if image.width != baseWidth or image.height != baseHeight else image for image in imageList]

width = nCols * baseWidth
height = nRows * baseHeight

combinedImage = Image.new('RGBA', (width, height), 'white')

widthPos = 0
heightPos = 0
startIndex = 0
endIndex = nCols
for row in range(nRows):
	for col in range(startIndex, endIndex):
		if col >= nImages:
			break
		combinedImage.paste(imageListResize[col], (widthPos, heightPos))
		widthPos += baseWidth
	heightPos += baseHeight
	widthPos = 0
	startIndex = endIndex
	endIndex += nCols

combinedImage.save(f'{outputDirectory}/{outputName}', 'PNG')