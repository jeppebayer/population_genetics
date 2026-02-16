#!/bin/env python3
from PIL import Image
import sys, os

# Get positional arguments
outputDirectory = os.path.abspath(sys.argv[1])
outputName = sys.argv[2]
nCols = int(sys.argv[3])
images = sys.argv[4:]

# Store images in list
imageList = [Image.open(image) for image in images]
# Get number of images
nImages = len(imageList)
# Determine number of rows
nRows = 1 + ((nImages - 1) // nCols)
# Set a base width and height for each images based on the first image in the list
baseWidth, baseHeight = imageList[0].size
# Resize any image that differs from the base width and height
imageListResize = [image.resize((baseWidth, baseHeight)) if image.width != baseWidth or image.height != baseHeight else image for image in imageList]

# Set width and height of concatenated image
if nImages < nCols:
	width = nImages * baseWidth
else:
	width = nCols * baseWidth
height = nRows * baseHeight

# Create empty image
combinedImage = Image.new('RGBA', (width, height), 'white')

# Loops through images concatenating them in columns and rows
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

# Saves new combined image
combinedImage.save(f'{outputDirectory}/{outputName}', 'PNG')