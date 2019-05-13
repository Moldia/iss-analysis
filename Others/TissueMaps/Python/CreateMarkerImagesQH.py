# File: CreateMarkerImage.py
# Generate full size images with marker symbols
# Petter Ranefall 2014

# Import neccessary modules:
import os
import shutil
from PIL import Image
import subprocess
import numpy as np
from ..VipsSystemWrapper import VipsSystemWrapper as vsw
from ..ImageReader import ImageReader as imageReader

class CreateMarkerImages(object):

    def __init__(self):
        pass

    def insertAllMarkersToNewImages(self, baseImageFileName, posFileName, tagListFileName, outputBaseFileName, iconFolder, createTiles = False, tileSize = 1024, includeNNNN = False, nnnnSymbol = 'redStar', qualityT = 0.5, markerScale = 2):
        success = True
        try:
            statsFp = open(outputBaseFileName + '_stats.csv', 'w')
            baseImageFileName = baseImageFileName.replace('"', '')
            baseImageFileName = baseImageFileName.replace('\\', '/')
            
            posFileName = posFileName.replace('"', '')
            posFileName = posFileName.replace('\\', '/')
            
            tagListFileName = tagListFileName.replace('"', '')
            tagListFileName = tagListFileName.replace('\\', '/')
            
            outputBaseFileName = outputBaseFileName.replace('"', '')
            outputBaseFileName = outputBaseFileName.replace('\\', '/')

            baseFolder = os.path.split(outputBaseFileName)[0]
            if not baseFolder.endswith('/'):
                baseFolder = baseFolder + '/'
            
            iconFolder = iconFolder.replace('\\', '/')
            iconFolder = iconFolder.replace('"', '')
            if not iconFolder.endswith('/'):
                iconFolder = iconFolder + '/'
            
            #read the posFile
            positions = np.genfromtxt(posFileName, delimiter=',', dtype=None)
            lenPositions = len(positions[:,0])
            positions = positions[1:lenPositions,:] # Remove header
            allTags = positions[:,1]
            usedIndices = [index for index, value in enumerate(allTags) if value != 'NNNN']
            usedPositions = positions[usedIndices,:]
            usedLetters = usedPositions[:,0]
            
            uniqueLetters = np.unique(usedLetters)
            #read the taglist file
            taglist = np.genfromtxt(tagListFileName, delimiter=',', dtype=None)
            usedTagIndices = [index for index, value in enumerate(taglist) if value[0] in uniqueLetters]
            usedTags = taglist[usedTagIndices]
            if createTiles:
                vswInstance = vsw.VipsSystemWrapper()

            imageReaderInstance = imageReader.ImageReader()
            imageSize = imageReaderInstance.getImageSize(baseImageFileName)
            iconScaleFactor = markerScale
            for tag in usedTags:
                symbol = tag[2]
                symbolFile = iconFolder + symbol + '.png'
                iconImage = Image.open(symbolFile)
                iconImageSize = iconImage.size
                halfIconImageSize = (iconImageSize[0] / 2, iconImageSize[1] / 2)
                iconImageSizeScaled = (int(iconImageSize[0] * iconScaleFactor), int(iconImageSize[1] * iconScaleFactor))
                iconImageScaled = iconImage.resize(iconImageSizeScaled)
                halfIconImageSizeScaled = (iconImageSizeScaled[0] / 2, iconImageSizeScaled[1] / 2)
                name = tag[1]
                currentPositionIndices = [index for index, value in enumerate(usedLetters) if value == tag[0]]
                currentPositions = usedPositions[currentPositionIndices,2:6]
                x = map(int, currentPositions[:,0])
                y = map(int, currentPositions[:,1])
                q = map(float, currentPositions[:,3])
                lenxy = len(x)
                newImage = Image.new("RGB", imageSize)
                statsFp.write(symbol + "," + name + "," + str(lenxy) + "\n")
                for i in range(0, lenxy):
                    if (q[i] >= qualityT):
                        pos = (x[i] - halfIconImageSizeScaled[0], y[i] - halfIconImageSizeScaled[1])
                        newImage.paste(iconImageScaled, pos)
                    
                #create mask and add alpha channel
                lut = [255 if v > 0 else 0 for v in range(256)]
                mask = newImage.convert("L") 
                mask = mask.point(lut) 
                newImage.putalpha(mask)
                workingFile = outputBaseFileName + '_' + symbol + '_' + name + '.png'
                newImage.save(workingFile)
                if createTiles:
                    workingFolder = outputBaseFileName + '_' + symbol + '_' + name
                    if imageSize[0] > tileSize or imageSize[1] > tileSize:
                        vswInstance.createTiles(inputFile = workingFile, outputFolder = workingFolder, tileSize = tileSize, overlap = 0, suffix = '.png')
                    else:
                        outputFolderName = outputBaseFileName + '_' + symbol + '_' + name + '_files'
                        if os.path.isdir(outputFolderName):
                            shutil.rmtree(outputFolderName)
                        os.mkdir(outputFolderName)
                        os.mkdir(outputFolderName + '/0')
                        newImage.save(outputFolderName + '/0/0_0.png')
                        
                    
            if includeNNNN:
                nnnnIndices = [index for index, value in enumerate(allTags) if value == 'NNNN']
                nnnnPositions = positions[nnnnIndices,:]
                nnnnTagIndices = [index for index, value in enumerate(taglist) if value[0] == 'NNNN']
                nnnnTags = taglist[nnnnTagIndices]
                if (len(nnnnTags) > 0):
                    symbol = nnnnTags[0][2]
                else:
                    symbol = nnnnSymbol
                symbolFile = iconFolder + symbol + '.png'
                iconImage = Image.open(symbolFile)
                iconImageSize = iconImage.size
                halfIconImageSize = (iconImageSize[0] / 2, iconImageSize[1] / 2)
                iconImageSizeScaled = (int(iconImageSize[0] * iconScaleFactor), int(iconImageSize[1] * iconScaleFactor))
                iconImageScaled = iconImage.resize(iconImageSizeScaled)
                halfIconImageSizeScaled = (iconImageSizeScaled[0] / 2, iconImageSizeScaled[1] / 2)
                name = 'NNNN'
                currentPositions = nnnnPositions[:,2:4]
                x = map(int, currentPositions[:,0])
                y = map(int, currentPositions[:,1])
                lenxy = len(x)
                newImage = Image.new("RGB", imageSize)
                statsFp.write(symbol + "," + name + "," + str(lenxy) + "\n")
                for i in range(0, lenxy):
                    pos = (x[i] - halfIconImageSizeScaled[0], y[i] - halfIconImageSizeScaled[1])
                    newImage.paste(iconImageScaled, pos)
                    
                #create mask and add alpha channel
                lut = [255 if v > 0 else 0 for v in range(256)]
                mask = newImage.convert("L") 
                mask = mask.point(lut) 
                newImage.putalpha(mask)
                workingFile = outputBaseFileName + '_' + symbol + '_' + name + '.png'
                newImage.save(workingFile)
                if createTiles:
                    workingFolder = outputBaseFileName + '_' + symbol + '_' + name
                    if imageSize[0] > tileSize or imageSize[1] > tileSize:
                        vswInstance.createTiles(inputFile = workingFile, outputFolder = workingFolder, tileSize = tileSize, overlap = 0, suffix = '.png')
                    else:
                        outputFolderName = outputBaseFileName + '_' + symbol + '_' + name + '_files'
                        if os.path.isdir(outputFolderName):
                            shutil.rmtree(outputFolderName)
                        os.mkdir(outputFolderName)
                        os.mkdir(outputFolderName + '/0')
                        newImage.save(outputFolderName + '/0/0_0.png')
                
            statsFp.close()
              
        except Exception as e:
            print("Exception: %s" % e)
            success = False
            pass
        
        return success
    

    def qualityHeatMap(self, baseImageFileName, posFileName, tagListFileName, outputFileName, createTiles = False, tileSize = 1024, includeNNNN = False):
        success = True
        try:
            baseImageFileName = baseImageFileName.replace('"', '')
            baseImageFileName = baseImageFileName.replace('\\', '/')
            
            posFileName = posFileName.replace('"', '')
            posFileName = posFileName.replace('\\', '/')
            
            tagListFileName = tagListFileName.replace('"', '')
            tagListFileName = tagListFileName.replace('\\', '/')
            
            outputFileName = outputBaseFileName.replace('"', '')
            outputFileName = outputBaseFileName.replace('\\', '/')

            baseFolder = os.path.split(outputFileName)[0]
            if not baseFolder.endswith('/'):
                baseFolder = baseFolder + '/'
            
            #read the posFile
            positions = np.genfromtxt(posFileName, delimiter=',', dtype=None)
            lenPositions = len(positions[:,0])
            positions = positions[1:lenPositions,:] # Remove header
            allTags = positions[:,1]
            usedIndices = [index for index, value in enumerate(allTags) if value != 'NNNN']
            usedPositions = positions[usedIndices,:]
            usedLetters = usedPositions[:,0]
            
            uniqueLetters = np.unique(usedLetters)
            #read the taglist file
            taglist = np.genfromtxt(tagListFileName, delimiter=',', dtype=None)
            usedTagIndices = [index for index, value in enumerate(taglist) if value[0] in uniqueLetters]
            usedTags = taglist[usedTagIndices]
            if createTiles:
                vswInstance = vsw.VipsSystemWrapper()

            imageReaderInstance = imageReader.ImageReader()
            imageSize = imageReaderInstance.getImageSize(baseImageFileName)
            newImage = Image.new("B", imageSize)
            nTiles = (int(imageSize.x/tileSize), int(imageSize.y/tileSize))
            tileMeans = np.zeros(nTiles)
            tileNs = np.zeros(nTiles)
            #TODO: Create 2D arrays for each tile: sum, count => compute mean
            for tag in usedTags:
                name = tag[1]
                currentPositionIndices = [index for index, value in enumerate(usedLetters) if value == tag[0]]
                currentPositions = usedPositions[currentPositionIndices,2:6]
                x = map(int, currentPositions[:,0])
                y = map(int, currentPositions[:,1])
                q = map(float, currentPositions[:,3])
                lenxy = len(x)
                for i in range(0, lenxy):
                    tileX = int(x[i]/tileSize)
                    tileY = int(y[i]/tileSize)
                    tileMeans[tileX, tileY] += q[i]
                    tileNs[tileX, tileY] = tileNs[tileX, tileY] + 1
                        
                    
##            if includeNNNN:
##                nnnnIndices = [index for index, value in enumerate(allTags) if value == 'NNNN']
##                nnnnPositions = positions[nnnnIndices,:]
##                nnnnTagIndices = [index for index, value in enumerate(taglist) if value[0] == 'NNNN']
##                nnnnTags = taglist[nnnnTagIndices]
##                if (len(nnnnTags) > 0):
##                    symbol = nnnnTags[0][2]
##                else:
##                    symbol = nnnnSymbol
##                symbolFile = iconFolder + symbol + '.png'
##                iconImage = Image.open(symbolFile)
##                iconImageSize = iconImage.size
##                halfIconImageSize = (iconImageSize[0] / 2, iconImageSize[1] / 2)
##                iconImageSizeScaled = (int(iconImageSize[0] * iconScaleFactor), int(iconImageSize[1] * iconScaleFactor))
##                iconImageScaled = iconImage.resize(iconImageSizeScaled)
##                halfIconImageSizeScaled = (iconImageSizeScaled[0] / 2, iconImageSizeScaled[1] / 2)
##                name = 'NNNN'
##                currentPositions = nnnnPositions[:,2:4]
##                x = map(int, currentPositions[:,0])
##                y = map(int, currentPositions[:,1])
##                lenxy = len(x)
##                newImage = Image.new("RGB", imageSize)
##                statsFp.write(symbol + "," + name + "," + str(lenxy) + "\n")
##                for i in range(0, lenxy):
##                    pos = (x[i] - halfIconImageSizeScaled[0], y[i] - halfIconImageSizeScaled[1])
##                    newImage.paste(iconImageScaled, pos)
##                    
##                #create mask and add alpha channel
##                lut = [255 if v > 0 else 0 for v in range(256)]
##                mask = newImage.convert("L") 
##                mask = mask.point(lut) 
##                newImage.putalpha(mask)
##                workingFile = outputBaseFileName + '_' + symbol + '_' + name + '.png'
##                newImage.save(workingFile)
##                if createTiles:
##                    workingFolder = outputBaseFileName + '_' + symbol + '_' + name
##                    if imageSize[0] > tileSize or imageSize[1] > tileSize:
##                        vswInstance.createTiles(inputFile = workingFile, outputFolder = workingFolder, tileSize = tileSize, overlap = 0, suffix = '.png')
##                    else:
##                        outputFolderName = outputBaseFileName + '_' + symbol + '_' + name + '_files'
##                        if os.path.isdir(outputFolderName):
##                            shutil.rmtree(outputFolderName)
##                        os.mkdir(outputFolderName)
##                        os.mkdir(outputFolderName + '/0')
##                        newImage.save(outputFolderName + '/0/0_0.png')



            tileNs = max(tileNs, 1)
            tileMeans = map(int, (255.0 * tileMeans)/tileNs + 0.5)

            for ix in range(0, nTiles[0]):
                for iy in range(0, nTiles[1]):
                    if tileMeans[ix, iy] > 0:
                        im = Image.new("B", (tileSize, tileSize), tileMeans[ix, iy]);
                        newImage.paste(im, (ix * tileSize, iy * tileSize))
            

            #TODO: Create filled tiles and insert into image
            
            newImage.save(outputFileName)
            
                    

##                if createTiles:
##                    workingFolder = outputBaseFileName + '_' + symbol + '_' + name
##                    if imageSize[0] > tileSize or imageSize[1] > tileSize:
##                        vswInstance.createTiles(inputFile = workingFile, outputFolder = workingFolder, tileSize = tileSize, overlap = 0, suffix = '.png')
##                    else:
##                        outputFolderName = outputBaseFileName + '_' + symbol + '_' + name + '_files'
##                        if os.path.isdir(outputFolderName):
##                            shutil.rmtree(outputFolderName)
##                        os.mkdir(outputFolderName)
##                        os.mkdir(outputFolderName + '/0')
##                        newImage.save(outputFolderName + '/0/0_0.png')
              
        except Exception as e:
            print("Exception: %s" % e)
            success = False
            pass
        
        return success
