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
import CreatePyramidFromTiles as cpft

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

            outputBaseFolder = os.path.split(outputBaseFileName)[0]
            if not outputBaseFolder.endswith('/'):
                outputBaseFolder = outputBaseFolder + '/'
            
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
    

    def insertAllMarkersToNewImagePyramids(self, baseImageFileName, posFileName, tagListFileName, outputBaseFolder, iconFolder, includeNNNN = False, nnnnSymbol = 'redStar', qualityT = 0.5, markerScale = 2):
        success = True
        try:
            baseImageFileName = baseImageFileName.replace('"', '')
            baseImageFileName = baseImageFileName.replace('\\', '/')
            
            vswInstance = vsw.VipsSystemWrapper()
            inputFolder, fileName = os.path.split(baseImageFileName)
            name, ext = os.path.splitext(fileName)
            inputBaseFolder = inputFolder + '/' + name + '_files'

            tileSize, nTilesX, nTilesY, minLevel, maxLevel = vswInstance.getMinMaxLevels(folderName = inputBaseFolder)
            
            posFileName = posFileName.replace('"', '')
            posFileName = posFileName.replace('\\', '/')
            
            tagListFileName = tagListFileName.replace('"', '')
            tagListFileName = tagListFileName.replace('\\', '/')
            
            outputBaseFolder = outputBaseFolder.replace('"', '')
            outputBaseFolder = outputBaseFolder.replace('\\', '/')

            if not outputBaseFolder.endswith('/'):
                outputBaseFolder = outputBaseFolder + '/'
            if os.path.isdir(outputBaseFolder):
                shutil.rmtree(outputBaseFolder)
            os.mkdir(outputBaseFolder)
            
            statsFp = open(outputBaseFolder + 'marker_stats.csv', 'w')
            
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
            usedNames = usedPositions[:,1]

            uniqueLetters = np.unique(usedLetters)
            uniqueNames = np.unique(usedNames)
            #read the taglist file
            taglist = np.genfromtxt(tagListFileName, delimiter=',', dtype=None)
            usedTagIndices = [index for index, value in enumerate(taglist) if value[1] in uniqueNames]
            usedTags = taglist[usedTagIndices]
            cpftInstance = cpft.CreatePyramidFromTiles()

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
                suffix = '.png'
                currentPositionIndices = [index for index, value in enumerate(usedNames) if value == name]
                currentPositions = usedPositions[currentPositionIndices,2:6]
                x = map(int, currentPositions[:,0])
                y = map(int, currentPositions[:,1])
                q = map(float, currentPositions[:,3])
                lenxy = len(x)
                statsFp.write(symbol + "," + name + "," + str(lenxy) + "\n")
                symbolDir = outputBaseFolder + symbol + '_' + name + '_files'
                if os.path.isdir(symbolDir):
                    shutil.rmtree(symbolDir)
                os.mkdir(symbolDir)
                symbolMaxDir = symbolDir + '/' + str(maxLevel)
                os.mkdir(symbolMaxDir)
                for ix in range(0, nTilesX):
                    for iy in range(0, nTilesY):
                        newImage = Image.new("RGB", (tileSize, tileSize))
                        loX = ix * tileSize
                        hiX = (ix + 1) * tileSize
                        loY = iy * tileSize
                        hiY = (iy + 1) * tileSize
                        for i in range(0, lenxy):
                            if (q[i] >= qualityT and x[i] >= loX and x[i] < hiX and y[i] >= loY and y[i] < hiY):
                                pos = (x[i] - halfIconImageSizeScaled[0] - loX, y[i] - halfIconImageSizeScaled[1] - loY)
                                newImage.paste(iconImageScaled, pos)
                        workingFile = symbolMaxDir + '/' + str(ix) + '_' + str(iy) + suffix
                        newImage.save(workingFile)
                        
                cpftInstance.createPyramid(symbolDir, suffix)
                
                        
            if includeNNNN:
                nnnnIndices = [index for index, value in enumerate(allTags) if value == 'NNNN']
                nnnnPositions = positions[nnnnIndices,2:6]
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
                suffix = '.png'
                x = map(int, nnnnPositions[:,0])
                y = map(int, nnnnPositions[:,1])
                q = map(float, nnnnPositions[:,3])
                lenxy = len(x)
                statsFp.write(symbol + "," + name + "," + str(lenxy) + "\n")
                symbolDir = outputBaseFolder + symbol + '_' + name + '_files'
                if os.path.isdir(symbolDir):
                    shutil.rmtree(symbolDir)
                os.mkdir(symbolDir)
                symbolMaxDir = symbolDir + '/' + str(maxLevel)
                os.mkdir(symbolMaxDir)
                for ix in range(0, nTilesX):
                    for iy in range(0, nTilesY):
                        newImage = Image.new("RGB", (tileSize, tileSize))
                        loX = ix * tileSize
                        hiX = (ix + 1) * tileSize
                        loY = iy * tileSize
                        hiY = (iy + 1) * tileSize
                        for i in range(0, lenxy):
                            if (q[i] >= qualityT and x[i] >= loX and x[i] < hiX and y[i] >= loY and y[i] < hiY):
                                pos = (x[i] - halfIconImageSizeScaled[0] - loX, y[i] - halfIconImageSizeScaled[1] - loY)
                                newImage.paste(iconImageScaled, pos)
                        workingFile = symbolMaxDir + '/' + str(ix) + '_' + str(iy) + suffix
                        newImage.save(workingFile)
                        
                cpftInstance.createPyramid(symbolDir, suffix)

            statsFp.close()
              
        except Exception as e:
            print("Exception: %s" % e)
            success = False
            pass
        
        return success
    

    def qualityHeatmap(self, baseImageFileName, posFileName, suffix = '.png'):
        success = True
        try:
            baseImageFileName = baseImageFileName.replace('"', '')
            baseImageFileName = baseImageFileName.replace('\\', '/')
            
            vswInstance = vsw.VipsSystemWrapper()
            inputFolder, fileName = os.path.split(baseImageFileName)
            name, ext = os.path.splitext(fileName)
            inputBaseFolder = inputFolder + '/' + name + '_files'

            tileSize, nTilesX, nTilesY, minLevel, maxLevel = vswInstance.getMinMaxLevels(folderName = inputBaseFolder)
            
            posFileName = posFileName.replace('"', '')
            posFileName = posFileName.replace('\\', '/')
            
            #read the posFile
            positions = np.genfromtxt(posFileName, delimiter=',', dtype=None)
            lenPositions = len(positions[:,0])
            positions = positions[1:lenPositions,2:6] # Remove header

            tileMeans = np.zeros((nTilesX, nTilesY))
            tileNs = np.zeros((nTilesX, nTilesY))
            x = map(int, positions[:,0])
            y = map(int, positions[:,1])
            q = map(float, positions[:,3])
            lenxy = len(x)
            for i in range(0, lenxy):
                tileX = int(x[i]/tileSize)
                tileY = int(y[i]/tileSize)
                tileMeans[tileX, tileY] += q[i]
                tileNs[tileX, tileY] = tileNs[tileX, tileY] + 1

            outputMaxDir = inputBaseFolder + '/' + str(maxLevel)
            if not os.path.isdir(outputMaxDir):
                os.mkdir(outputMaxDir)
                
##            palette = []
##            palette.extend((0, 0, 0))
##            for i in range(1,64):
##                palette.extend((0, i*4, 255))
##            palette.extend((0, 255, 255))
##            for i in range(1,64):
##                palette.extend((0, 255, 255-i*4))
##            palette.extend((0, 255, 0))
##            for i in range(1,64):
##                palette.extend((i*4, 255, 0))
##            palette.extend((255, 255, 0))
##            for i in range(1,64):
##                palette.extend((255, 255-i*4, 0))            
##            palette.extend((255, 0, 0))
                
            palette = []
            palette.extend((0, 0, 0))
            for i in range(1,64):
                palette.extend((0, 0, 0))
            palette.extend((0, 0, 255))
            for i in range(1,16):
                palette.extend((0, i*16, 255))
            palette.extend((0, 255, 255))
            for i in range(1,16):
                palette.extend((0, 255, 255 - i*16))
            palette.extend((0, 255, 0))
            for i in range(1,16):
                palette.extend((i*16, 255, 0))
            palette.extend((255, 255, 0))
            for i in range(1,16):
                palette.extend((255, 255 - i*16, 0))
            palette.extend((255, 0, 0))
            for i in range(1,128):
                palette.extend((255, i*2, i*2))
            palette.extend((255, 255, 255))
            
            im = Image.new("P", (tileSize, tileSize))
            im.putpalette(palette)
            pixelMap = im.load()
            border = int(tileSize / 20)
            halfBorder = int(border / 2)
            step = (tileSize - border) / 256.0
            for i in range(halfBorder + 1, tileSize - halfBorder):
                for j in range(0, border + 1):
                    pixelMap[i,j] = int((i-halfBorder)/step)
            for ix in range(0, nTilesX):
                for iy in range(0, nTilesY):
                    if (tileNs[ix, iy] > 0):
                        tileMeans[ix, iy] = int((255.0 * tileMeans[ix, iy])/tileNs[ix, iy] + 0.5)
                    for i in range(halfBorder + 1, tileSize - halfBorder):
                        for j in range(border + 1, tileSize - border):
                            pixelMap[i,j] = tileMeans[ix, iy]
                    workingFile = outputMaxDir + '/' + str(ix) + '_' + str(iy) + suffix
                    im.save(workingFile)
              
        except Exception as e:
            print("Exception: %s" % e)
            success = False
            pass
        
        return success
