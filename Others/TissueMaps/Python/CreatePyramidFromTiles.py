# File: CreatePyramidFromTiles.py
# Generate resolution pyramid from a set of tiles
# Petter Ranefall 2013

# Import neccessary modules:
import os
import shutil
from PIL import Image
import subprocess
from ..VipsSystemWrapper import VipsSystemWrapper as vsw

class CreatePyramidFromTiles(object):

    def __init__(self):
        pass

    def getTilesInfo(self, inputFolderName, suffix = ".png"):
        if inputFolderName.endswith('/'):
            inputFolderName = os.path.split(inputFolderName)[0]
        folderName = os.path.split(inputFolderName)[1]
        
        baseFolder = inputFolderName + '/'
        maxX = -1
        maxY = -1
        
        image0_0 = baseFolder + '0_0' + suffix
        tileSize = 0
        if os.path.isfile(image0_0):
            pilIm1 = Image.open(image0_0)
            tileSize = max(pilIm1.size[0], pilIm1.size[1])
                
        for fileName in os.listdir(baseFolder):
            if os.path.isfile(baseFolder + fileName):
                if fileName.endswith(suffix):
                    x_y = fileName.split(suffix)[0]
                    xstr, ystr = x_y.split('_')
                    x = int(xstr)
                    y = int(ystr)
                    if x > maxX:
                        maxX = x
                    if y > maxY:
                        maxY = y

        nTilesX = maxX + 1
        nTilesY = maxY + 1
        return tileSize, nTilesX, nTilesY

    def createPyramid(self, baseFolder, suffix = '.png'):
        baseFolder = baseFolder.replace('\\', '/')
        baseFolder = baseFolder.replace('"', '')
        if not baseFolder.endswith('/'):
            baseFolder = baseFolder + '/'
        vswInstance = vsw.VipsSystemWrapper()
        tileSize, nTilesX, nTilesY, minLevel, maxLevel = vswInstance.getMinMaxLevels(folderName = baseFolder, suffix = suffix)
        inputFolder = baseFolder + str(maxLevel)
        folderName = str(maxLevel)
        firstTileSize = 0
        nTilesXY = (nTilesX, nTilesY)
        firstNTilesY = 0
        parentFolder = baseFolder
        lut = [255 if v > 0 else 0 for v in range(256)]
        for y in range(0, nTilesY):
            for x in range(0, nTilesX):
                fileName = inputFolder + '/' + str(x)+ '_' + str(y) + suffix
                try:
                    pilIm = Image.open(fileName)
                    mask = pilIm.convert("L") 
                    mask = mask.point(lut) 
                    pilIm.putalpha(mask)
                    pilIm.save(fileName)
                except Exception as e:
                    print("Exception: %s" % e)
                    pass
        
        for level in range(maxLevel - 1, - 1, -1):
            folder = parentFolder + str(level)
            if not os.path.isdir(folder):
                os.mkdir(folder)
            prevLevel = level + 1
            prevFolder = parentFolder + str(prevLevel)
            tileSize, nTilesX, nTilesY = self.getTilesInfo(prevFolder, suffix)
            if prevLevel == maxLevel:
                firstTileSize = tileSize
                nTilesXY = (nTilesX, nTilesY)
                
            for y in range(1, nTilesY, 2):
                for x in range(1, nTilesX, 2):
                    try:
                        prevFileName00 = str(x-1)+ '_' + str(y-1) + suffix
                        prevFileName01 = str(x-1)+ '_' + str(y) + suffix
                        prevFileName10 = str(x)+ '_' + str(y-1) + suffix
                        prevFileName11 = str(x)+ '_' + str(y) + suffix
                        pilIm00 = Image.open(prevFolder + '/' + prevFileName00)
                        pilIm01 = Image.open(prevFolder + '/' + prevFileName01)
                        pilIm10 = Image.open(prevFolder + '/' + prevFileName10)
                        pilIm11 = Image.open(prevFolder + '/' + prevFileName11)
                        halfSize00 = (int((1+pilIm00.size[0])/2), int((1+pilIm00.size[1])/2))
                        halfSize01 = (int((1+pilIm01.size[0])/2), int((1+pilIm01.size[1])/2))
                        halfSize10 = (int((1+pilIm10.size[0])/2), int((1+pilIm10.size[1])/2))
                        halfSize11 = (int((1+pilIm11.size[0])/2), int((1+pilIm11.size[1])/2))
                        pilIm00 = pilIm00.resize(halfSize00, Image.NEAREST)
                        pilIm01 = pilIm01.resize(halfSize01, Image.NEAREST)
                        pilIm10 = pilIm10.resize(halfSize10, Image.NEAREST)
                        pilIm11 = pilIm11.resize(halfSize11, Image.NEAREST)
                        newSize = (halfSize00[0] + halfSize10[0], halfSize01[1] + halfSize11[1])
                        newImage = Image.new(pilIm00.mode, newSize)
                        box00 = (0, 0)
                        box01 = (0, halfSize00[1])
                        box10 = (halfSize00[0], 0)
                        box11 = (halfSize00[0], halfSize00[1])
                        newImage.paste(pilIm00, box00)
                        newImage.paste(pilIm01, box01)
                        newImage.paste(pilIm10, box10)
                        newImage.paste(pilIm11, box11)
                        newFileName = str(int(x / 2)) + '_' + str(int(y / 2)) + suffix
                        newImage.save(folder + '/' + newFileName)
                    except Exception as e:
                        print("Exception: %s" % e)
                        pass
                if nTilesX % 2 == 1:
                    try:
                        x = nTilesX
                        prevFileName00 = str(x-1)+ '_' + str(y-1) + suffix
                        prevFileName01 = str(x-1)+ '_' + str(y) + suffix
                        pilIm00 = Image.open(prevFolder + '/' + prevFileName00)
                        pilIm01 = Image.open(prevFolder + '/' + prevFileName01)
                        halfSize00 = (int((1+pilIm00.size[0])/2), int((1+pilIm00.size[1])/2))
                        halfSize01 = (int((1+pilIm01.size[0])/2), int((1+pilIm01.size[1])/2))
                        pilIm00 = pilIm00.resize(halfSize00, Image.NEAREST)
                        pilIm01 = pilIm01.resize(halfSize01, Image.NEAREST)
                        newSize = (halfSize00[0], halfSize00[1]+halfSize01[1])
                        newImage = Image.new(pilIm00.mode, newSize)
                        box00 = (0, 0)
                        box01 = (0, halfSize00[1])
                        newImage.paste(pilIm00, box00)
                        newImage.paste(pilIm01, box01)
                        newFileName = str(int(x / 2)) + '_' + str(int(y / 2)) + suffix
                        newImage.save(folder + '/' + newFileName)
                    except Exception as e:
                        print("Exception: %s" % e)
                        pass
            if nTilesY % 2 == 1:
                y = nTilesY
                for x in range(1, nTilesX, 2):
                    try:
                        prevFileName00 = str(x-1)+ '_' + str(y-1) + suffix
                        prevFileName10 = str(x)+ '_' + str(y-1) + suffix
                        pilIm00 = Image.open(prevFolder + '/' + prevFileName00)
                        pilIm10 = Image.open(prevFolder + '/' + prevFileName10)
                        halfSize00 = (int((1+pilIm00.size[0])/2), int((1+pilIm00.size[1])/2))
                        halfSize10 = (int((1+pilIm10.size[0])/2), int((1+pilIm10.size[1])/2))
                        pilIm00 = pilIm00.resize(halfSize00, Image.NEAREST)
                        pilIm10 = pilIm10.resize(halfSize10, Image.NEAREST)
                        newSize = (halfSize00[0] + halfSize10[0], halfSize00[1])
                        newImage = Image.new(pilIm00.mode, newSize)
                        box00 = (0, 0)
                        box10 = (halfSize00[0], 0)
                        newImage.paste(pilIm00, box00)
                        newImage.paste(pilIm10, box10)
                        newFileName = str(int(x / 2)) + '_' + str(int(y / 2)) + suffix
                        newImage.save(folder + '/' + newFileName)
                    except Exception as e:
                        print("Exception: %s" % e)
                        pass
                if nTilesX % 2 == 1:
                    try:
                        x = nTilesX
                        prevFileName00 = str(x-1)+ '_' + str(y-1) + suffix
                        pilIm00 = Image.open(prevFolder + '/' + prevFileName00)
                        halfSize00 = (int((1+pilIm00.size[0])/2), int((1+pilIm00.size[1])/2))
                        pilIm00 = pilIm00.resize(halfSize00, Image.NEAREST)
                        newFileName = str(int(x / 2)) + '_' + str(int(y / 2)) + suffix
                        pilIm00.save(folder + '/' + newFileName)
                    except Exception as e:
                        print("Exception: %s" % e)
                        pass            
        return nTilesXY, firstTileSize, maxLevel
