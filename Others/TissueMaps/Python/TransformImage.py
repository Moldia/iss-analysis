# File: TransformImage.py
# Transform an image
# Petter Ranefall 2013

# Import neccessary modules:
import sys
import os
import numpy as np
import scipy.ndimage as nd
import time
from ..Tifffile import tifffile
from ..VipsSystemWrapper import VipsSystemWrapper as vsw
from ..ImageReader import ImageReader as imageReader
from PIL import Image
import ctypes

class TransformImage(object):
    """
    Read, register and tile large images from e.g. slide scanner
    """
    
    def __init__(self):
        pass
    
    def run(self, inputFileName, outputFileName, \
            referencePoints, floatingPoints, filter = 'NEAREST', \
            createTiles = True, \
            tileSize = 1024, \
            tilesSuffix = '.png', \
            verbose = True):
        
        """
        Affine transform of an image. The data points parameters should be strings on the form x1,y1,x2,y2,x3,y3,...
        """

        start = time.clock()
        prev = start
      
        imageReaderInstance = imageReader.ImageReader()
        
        if verbose:
            print("Start opening images") 
        # Read image
        print(inputFileName)
        im1, pilIm1 = imageReaderInstance.openAsNumpyArrayAndPilImage(inputFileName)
        if verbose:
            stop = time.clock()
            print("Image1 opened. Time %s s" % str(stop-prev))
            print(inputFileName)
            prev = stop
        
        # Apply the transformation
        if filter == 'BILINEAR':
            resample = Image.BILINEAR
        elif filter == 'BICUBIC':
            resample = Image.BICUBIC
        else:
            resample = Image.NEAREST

        try:
            referencePoints = referencePoints.replace('(', '')
            referencePoints = referencePoints.replace(')', '')
            refPointsList = referencePoints.split(',')
            floatingPoints = floatingPoints.replace('(', '')
            floatingPoints = floatingPoints.replace(')', '')
            floatingPointsList = floatingPoints.split(',')
            lr = len(refPointsList)
            lf = len(floatingPointsList)
            assert lr == lf, "There must be the same number of points"
            assert lr >= 6, "There must be at least three points"
            assert lr % 2 == 0, "Missing y-coordinate"
            aList = [[float(refPointsList[0]), float(refPointsList[1]), 1]]
            bList = [[float(floatingPointsList[0]), float(floatingPointsList[1]), 1]]
            for i in range(2, lr, 2):
                aList.append([float(refPointsList[i]), float(refPointsList[i+1]), 1])
                bList.append([float(floatingPointsList[i]), float(floatingPointsList[i+1]), 1])
            a = np.asarray(aList)
            b = np.asarray(bList)
            if verbose:
                print('a: %s' % str(a))
                print('b: %s' % str(b))
            t=np.linalg.lstsq(a, b)[0]
            tData = (t[0][0], t[1][0], t[2][0], t[0][1], t[1][1], t[2][1])
            if verbose:
                print('t: %s' % str(t))
                print('tData: %s' % str(tData))
            transformed = pilIm1.transform(pilIm1.size, Image.AFFINE, tData, resample=resample)
            
            outputFolder, fileName = os.path.split(outputFileName)
            if not os.path.isdir(outputFolder):
                os.mkdir(outputFolder)
            transformed.save(outputFileName)
            if createTiles:
                vswInstance = vsw.VipsSystemWrapper()
                inputFolder, fileName = os.path.split(inputFileName)
                parentFolder, baseName = os.path.split(inputFolder)
                name, ext = os.path.splitext(fileName)
                if baseName == 'Transformed':
                    inputFolder = parentFolder
                tilesFolder = inputFolder + '/' + name
                outputFilesFolder = tilesFolder + "_files"
                vswInstance.createTiles(inputFile = outputFileName, outputFolder = tilesFolder, tileSize = tileSize, overlap = 0, suffix = tilesSuffix)
                    
        except Exception as ex:
            print("Failed: %s" % str(ex))

        # Done!
        if verbose:
            stop = time.clock()
            print("Total time %s s" % str(stop-start))

    
    def reset(self, inputFileName, outputFileName, \
            createTiles = True, \
            tileSize = 1024, \
            tilesSuffix = '.png', \
            verbose = True):
        
        """
        Reset the transformation by copying from input to Transformed.
        """
        start = time.clock()
        prev = start
      
        try:
            imageReaderInstance = imageReader.ImageReader()
            
            if verbose:
                print("Start opening images") 
            # Read image
            print(inputFileName)
            pilIm1 = imageReaderInstance.openAsNumpyArrayAndPilImage(inputFileName)[1]
            if verbose:
                stop = time.clock()
                print("Image1 opened. Time %s s" % str(stop-prev))
                print(inputFileName)
                prev = stop
                
            outputFolder, fileName = os.path.split(outputFileName)
            if not os.path.isdir(outputFolder):
                os.mkdir(outputFolder)
            pilIm1.save(outputFileName)
            if createTiles:
                vswInstance = vsw.VipsSystemWrapper()
                inputFolder, fileName = os.path.split(inputFileName)
                parentFolder, baseName = os.path.split(inputFolder)
                name, ext = os.path.splitext(fileName)
                if baseName == 'Transformed':
                    inputFolder = parentFolder
                tilesFolder = inputFolder + '/' + name
                outputFilesFolder = tilesFolder + "_files"
                vswInstance.createTiles(inputFile = outputFileName, outputFolder = tilesFolder, tileSize = tileSize, overlap = 0, suffix = tilesSuffix)
                    
        except Exception as ex:
            print("Failed: %s" % str(ex))

        # Done!
        if verbose:
            stop = time.clock()
            print("Total time %s s" % str(stop-start))
    
    def __scaleRotateTranslate(self, pilImage, angle, rotationCenter = None, translation = None, scale = None):
        """
        Apply transformation using PIL
        """
        
        if rotationCenter is None:
            y, x = pilImage.size[1]/2, pilImage.size[0]/2
        else:
            y, x = rotationCenter
        if translation is None:
            dy = dx = 0
        else:
            dy, dx = translation
        angle = -angle/180.0*np.pi
        sx=sy=1.0
        if scale:
            (sx,sy) = scale
        cosine = np.cos(angle)
        sine = np.sin(angle)
        a = cosine/sx
        b = sine/sx
        c = x - x * a - y * b + dx
        d = -sine/sy
        e = cosine/sy
        f = y - y * a + x * b + dy
        return pilImage.transform(pilImage.size, Image.AFFINE, (a,b,c,d,e,f), resample=Image.NEAREST)
