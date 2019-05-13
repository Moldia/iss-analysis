# File: MakeTilesAndCreateCsv.py
# Regisering large slide scanner images and create csv file
# Petter Ranefall 2013

# Import neccessary modules:
import sys
import numpy as np
import scipy.ndimage as nd
import time
from PIL import Image
import os
from . import TissueMapsRegistrationVips as ssr
from ..ImageReader import ImageReader as imageReader

class MakeTilesAndCreateCsv(object):
    """
    Read, register and tile large images from e.g. slide scanner, and create a .csv file as input for CellProfiler
    """
    
    def __init__(self):
        pass

    def run(self, inputFiles, \
            createCsv = True, \
            csvFileName = None, \
            additionalFiles = None, \
            align = True, \
            createTiles = True, \
            suffix = ".png", \
            outputTileSize = 1024, \
            flipToFit = False, \
            verbose = False, \
            inputTileSize = 512,
            refBase = 1,
            createRGImage = False):
        """
        Input parameters:
        csvFileName: name (including path) of output csv file
        inputFiles: string with comma separated file names (including path) to images that should be aligned AND/OR tiled.
        There should be 6 images for each hyb, in the order A, C, T, G, General_blob, Nuclei. The first hyb is considered as the reference to which the others are aligned. (NEEDED)
        additionalFiles: extra files to be tiled, but not used for analysis. There must be the same number of extra files for all hybs.
        align: True if automatically align reference and floating
        createTiles: True if create tiles for input files
        inputTileSize: size of tiles in the input images. This is used when extracting the largest non-empty square. (DEFAULT = 512)
        createTiles: True if tiles should be saved. (DEFAULT = True)
        outputTileSize: Target tile size. (DEFAULT = 1024)
        flipToFit: True if flip to same orientation as the screen. (DEFAULT = False)
        verbose: True if progress strings should be written. (DEFAULT = True)
        refBase: Which of the bases that should be the reference in alignment. First index is 1. (DEFAULT = 1)

        Output:
        csv file
        tiles will be saved in a subfolder to where the corresponding image lies, i.e. inputFolder + '/' + name + '_files/'
        full size aligned images will be saved in a subfolder to where the corresponding image lies, i.e. inputFolder + '/Transformed/'  

        Prerequisites:
        Pyelastix\pyelastix.py + __init__.py
        Elastix must be installed
        Vips must be installed
        OpenSeadragon must be installed
        """

        refBase0 = refBase - 1;
        start = time.clock()
        prev = start

        inputFileList = inputFiles.split(',')
        nInputFiles = len(inputFileList)
        nHybs = nInputFiles / 6

        assert nInputFiles > 0, "Input file list must not be empty."
        
        if createCsv:
            assert csvFileName is not None, "Output Csv file name must not be empty"
            assert nInputFiles % 6 == 0, "Wrong number of input files. There must be exactly 6 images per hyb. nInput: %s nHybs %s %s" % (str(nInputFiles), str(nHybs), inputFiles)
        
        nAdditionalFilesPerHyb = 0
        if additionalFiles is not None:
            additionalFiles = additionalFiles.strip()
            if len(additionalFiles) > 0:
                additionalFilesList = additionalFiles.split(',')
                nAdditionalFiles = len(additionalFilesList)
                assert nAdditionalFiles % nHybs == 0, "Wrong number of additional files. There must be the same number of additional images per hyb. nAdditional: %s nHybs %s %s" % (str(nAdditionalFiles), str(nHybs), additionalFiles)
                nAdditionalFilesPerHyb = nAdditionalFiles / nHybs  

        referenceFileName = None
        nFirstInputFiles = nInputFiles
        if nInputFiles >= 6:
            nFirstInputFiles = 6
            referenceFileName = inputFileList[refBase0 * nFirstInputFiles + 5]
            referenceFileName = referenceFileName.strip()
            
        ssrInstance = ssr.TissueMapsRegistrationVips()

        #Create tiles from first hyb
        currentFileList = inputFileList[refBase0 * nFirstInputFiles : refBase0 * nFirstInputFiles + nFirstInputFiles]
        for j in range(0, nAdditionalFilesPerHyb):
            currentFile = additionalFilesList[refBase0 * nAdditionalFilesPerHyb + j]
            currentFile = currentFile.strip()
            if len(currentFile) > 0:
                currentFileList.append(currentFile)
            
        nTilesXY, outputTileSize, maxLevel = ssrInstance.run(currentFileList, \
                                                         align = False, 
                                                         referenceFileName = referenceFileName, \
                                                         createTiles = createTiles, \
                                                         suffix = suffix, \
                                                         outputTileSize = outputTileSize, \
                                                         flipToFit = flipToFit, \
                                                         verbose = verbose)
        
        if createRGImage and nInputFiles >= 6:
            imageReaderInstance = imageReader.ImageReader()
            path, ext = os.path.splitext(inputFileList[refBase0 * nFirstInputFiles + 5])
            outputRGName = path + "_RG" + ext
            imageReaderInstance.mergeTilesToRG(inputFileList[refBase0 * nFirstInputFiles + 5], inputFileList[refBase0 * nFirstInputFiles + 4], outputRGName)
            currentFileList.append(outputRGName)
            
        refNTilesXY = nTilesXY
        refTileSize = outputTileSize

        #Create tiles from the rest of the hybs
        for i in range(0, nInputFiles, nFirstInputFiles):
            base = i / nFirstInputFiles
            if (base != refBase0):
                currentFileList = inputFileList[i : i + 6]
                for j in range(0, nAdditionalFilesPerHyb):
                    currentFile = additionalFilesList[base * nAdditionalFilesPerHyb + j]
                    currentFile = currentFile.strip()
                    if len(currentFile) > 0:
                        currentFileList.append(currentFile)
                        
                nTilesXY, outputTileSize, maxLevel = ssrInstance.run(currentFileList, \
                                                                     align = align,
                                                                     referenceFileName = referenceFileName, floatingFileName = inputFileList[i+5], \
                                                                     createTiles = createTiles, \
                                                                     suffix = suffix, \
                                                                     outputTileSize = outputTileSize, \
                                                                     flipToFit = flipToFit, \
                                                                     verbose = verbose)
                
                if createRGImage and nInputFiles >= 6:
                    path, ext = os.path.splitext(inputFileList[i * nFirstInputFiles + 5])
                    outputRGName = path + "_RG" + ext
                    imageReaderInstance.mergeTilesToRG(inputFileList[i * nFirstInputFiles + 5], inputFileList[i * nFirstInputFiles + 4], outputRGName)
                    currentFileList.append(outputRGName)
                    
            ##assert outputTileSize == refTileSize and nTilesXY == refNTilesXY, "All images must have the same size."



        if createCsv:
            #Create CSV file
            csvFolder = os.path.dirname(csvFileName)
            if not os.path.exists(csvFolder):
                os.makedirs(csvFolder)
            csvFile = open(csvFileName, "w")

            headings = 'Metadata_position,Metadata_Tile_xPos,Metadata_Tile_yPos,Metadata_HybStep'
            headings += ',Image_PathName_A,Image_FileName_A'
            headings += ',Image_PathName_C,Image_FileName_C'
            headings += ',Image_PathName_G,Image_FileName_G'
            headings += ',Image_PathName_T,Image_FileName_T'
            headings += ',Image_PathName_General_blob,Image_FileName_General_blob'
            headings += ',Image_PathName_Spec_blob,Image_FileName_Spec_blob'
            headings += ',Image_PathName_Nuclei,Image_FileName_Nuclei'
            headings += ',Image_PathName_Spec_Nuclei,Image_FileName_Spec_Nuclei'
            
            csvFile.write("%s" % headings)

            metadataPosition = 1
            if verbose:
                print("outputTileSize %s" % str(outputTileSize))
            
            for y in range(0, nTilesXY[1]):
                for x in range(0, nTilesXY[0]):
                    #Only use tiles with same size as outputTileSize, i.e. not tiles at the bottom or right borders
                    tileSizeOK = True
                    inputFolder, fileName = os.path.split(inputFileList[refBase0 * nFirstInputFiles + 5])
                    name, ext = os.path.splitext(fileName)
                    tilesOutputFolder = inputFolder + '/' + name + '_files/' + str(maxLevel) + '/'
                    tilesOutputFolder = tilesOutputFolder.strip()
                    tilesOutputFileName = "%s_%s%s" % (str(x), str(y), suffix)
                    tilesOutputFileName = tilesOutputFileName.strip()
                    if os.path.isfile(tilesOutputFolder + tilesOutputFileName):
                        pilIm1 = Image.open(tilesOutputFolder + tilesOutputFileName)
                        minTileSize = min(pilIm1.size[0], pilIm1.size[1])
                        if not minTileSize == outputTileSize:
                            tileSizeOK = False
                    else:
                        tileSizeOK = False
                    if tileSizeOK:
                        for h in range(0, nHybs):
                            line = str(metadataPosition)
                            line += ',' + str(x * outputTileSize) + ',' + str(y * outputTileSize)
                            line += ',hyb' + str(h)
                            for s in range(0, 4): #A,C,T,G,
                                inputFolder, fileName = os.path.split(inputFileList[h * 6 + s])
                                name, ext = os.path.splitext(fileName)
                                inputFolder = inputFolder.replace('\\','/')
                                tilesOutputFolder = inputFolder + '/' + name + '_files/' + str(maxLevel) + '/'
                                tilesOutputFolder = tilesOutputFolder.strip()
                                tilesOutputFileName = "%s_%s%s" % (str(x), str(y), suffix)
                                line += ',' + tilesOutputFolder
                                line += ',' + tilesOutputFileName
                            for s in range(4, 6): #General_blob, Spec_blob, Nuclei from first hyb, spec_nuclei
                                #General
                                inputFolder, fileName = os.path.split(inputFileList[refBase0 * 6 + s]) #Note: Not h * 6 + s
                                name, ext = os.path.splitext(fileName)
                                tilesOutputFolder = inputFolder + '/' + name + '_files/' + str(maxLevel) + '/'
                                tilesOutputFolder = tilesOutputFolder.strip()
                                tilesOutputFileName = "%s_%s%s" % (str(x), str(y), suffix)
                                tilesOutputFileName = tilesOutputFileName.strip()
                                line += ',' + tilesOutputFolder
                                line += ',' + tilesOutputFileName
                                #Spec
                                inputFolder, fileName = os.path.split(inputFileList[h * 6 + s]) #Note: h * 6 + s
                                name, ext = os.path.splitext(fileName)
                                tilesOutputFolder = inputFolder + '/' + name + '_files/' + str(maxLevel) + '/'
                                tilesOutputFolder = tilesOutputFolder.strip()
                                line += ',' + tilesOutputFolder
                                line += ',' + tilesOutputFileName 

                            csvFile.write("\n%s" % line)
                    metadataPosition += 1                 
            
            csvFile.close()
        
        # Done!
        if verbose:
            stop = time.clock()
            print("Total time %s s" % str(stop-start))
        pass
