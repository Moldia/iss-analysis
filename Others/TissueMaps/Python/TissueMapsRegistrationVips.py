# File: TissueMapsRegistrationVips.py
# Regisering large slide scanner images
# Petter Ranefall 2013

# Import neccessary modules:
import sys
import numpy as np
import scipy.ndimage as nd
import time
from ..Tifffile import tifffile
from PIL import Image
import os
from ..Pyelastix import pyelastix
from ..VipsSystemWrapper import VipsSystemWrapper as vsw
from ..ImageReader import ImageReader as imageReader
import ctypes
import Tkinter

class TissueMapsRegistrationVips(object):
    """
    Read, register and tile large images from e.g. slide scanner
    """
    
    def __init__(self):
        pass
    
    def run(self, inputFileList, \
            align = True, referenceFileName = None, floatingFileName = None, \
            inputTileSize = 512, \
            createTiles = False, \
            suffix = ".png", \
            outputTileSize = 1024, 
            flipToFit = False, \
            verbose = True):
        
        """
        Run the following sequence of operations:
        1) Find the largest non-empty square in the reference image and extract corresponding subimages
        2) Compute rigid registration using Elastix
        3) Apply transformation to floating image
        4) Save transformed image
        5) Tile reference and floating images and save

        Input parameters:
        inputFileList: list of file names (including path) to images that should be aligned AND/OR tiled. (NEEDED)
        align: True if automatically align reference and floating
        referenceFileName: file name (including path) of the input reference image. (NEEDED if align)
        floatingFileName: file name (including path) of the input floating image. (NEEDED if align)
        createTiles: True if create tiles for input files
        inputTileSize: size of tiles in the input images. This is used when extracting the largest non-empty square. (DEFAULT = 512)
        createTiles: True if tiles should be saved. (DEFAULT = True)
        outputTileSize: Target tile size, max = 1024. (DEFAULT = 1024)
        flipToFit: True if flip to same orientation as the screen. (DEFAULT = False)
        verbose: True if progress strings should be written. (DEFAULT = True)

        Output:
        tiles will be saved in a subfolder to where the corresponding image lies, i.e. inputFolder + '/' + name + '_files/'
        full size aligned images will be saved in a subfolder to where the corresponding image lies, i.e. inputFolder + '/Transformed/'
        
        Prerequisites:
        Tifffile\tifffile.py + __init__.py
        Pyelastix\pyelastix.py + __init__.py
        Elastix must be installed
        Vips must be installed
        OpenSeadragon must be installed
        """
        
        imageReaderInstance = imageReader.ImageReader()
        
        if verbose:
            print("Reference file: %s" % referenceFileName)
            print("Floating file: %s" % floatingFileName)

        if referenceFileName is not None:
            referenceFileName = referenceFileName.strip()
            referenceFileName = referenceFileName.replace('\\', '/')
        if floatingFileName is not None:
            floatingFileName = floatingFileName.strip()
            floatingFileName = floatingFileName.replace('\\', '/')

        assert align is False or (referenceFileName is not None and os.path.isfile(referenceFileName)), "Missing reference file"
        assert align is False or (floatingFileName is not None and os.path.isfile(floatingFileName)), "Missing floating file name"
        start = time.clock()
        prev = start
      
        flip = False 
        root = Tkinter.Tk()
        root.withdraw()
        screensize = root.winfo_screenwidth(), root.winfo_screenheight()
        if verbose:
            print("screensize %s" % (str(screensize)))

        if referenceFileName is not None:
            if verbose:
                print("Start opening images") 
            # Read image 1
            referenceFileName = referenceFileName.strip()
            im1, pilIm1 = imageReaderInstance.openAsNumpyArrayAndPilImage(referenceFileName)
            if verbose:
                stop = time.clock()
                print("Image1 opened. Time %s s" % str(stop-prev))
                print(referenceFileName)
                print("im1.shape %s pilIm1.size %s" % (str(im1.shape), str(pilIm1.size)))
                prev = stop              

            nTiles = (int(np.trunc(im1.shape[0]/outputTileSize)), int(np.trunc(im1.shape[1]/outputTileSize)))
            newShape = (outputTileSize * nTiles[0], outputTileSize * nTiles[1])
            if verbose:
                print("nTiles %s newShape %s" % (str(nTiles), str(newShape)))
            
            if newShape != im1.shape:
                #Crop to new shape
                pilIm1 = pilIm1.crop((0, 0, newShape[1], newShape[0]))
                im1 = np.array(pilIm1)
                
            if flipToFit:
                if (screensize[0] > screensize[1] and pilIm1.size[0] < pilIm1.size[1]) or \
                   (screensize[0] < screensize[1] and pilIm1.size[0] > pilIm1.size[1]):
                    flip = True               
                    if verbose:
                        print("Flip")

            if flip:
                pilIm1 = pilIm1.transpose(Image.ROTATE_90)
                im1 = np.array(pilIm1)
                
        if align:            
            # Read image 2
            floatingFileName = floatingFileName.strip() 
            im2, pilIm2 = imageReaderInstance.openAsNumpyArrayAndPilImage(floatingFileName)

            if flip:
                pilIm2 = pilIm2.transpose(Image.ROTATE_90)
                im2 = np.array(pilIm2)
                
            if im2.shape != im1.shape:
                #Crop/Pad to same size as im1
                pilIm2 = pilIm2.crop((0, 0, im1.shape[1], im1.shape[0]))
                im2 = np.array(pilIm2)
            
            if verbose:
                stop = time.clock()
                print("Image2 opened. Time %s s" % str(stop-prev))
                prev = stop
            
            if verbose:
                print("Start extracting subimages")
            # Get parameters for largest non-empty square
            corner, size = self.__largestNonEmptySquare(im1, inputTileSize)
            oppositeCorner = (corner[0] + size[0] - 1, corner[1] + size[1] - 1)

            # Adjust size if image 2 is too small
            modifiedSize = False
            if oppositeCorner[0] >= im2.shape[0]:
                oppositeCorner = (im2.shape[0] - 1, oppositeCorner[1])
                modifiedSize = True
            if oppositeCorner[1] >= im2.shape[1]:
                oppositeCorner = (oppositeCorner[0], im2.shape[1] - 1)
                modifiedSize = True
            if modifiedSize:
                size = (1 + oppositeCorner[0] - corner[0], 1 + oppositeCorner[1] - corner[1])
                minSubImageSize = np.min(size[0], size[1])
                size = (minSubImageSize, minSubImageSize)

            # Extract subimages
            subIm1 = im1[corner[0] : corner[0] + size[0], corner[1] : corner[1] + size[1]]
            subIm2 = im2[corner[0] : corner[0] + size[0], corner[1] : corner[1] + size[1]]
            if verbose:
                stop = time.clock()
                print("Subimages extracted. Time %s s" % str(stop-prev))
                prev = stop
                
            if verbose:
                print("Start computing transformation")
            # Compute transformation
            rotation, rotationCenter, translation = self.__computeRigidRegistration(subIm1, subIm2)
            
            # Adjust rotation center
            rotationCenter = (corner[0] + rotationCenter[0], corner[1] + rotationCenter[1])
            
            # Save parameters to file in output folder
            
            outputFolder = os.path.split(floatingFileName)[0] + '/Transformed/'
            outputFolder = outputFolder.replace('\\', '/')
            if not os.path.exists(outputFolder):
                os.makedirs(outputFolder)
            parameterFile = open(outputFolder + "transformationParameters.txt", "w")
            parameterFile.write("ReferenceFile: %s\n" % referenceFileName)
            parameterFile.write("FloatingFile: %s\n" % floatingFileName)
            parameterFile.write("Rotation: %s (deg)\n" % str(rotation))
            parameterFile.write("RotationCenter: %s (pix)\n" % str(rotationCenter))
            parameterFile.write("Translation: %s (pix)\n" % str(translation))
            parameterFile.close();
            
            if verbose:
                stop = time.clock()
                print("Transformation computed. Time %s s" % str(stop-prev))
                prev = stop
            
        margin = (200, 130)
        windowWidth = screensize[0] - margin[0]
        windowHeight = screensize[1] - margin[1]

        nTilesXY = (0, 0)
        maxLevel = 0
        vswInstance = vsw.VipsSystemWrapper()
        
        # Apply the transformation
        first = True
        for filePath in inputFileList:
            workingFile = filePath
            filePath = filePath.strip()
            inputFolder, fileName = os.path.split(filePath)
            inputFolder = inputFolder.replace('\\', '/')
            name, ext = os.path.splitext(fileName)
            if verbose:
                print('filePath: %s' % filePath)
                print('inputFolder: %s' % inputFolder)
                print('fileName: %s' % fileName)
                
            try:
                im, pilIm = imageReaderInstance.openAsNumpyArrayAndPilImage(filePath)

                if align is False and flipToFit:
                    if (screensize[0] > screensize[1] and pilIm.size[0] < pilIm.size[1]) or \
                       (screensize[0] < screensize[1] and pilIm.size[0] > pilIm.size[1]):
                        flip = True
                        
                if flip:
                    pilIm = pilIm.transpose(Image.ROTATE_90)
                    im = np.array(pilIm)
                    if verbose:
                        stop = time.clock()
                        print("Flipped. Time %s s" % str(stop-prev))
                        prev = stop
                        
                if referenceFileName is not None: 
                    if pilIm.size != pilIm1.size:
                        #Crop/Pad to same size as im1
                        pilIm = pilIm.crop((0, 0, pilIm1.size[0], pilIm1.size[1]))
                else:           
                    nTiles = (int(np.trunc(im.shape[0]/outputTileSize)), int(np.trunc(im.shape[1]/outputTileSize)))
                    newShape = (outputTileSize * nTiles[0], outputTileSize * nTiles[1])
                    if verbose:
                        print("nTiles %s newShape %s" % (str(nTiles), str(newShape)))
                    
                    if newShape != im.shape:
                        #Crop to new shape
                        pilIm = pilIm.crop((0, 0, newShape[1], newShape[0]))
                        im = np.array(pilIm)
                        outputFolder = inputFolder + '/Transformed/'            
                        if not os.path.exists(outputFolder):
                            os.makedirs(outputFolder)
                        workingFile = outputFolder + fileName
                        pilIm.save(workingFile)
  
                if align:           
                    pilIm = self.__scaleRotateTranslate(pilIm, rotation, rotationCenter = rotationCenter, translation = translation )
                    if verbose:
                        stop = time.clock()
                        print("Transformation applied. Time %s s" % str(stop-prev))
                        prev = stop            
                    if verbose:
                        print("Start save transformed image")
                    # Save the transformed image
                    outputFolder = inputFolder + '/Transformed/'            
                    if not os.path.exists(outputFolder):
                        os.makedirs(outputFolder)
                    workingFile = outputFolder + fileName
                    pilIm.save(workingFile)
                    if verbose:
                        stop = time.clock()
                        print("Transformed Image saved. Time %s s" % str(stop-prev))
                        prev = stop
                else:
                    outputFolder = inputFolder + '/Transformed/'           
                    if not os.path.exists(outputFolder):
                        os.makedirs(outputFolder)
                    workingFile = outputFolder + fileName
                    pilIm.save(workingFile)
                    
                tilesFolder = inputFolder + '/' + name
                outputFilesFolder = tilesFolder + "_files"
                if createTiles:
                    outputTileSize, nTilesWidth, nTilesHeight, minLevel, maxLevel = vswInstance.createTiles(inputFile = workingFile, outputFolder = tilesFolder, tileSize = outputTileSize, overlap = 0, suffix = suffix)
                    if verbose:
                        stop = time.clock()
                        print("Tiles generated. Time %s s" % str(stop-prev))
                        prev = stop
                else:
                    outputTileSize, nTilesWidth, nTilesHeight, minLevel, maxLevel = vswInstance.getMinMaxLevels(outputFilesFolder, suffix)
                nTilesXY = (nTilesWidth, nTilesHeight)
                if first:
                    firstNTilesXY = nTilesXY
                    firstTileSize = outputTileSize
                    first = False
                else:
                    assert outputTileSize == firstTileSize and nTilesXY == firstNTilesXY, "All images must have the same size."
                    
            except Exception as ex:
                print("Failed for %s: %s" % (fileName, str(ex)))

        # Done!
        if verbose:
            stop = time.clock()
            print("Total time %s s" % str(stop-start))

        return (nTilesXY, outputTileSize, maxLevel)
    
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

    def __largestNonEmptySquare(self, im, origTileSize):
        """
        Find largest square of non-empty image data
        """
        
        # Subsample to a size < n x n
        maxAllowed = 256
        maxSize1 = np.amax(im.shape)
        subSamplingFactor = np.ceil(maxSize1 / maxAllowed)
        imSS = im[ : : subSamplingFactor, : : subSamplingFactor]

        # Threshold 0-0 to get background
        imSSThr0 = imSS == 0

        # Label
        imSSLabel, imNLabels = nd.measurements.label(imSSThr0)

        # Compute histogram
        imSSLabelHist = nd.measurements.histogram(imSSLabel, 0, imNLabels, imNLabels + 1)

        # Loop through labelled image and remove pixels corresponding to small objects
        subSampledTileSize = origTileSize / subSamplingFactor
        subSampledTileSize2 = subSampledTileSize * subSampledTileSize
        sizeThr = subSampledTileSize2 / 2
        for r in range(0, imSS.shape[0] - 1):
            for c in range(0, imSS.shape[1] - 1):
                if imSSLabelHist[imSSLabel[r, c]] < sizeThr:
                    imSSThr0[r, c] = False

        # Compute a distance transform on the non-background
        imSSThr0 = imSSThr0 == False
        imSSDT = nd.morphology.distance_transform_cdt(imSSThr0, metric='chessboard') #Yes, it should be chessboard to get at square

        # Loop through the distance transform to find the maximum value and corresponding position
        maxDist = 0
        maxDistPos = (0, 0)
        for r in range(0, imSS.shape[0] - 1):
            for c in range(0, imSS.shape[1] - 1):
                if imSSDT[r, c] > maxDist:
                    maxDist = imSSDT[r, c]
                    maxDistPos = (r, c)

        # Compute parameters for sub image
        corner = (maxDistPos[0] - maxDist, maxDistPos[1] - maxDist)
        subImageSize = (2 * maxDist, 2 * maxDist)    
        cornerUpsampled = (corner[0] * subSamplingFactor, corner[1] * subSamplingFactor)
        subImageSizeUpsampled = (subImageSize[0] * subSamplingFactor, subImageSize[1] * subSamplingFactor)

        return (cornerUpsampled, subImageSizeUpsampled)

    def __computeRigidRegistration(self, im1, im2):
        """
        Compute rigid registration using Elastix.
        """
        
        reg = pyelastix.Elastix()
        # Get params and change a few values
        params = reg.get_default_params('rigid')
        params.MaximumNumberOfIterations = 200
        params.FinalGridSpacingInVoxels = 1
        params.NumberOfResolutions = 4
        params.ImagePyramidSchedule = [8, 8,  4, 4,  2, 2,  1, 1]

        # Apply the registration (im1 and im2 can be 2D or 3D)
        t, c = reg.computeRegistration(im1, im2, params, verbose = 0)
        rotation = -np.rad2deg(t[0]) #Elastix uses radians with inverted orientation
        translation = (-t[2], -t[1]) #Elastix has other x-y order
        rotationCenter = (t[1] + c[0], t[2] + c[1])
        return (rotation, rotationCenter, translation)
