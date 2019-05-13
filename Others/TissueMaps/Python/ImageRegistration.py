# File: ImageRegistration.py
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
import ctypes

class ImageRegistration(object):
    """
    Read, register and tile large images from e.g. slide scanner
    """
    
    def __init__(self):
        pass
    
    def run(self, inputFolder, \
            referenceFileName = None, \
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
        if verbose:
            print("Reference file: %s" % referenceFileName)

        if referenceFileName is not None:
            referenceFileName = referenceFileName.strip()
            referenceFileName = referenceFileName.replace('\\', '/')

        start = time.clock()
        prev = start

        if referenceFileName is not None:
            if verbose:
                print("Start opening images") 
            # Read image 1
            referenceFileName = referenceFileName.strip()
            path, ext = os.path.splitext(referenceFileName)
            if ext == ".tif":
                im1 = tifffile.imread(referenceFileName)
                pilIm1 = Image.fromarray(im1)
            else:
                pilIm1 = Image.open(referenceFileName)
                im1 = np.array(pilIm1)
            if verbose:
                stop = time.clock()
                print("Image1 opened. Time %s s" % str(stop-prev))
                print(referenceFileName)
                print("im1.shape %s pilIm1.size %s" % (str(im1.shape), str(pilIm1.size)))
                prev = stop
                   
            for fileName in os.listdir(inputFolder):
                try:
                    if not inputFolder.endswith('/'):
                        inputFolder = inputFolder + '/'
                    floatingFileName = inputFolder + fileName
                    floatingFileName = floatingFileName.strip()
                    print(floatingFileName)
                    if os.path.isfile(floatingFileName):
                        print("Start open image")
                        # Read image 2 
                        path, ext = os.path.splitext(floatingFileName)
                        if ext == ".tif":
                            im2 = tifffile.imread(floatingFileName)
                            pilIm2 = Image.fromarray(im2)
                        else:
                            pilIm2 = Image.open(floatingFileName)
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
                            print("Start computing transformation")
                        # Compute transformation
                        rotation, rotationCenter, translation = self.__computeRigidRegistration(im1, im2)
                                            
                        if verbose:
                            stop = time.clock()
                            print("Transformation computed. Time %s s" % str(stop-prev))
                            print("Transformation parameters: rotation %s rotationCenter %s translation %s" % (str(rotation), str(rotationCenter), str(translation)))
                            prev = stop
                            
                        pilIm2 = self.__scaleRotateTranslate(pilIm2, rotation, rotationCenter = rotationCenter, translation = translation )
                        if verbose:
                            stop = time.clock()
                            print("Transformation applied. Time %s s" % str(stop-prev))
                            prev = stop            
                        if verbose:
                            print("Start save transformed image")
                        # Save the transformed image
                        outputFolder = inputFolder + 'Transformed/'            
                        if not os.path.exists(outputFolder):
                            os.makedirs(outputFolder)
                        workingFile = outputFolder + fileName
                        pilIm2.save(workingFile)
                        if not os.path.isfile(workingFile):
                            print("No file saved: %s" % workingFile)
                        if verbose:
                            stop = time.clock()
                            print("Transformed Image saved. Time %s s" % str(stop-prev))
                            prev = stop
                except Exception as e:
                    print("Exception: %s" % e)
                    pass

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
        t, c = reg.computeRegistration(im1, im2, params, 0, 0)
        rotation = -np.rad2deg(t[0]) #Elastix uses radians with inverted orientation
        translation = (-t[2], -t[1]) #Elastix has other x-y order
        rotationCenter = (t[1] + c[0], t[2] + c[1])
        return (rotation, rotationCenter, translation)
