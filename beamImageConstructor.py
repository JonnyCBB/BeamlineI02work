import numpy as np
import matplotlib.image as mpimg
import matplotlib.pyplot as plt

def generateBeamFromPNG(pngImage, redWeightValue, greenWeightValue, blueWeightValue):
    rgbBeamImage = mpimg.imread(pngImage)
    grayscaleBeamImage = rgb2Grayscale(rgbBeamImage, redWeightValue, greenWeightValue, blueWeightValue)
    beamImage = scaleArray(grayscaleBeamImage)
    return beamImage

def rgb2Grayscale(imageFile,redWeight,greenWeight,blueWeight):

    redPixels, greenPixels, bluePixels = imageFile[:,:,0], imageFile[:,:,1], imageFile[:,:,2]
    grayscaleImage = (redWeight * redPixels) + (greenWeight * greenPixels) + (blueWeight * bluePixels)

    return grayscaleImage

def scaleArray(array):
    scalingValue = 255/array.max()
    scaledArrayAsFloats = np.around(array * scalingValue)
    scaledArray = scaledArrayAsFloats.astype(int)
    return scaledArray

beam2 = generateBeamFromPNG('beam_z20_g1_e70_s100_60_t100.png', 0.3, 0.6, 0.1)
