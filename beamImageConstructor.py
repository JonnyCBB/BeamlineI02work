import numpy as np
import matplotlib.image as mpimg
import matplotlib.pyplot as plt

def generateBeamFromPNG(pngImage, redWeightValue, greenWeightValue, blueWeightValue):
    """Function that generates a beam array from a PNG file.

    INPUTS:
        pngImage           -The path to a png file that contains an image of an X-ray beam as a string.
        redWeightValue     -A scalar (float) value giving the weight of the red pixels in the png image
                            for the conversion to grayscale.
        greenWeightValue   -A scalar (float) value giving the weight of the green pixels in the png image
                            for the conversion to grayscale.
        blueWeightValue    -A scalar (float) value giving the weight of the blue pixels in the png image
                            for the conversion to grayscale.

    OUTPUTS:
        beamArray          -A 2D numpy array of integer values as a spatially resolved representation
                            of relative intensities. The values should be between 0 and 255 to be
                            compatible with the .pgm file format.
    """

    rgbBeamImage = mpimg.imread(pngImage)     #read the png file
    grayscaleBeamImage = rgb2Grayscale(rgbBeamImage, redWeightValue, greenWeightValue, blueWeightValue)     #convert from rgb values to grayscale
    beamArray = scaleArray(grayscaleBeamImage)     #scale the array whilst converting floats to integer values
    return beamArray     #return beam array

def rgb2Grayscale(imageFile,redWeight,greenWeight,blueWeight):
    """Function to convert rgb images to grayscale

    INPUTS:
        imageFile              -An N x M x 3 (3D) array of floats representing the RGB values of an image.
        redWeight              -A scalar (float) value giving the weight of the red pixels in the png image
                                for the conversion to grayscale.
        greenWeight            -A scalar (float) value giving the weight of the green pixels in the png image
                                for the conversion to grayscale.
        blueWeight             -A scalar (float) value giving the weight of the blue pixels in the png image
                                for the conversion to grayscale.

    OUTPUTS:
        grayscaleimage         -An N x M (2D) array of floats representing the weighted average grayscale values
                                of the image.
    """

    #Extract the red, green and blues pixels from the image
    redPixels, greenPixels, bluePixels = imageFile[:,:,0], imageFile[:,:,1], imageFile[:,:,2]
    #Calculate the weighted average
    grayscaleImage = (redWeight * redPixels) + (greenWeight * greenPixels) + (blueWeight * bluePixels)

    return grayscaleImage     #return the grayscale image

def scaleArray(array):
    """ Function to scale the values of the array to range between 0 and 255.
    This function also converts all values to integer values so it's compatible
    with the pgm image file format.

    INPUTS:
        array           -A numpy array of floats/ints

    OUTPUTS:
        scaledArray     -A scaled version of the input array. Values in the array are
                            converted to integers.
    """
    scalingValue = 255/array.max()
    scaledArrayAsFloats = np.around(array * scalingValue)
    scaledArray = scaledArrayAsFloats.astype(int)
    return scaledArray

def parseBeamImageFilename(imageFilename):
    """Parse the image file name to get the information about the camera settings used to capture the image

    INPUTS:
        imageFilename         -String pointing to the png image file.

    OUTPUTS:
        (Tuple of values)     -The output is a tuple of values giving various camera settings set during
                                image capture. These are:
                                zoom
                                gain
                                exposureTime
                                transmission
                                slits
    """

    splitFilename = imageFilename.split("_")     #Split the filename by underscores
    slits = np.zeros(2, dtype=np.float)          #Preallocate array to store the slit dimensions
    #Loop through each split entry of the filename
    for entry in splitFilename:
        if entry != "beam":     #ignore the "beam" entry
            if ".png" in entry:             #check if the entry contains the file extension. If so then remove it.
                splitEntry = entry.split(".")
                entry = splitEntry[0]
            if entry.isdigit():             #check if the entry only contains numbers. If so then it's the vertical slit distance
                slits[1] = float(entry)
            else:                           #else check the different camera settings and assign them to the correct variable.
                firstLetter = entry[0]
                if firstLetter == 'z':
                    zoom = float(entry[1:])
                elif firstLetter == 'g':
                    gain = float(entry[1:])
                elif firstLetter == 'e':
                    exposureTime = float(entry[1:])
                elif firstLetter == 't':
                    transmission = float(entry[1:])
                elif firstLetter == 's':
                    slits[0] = float(entry[1:])

    return(zoom,gain,exposureTime,transmission,slits)     #Return tuple

############# Main script
beamFromPNG = generateBeamFromPNG('beam_z20_g1_e70_s100_60_t100.png', 0.3, 0.6, 0.1)
camera = parseBeamImageFilename('beam_z20_g1_e70_s100_60_t100.png')
