import numpy as np
import matplotlib.pyplot as plt
import math
from scipy import ndimage

def apMeasPreProcessManip(apMeas,processType):
    """Function that applies some processing to the aperture measurements
    prior to using them to create a beam.

    INPUTS:
        apMeas         - The i_pin measurement recordings from the aperture measurements
                        as a numpy array with float data types
        ProcessType    - The type of processing performed input as a string. There are 4
                        types of processing options available:
                        1) "positive" - uses only recordings where the i_pin reading were
                        positive. NOTE: This option is not recommended because it changes
                        the number of data points and it's likely that this option will
                        break the rest of the code. I haven't tested this but I haven't
                        had time to address the potential issue.
                        2) "shift" - If the minimum i_pin reading is negative than it will
                        add the absolute value of that reading to all values in the apMeas
                        array. This has the effect of shifting all the values up so every
                        value is above zero.
                        3)"threshold" - Any negative i_pin reading is set to zero.
                        4) "" - Any other string input will not do any preprocessing.

    OUTPUTS:
        apMeas         - A 1D numpy array of floats containing the processed i_pin readings
    """
    if processType == "positive":
        print 'Preprocessing type used on the i_pin measurements: "{}" '.format(processType)
        apMeas = apMeas[apMeas >= 0]     #Only use non-negative values (NOT RECOMMENDED)
    elif processType == "shift":
        print 'Preprocessing type used on the i_pin measurements: "{}" '.format(processType)
        minValue = apMeas.min()         #Find minimum value
        if minValue < 0:
            apMeas = math.fabs(minValue) + apMeas     #If minimum value is negative then shift all the values up so they are positive
    elif processType == "threshold":
        print 'Preprocessing type used on the i_pin measurements: "{}" '.format(processType)
        apMeas[apMeas < 0] = 0         #Set negative values to zero
    else:
        print 'No preprocessing has been performed on the i_pin measurements.'
    return apMeas

def beamScalingRowWise(rows,cols):
    """Scales the i_pin measurements in the horizontal direction according to the normalised
    values in the vertical directon to generate a 2D beam array.

    INPUTS:
        rows         -i_pin readings in the vertical direction as a numpy array of floats
        cols         -i_pin readings in the horizontal direction as a numpy array of floats

    OUTPUTS:
        beamArray    -2D numpy array of floats representing an expected set of i_pins values
                        if the readings were taken across the 2D area
    """
    print 'Performing row-wise scaling of i_pin measurements to generate temporary beam array'
    beamArray = np.zeros((rows.size,cols.size), dtype=np.float64)     #Preallocate beam array
    maxValue = rows.max()             #Get the maximum i_pin value in the vertical direction

    #For each row
    for i in xrange(0,rows.size):
        beamArray[i,:] = (rows[i]/maxValue) * cols  #Scale the horizontal i_pin measurements by the corresponding normalised vertical i_pin reading
    return beamArray

def beamScalingColWise(rows,cols):
    """Scales the i_pin measurements in the vertical direction according to the relative
    values in the horizontal directon to generate a 2D beam array.

    INPUTS:
        rows         -i_pin readings in the vertical direction as a numpy array of floats
        cols         -i_pin readings in the horizontal direction as a numpy array of floats

    OUTPUTS:
        beamArray    -2D numpy array of floats representing an expected set of i_pins values
                        if the readings were taken across the 2D area
    """
    print 'Performing column-wise scaling of i_pin measurements to generate temporary beam array'
    beamArray = np.zeros((rows.size,cols.size), dtype=np.float64)     #Preallocate beam array
    maxValue = cols.max()          #Get the maximum i_pin value in the horizontal direction

    #for each column
    for i in xrange(0,cols.size):
        beamArray[:,i] = (cols[i]/maxValue) * rows  #Scale the vertical i_pin measurements by the corresponding normalised horizontal i_pin reading
    return beamArray

def createAperturePSF(apertureDiameter,apertureStep):
    """Creates a matrix that acts as a Point Spread Function (PSF) of the aperture used to
    generate the i_pin readings.

    INPUTS:
        apertureDiameter        -The diameter of the aperture in microns
        apertureStep            -The incremental position at which consecutive i_pin readings
                                    are taken (in microns).

    OUTPUTS:
        psf                     - A 2D numpy array of integers, either 1's or 0's. 1's represent
                                    the area covered by the aperture.
    """
    print 'Generating the Point Spread Function of the {} micron diameter aperture with a measurement step of {} microns'.format(apertureDiameter,apertureStep)
    #Calculate the aperture radius from the diameter
    apertureRadius = apertureDiameter / 2.0

    ######################################################################
    ### Calculate from the radius how big the psf matrix should be and set all
    ### Values to zeros.
    ###
    #The matrix consists of all surrounding points from the centre of the
    #aperture that are within the area of the diameter - i.e. these matrix
    #elements could potentially contribute to the signal reading
    if apertureRadius%2 == 0:
        psf = np.zeros((apertureRadius + 1, apertureRadius + 1), dtype=np.int)
    else:
        psf = np.zeros((math.floor(apertureRadius), math.floor(apertureRadius)), dtype=np.int)

    ######################################################################
    #Get dimensions of the matrix (it should be a square matrix but we'll
    #take both dimensions anyway)
    dimensionsOfKernel = psf.shape
    #Get the index of the central point of the aperture psf
    indexOfPSFCentre = [dimensionsOfKernel[0] / 2, dimensionsOfKernel[1] / 2]
    #Set the central position in the psf to 1
    psf[indexOfPSFCentre[0], indexOfPSFCentre[1]] = 1

    ######################################################################
    ### Find the indices of the points that lie within the area of the aperture

    #Calculate the Euclidean distance of each point from the centre of the aperture
    distanceMatrix = ndimage.distance_transform_edt(psf==0,sampling=[apertureStep,apertureStep])
    #Find the points that lie within the area of the aperture
    booleanMatrix = distanceMatrix <= apertureRadius

    #Loop through matrix and set points that lie within the aperture to equal 1
    #and points that lie outside the aperture to equal 0
    for i in xrange(0,dimensionsOfKernel[0]):
        for j in xrange(0,dimensionsOfKernel[1]):
            if booleanMatrix[i,j] == True:
                psf[i,j] = 1

############ Input Variables ###########
beamApMeasXFilename = "20141216/20141216_Beam_profile_x_sma_ap.dat" #location of horizontal i_pin readings .dat file
beamApMeasYFilename = "20141216/20141216_Beam_profile_y_sma_ap.dat" #location of vertical i_pin readings .dat file
apMeasPreProcessType = "None"     #Processing type
apDiameter = 10     #Aperture diameter
apStep = 2         #Aperture step

#load the .dat file and store aperture data in a 2D numpy array.
#Then split the 2D array into two 1D arrays, one to store the i_pin readings.
#and one to store the aperture positions.
#Do this for both the vertical and horizontal directions.
apertureX = np.loadtxt(beamApMeasXFilename,skiprows=7)
apertureXPosition = apertureX[:,0]
apertureXMeasurement = apertureX[:,1]

apertureY = np.loadtxt(beamApMeasYFilename,skiprows=7)
apertureYPosition = apertureY[:,0]
apertureYMeasurement = apertureY[:,1]

#Apply preprocessing to the aperture measurements
apertureXMeasurement = apMeasPreProcessManip(apertureXMeasurement,apMeasPreProcessType)
apertureYMeasurement = apMeasPreProcessManip(apertureYMeasurement,apMeasPreProcessType)

#Create temporary 2D beam arrays using both row and column wise scaling
tempBeamArrayX = beamScalingRowWise(apertureYMeasurement,apertureXMeasurement)
tempBeamArrayY = beamScalingColWise(apertureYMeasurement,apertureXMeasurement)

#Take an average of the two temporary beam arrays to get a single convolved beam array
convolvedBeamArray = (tempBeamArrayX + tempBeamArrayY) / 2

plt.imshow(convolvedBeamArray,cmap='jet')

#Create the aperture Point Spread Function.
aperturePSF = createAperturePSF(apDiameter,apStep)
