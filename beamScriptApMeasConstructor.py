import numpy as np
import matplotlib.pyplot as plt
import math
from scipy import ndimage
from scipy import signal
from scipy.optimize import minimize
from skimage import restoration

def generateBeamFromApMeas(beamApMeasXFilename, beamApMeasYFilename, beamPostProcessingType, apDiameter, apStep):
    """Create a beam array from aperture scan measurements

    INPUTS:
        beamApMeasXFilename     -string with the location of the aperture scan measurements in the horizontal (x)
                                    direction.
        beamApMeasYFilename     -string with the location of the aperture scan measurements in the vertical (y)
                                    direction.
        beamPostProcessingType  -string giving the type of processing that should be carried out on the beam
                                    array once the deconvolution has taken place.
        apertureDiameter        -The diameter of the aperture in microns
        apertureStep            -The incremental position at which consecutive i_pin readings
                                    are taken (in microns).

    OUTPUTS:
        beamArray               -A 2D numpy array of integer values as a spatially resolved representation
                                    of relative intensities. The values should be between 0 and 255 to be
                                    compatible with the .pgm file format.
    """
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

    #Create temporary 2D beam arrays using both row and column wise scaling
    tempBeamArrayX = beamScalingRowWise(apertureYMeasurement,apertureXMeasurement)
    tempBeamArrayY = beamScalingColWise(apertureYMeasurement,apertureXMeasurement)

    #Take an average of the two temporary beam arrays to get a single convolved beam array
    convolvedBeamArray = (tempBeamArrayX + tempBeamArrayY) / 2

    #Create the aperture Point Spread Function.
    aperturePSF = createAperturePSF(apDiameter,apStep)

    #Deconvolve the beam image
    blurredBeamTuple = restoration.unsupervised_wiener(convolvedBeamArray, aperturePSF)
    blurredBeamArray = blurredBeamTuple[0]     #Get beam array

    initialBeamNoiseGuess = 1
    res = minimize(lambda beamNoise: deblurBeamObjectiveFunction(beamNoise,blurredBeamArray,aperturePSF,
                                                                apertureXMeasurement,apertureYMeasurement),
                   initialBeamNoiseGuess, method='nelder-mead', options={'xtol': 1e-8, 'disp': True})

    #Deblur the image
    deconvolvedBeamArray = signal.wiener(blurredBeamArray,aperturePSF.shape,res.x[0])

    #Apply post processing on the beam array
    processedBeamArray = beamPostProcessManip(deconvolvedBeamArray,beamPostProcessingType)

    #Transform all values so they lie between 0 and 255
    scalingValue = 255 / processedBeamArray.max()
    arrayAsFloats = np.around(processedBeamArray * scalingValue)
    beamArray = arrayAsFloats.astype(int)

    return beamArray

def beamPostProcessManip(beamArray,processType):
    """Function that applies some processing to the aperture measurements
    prior to using them to create a beam.

    INPUTS:
        beamArray      - The deconvolved beam array as a 2D numpy array of floats.
        ProcessType    - The type of processing performed input as a string. There are 3
                        types of processing options available:
                        1) "shift" - If the minimum beam array reading is negative than it will
                        add the absolute value of that reading to all values in the beamArray
                        array. This has the effect of shifting all the values up so every
                        value is above zero.
                        2)"threshold" - Any negative beam array reading is set to zero.
                        3) "" - Any other string input will not do any preprocessing.

    OUTPUTS:
        beamArray         - A 2D numpy array of floats containing the processed beam array values.
    """
    if processType == "positive":
        print 'Post processing type used on the beam array: "{}" '.format(processType)
        beamArray = beamArray[beamArray >= 0]     #Only use non-negative values (NOT RECOMMENDED)
    elif processType == "shift":
        print 'Post processing type used on the beam array: "{}" '.format(processType)
        minValue = beamArray.min()         #Find minimum value
        if minValue < 0:
            beamArray = math.fabs(minValue) + beamArray     #If minimum value is negative then shift all the values up so they are positive
    elif processType == "threshold":
        print 'Post processing type used on the beam array: "{}" '.format(processType)
        beamArray[beamArray < 0] = 0         #Set negative values to zero
    else:
        print 'No post processing has been performed on the beam array measurements.'
    return beamArray

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

    return psf     #Return the point spread function

##########################################################
def simulateApertureScans(beamArray,psf):
    """Function that uses a beam array and a point spread function representing the aperture
    to simulate the aperture scans carried out on the I02 beamline at Diamond Light Source.

    INPUTS:
        beamArray         -A 2D numpy array of floats that represents the theoretical deconvolved i_pin
                            readings at each point in the space.
        psf               -A 2D numpy array of floats that represents the aperture used to measure the
                            beam intensity.

    OUTPUTS:
        simApScanX        -1D numpy array of floats containing i_pin readings from the simulated
                            aperture scan in the horizontal direction.
        simApScanY        -1D numpy array of floats containing i_pin readings from the simulated
                            aperture scan in the vertical direction.
    """

    #Calculate total number of elements to be added in each dimension of beam matrix
    matrixBuffer = int(2*math.floor(psf.shape[0]/2))

    #Preallocate the buffered matrix. It's buffered with zeros around the outside
    bufferedBeamMatrix = np.zeros((beamArray.shape[0] + matrixBuffer, beamArray.shape[1] + matrixBuffer), dtype=np.float)

    #Check if the beam matrix has been cropped in the second dimension by the deconvolution process
    if bufferedBeamMatrix.shape[0] == bufferedBeamMatrix.shape[1]:
        beamArrayCrop = 0
    else:
        beamArrayCrop = bufferedBeamMatrix.shape[0] - bufferedBeamMatrix.shape[1]

    #Fill bufferedBeamMatrix with elements from the actual beam array.
    for i in xrange(matrixBuffer / 2, bufferedBeamMatrix.shape[0] - (matrixBuffer / 2)):
        for j in xrange(matrixBuffer / 2, bufferedBeamMatrix.shape[1] - (matrixBuffer / 2)):
            bufferedBeamMatrix[i,j] = beamArray[i - (matrixBuffer / 2), j - (matrixBuffer / 2)]

    #Preallocate arrays to contain simulated aperture scan measurements
    simApScanX = np.zeros(beamArray.shape[1], dtype=np.float)
    simApScanY = np.zeros(beamArray.shape[0], dtype=np.float)

    #Get the central element for horizontal (X) and the vertical (Y) scans
    apScanRow = int(math.floor(beamArray.shape[0] / 2) + matrixBuffer / 2)
    apScanCol = int(math.floor((beamArray.shape[1] + beamArrayCrop) / 2) + matrixBuffer / 2)

    #Simulate horizontal aperture scan
    for col in xrange(matrixBuffer / 2, bufferedBeamMatrix.shape[1] - (matrixBuffer / 2)):
        for i in xrange(0,psf.shape[0]):
            for j in xrange(0,psf.shape[1]):
                a = i - matrixBuffer / 2
                b = j - matrixBuffer / 2
                simApScanX[col - matrixBuffer] += bufferedBeamMatrix[apScanRow + a, col + b] * psf[i,j]

    #Simulate vertical aperture scan
    for row in xrange(matrixBuffer / 2, bufferedBeamMatrix.shape[0] - (matrixBuffer / 2)):
        for i in xrange(0,psf.shape[0]):
            for j in xrange(0,psf.shape[1]):
                a = i - matrixBuffer / 2
                b = j - matrixBuffer / 2
                simApScanY[row - matrixBuffer] += bufferedBeamMatrix[row + a, apScanCol + b] * psf[i,j]

    return simApScanX, simApScanY

def rootMeanSquaredDeviation(xPredicted,xMeasured):
    """Find the root mean squared deviation between two 1D numpy arrays.
        Note: both arrays must be the same size.

    INPUTS:
        xPredicted        -1D numpy array of floats. These are the predicted values
        xMeasured         -1D numpy array of floats. These are the measured values

    OUTPUTS:
        rmsd              -A scalar float value representing the root mean squared deviation
    """
    rmsd = np.sqrt(np.sum(np.square(xPredicted - xMeasured)))

    return rmsd


def deblurBeamObjectiveFunction(noiseRatio,beamArray,psf,actualApMeasurementX,actualApMeasurementY):
    """Objective function used to deblur the deconvolved beam image.
    An objective funtion is a function whos output is required to be optimal (in this case our optimal value
    is the minimal one) by some optimisation routine. One, or many arguments can be altered by the
    optimisation routine. In this case that parameter is the noiseRatio.

    INPUTS:
        noiseRatio            -A scalar value that represents the noise-power term in the wiener deconvolution
                                function
        beamArray             -A 2D numpy array of floats that represents the theoretical deconvolved i_pin
                                readings at each point in the space.
        psf                   -A 2D numpy array of floats that represents the aperture used to measure the
                                beam intensity.
        actualApMeasurementX  -Measured i_pin readings in the horizontal direction as a 1D numpy array of floats
        actualApMeasurementY  -Measured i_pin readings in the vertical direction as a 1D numpy array of floats

    OUTPUTS:
        totalRMSD             -The sum of the root mean squared deviations of the theoretical and measured
                                i_pin readings.
    """
    #Deblur the image
    deconvolvedBeamArray = signal.wiener(beamArray,psf.shape,noiseRatio)

    #simulate the aperture scans
    simulatedApMeasurementsX,simulatedApMeasurementsY = simulateApertureScans(deconvolvedBeamArray,psf)

    #Calculate the root mean squared deviations (rmsd)
    rmsdX = rootMeanSquaredDeviation(simulatedApMeasurementsX, actualApMeasurementX[0:-1])
    rmsdY = rootMeanSquaredDeviation(simulatedApMeasurementsY, actualApMeasurementY)

    #Add the rmsds
    totalRMSD = rmsdX + rmsdY

    return totalRMSD

########## Main program
beam = generateBeamFromApMeas("20141216/20141216_Beam_profile_x_sma_ap.dat",
                             "20141216/20141216_Beam_profile_y_sma_ap.dat",
                              "threshold",
                              10,
                              2)
plt.imshow(beam,cmap='gray')
plt.show()
