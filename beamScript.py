import numpy as np
import matplotlib.pyplot as plt
import math
from scipy import ndimage

def apMeasPreProcessManip(apMeas,processType):
    if processType == "positive":
        apMeas = apMeas[apMeas >= 0]
    elif processType == "shift":
        minValue = apMeas.min()
        if minValue < 0:
            apMeas = math.fabs(minValue) + apMeas
    elif processType == "threshold":
        apMeas[apMeas < 0] = 0
    return apMeas

def beamScalingRowWise(rows,cols):
    beamArray = np.zeros((rows.size,cols.size), dtype=np.float64)
    maxValue = rows.max()

    for i in xrange(0,rows.size):
        beamArray[i,:] = (rows[i]/maxValue) * cols
    return beamArray

def beamScalingColWise(rows,cols):
    beamArray = np.zeros((rows.size,cols.size), dtype=np.float64)
    maxValue = cols.max()

    for i in xrange(0,cols.size):
        beamArray[:,i] = (cols[i]/maxValue) * rows
    return beamArray

def createAperturePSF(apertureDiameter,apertureStep):
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



beamApMeasXFilename = "20141216/20141216_Beam_profile_x_sma_ap.dat"
beamApMeasYFilename = "20141216/20141216_Beam_profile_y_sma_ap.dat"
apMeasPreProcessType = "None"
apDiameter = 10
apStep = 2

apertureX = np.loadtxt(beamApMeasXFilename,skiprows=7)
apertureXPosition = apertureX[:,0]
apertureXMeasurement = apertureX[:,1]

apertureY = np.loadtxt(beamApMeasYFilename,skiprows=7)
apertureYPosition = apertureY[:,0]
apertureYMeasurement = apertureY[:,1]

apertureXMeasurement = apMeasPreProcessManip(apertureXMeasurement,apMeasPreProcessType)
apertureYMeasurement = apMeasPreProcessManip(apertureYMeasurement,apMeasPreProcessType)

tempBeamArrayX = beamScalingRowWise(apertureYMeasurement,apertureXMeasurement)
tempBeamArrayY = beamScalingColWise(apertureYMeasurement,apertureXMeasurement)

convolvedBeamArray = (tempBeamArrayX + tempBeamArrayY) / 2

plt.imshow(convolvedBeamArray,cmap='jet')

aperturePSF = createAperturePSF(apDiameter,apStep)
