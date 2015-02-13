import numpy as np
import matplotlib.pyplot as plt
import math

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
    beamArray = np.zeros((rows.size,cols.size))
    maxValue = rows.max()

    for i in xrange(0,rows.size):
        beamArray[i,:] = (rows[i]/maxValue) * cols
    return beamArray

def beamScalingColWise(rows,cols):
    beamArray = np.zeros((rows.size,cols.size))
    maxValue = cols.max()

    for i in xrange(0,cols.size):
        beamArray[:,i] = (cols[i]/maxValue) * rows
    return beamArray


beamApMeasXFilename = "20141216/20141216_Beam_profile_x_sma_ap.dat"
beamApMeasYFilename = "20141216/20141216_Beam_profile_y_sma_ap.dat"
apMeasPreProcessType = "None"

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

plt.imshow(convolvedBeamArray)
plt.show()
