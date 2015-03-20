import BeamModule
import RaddoseBeamRun
import numpy as np
import matplotlib.pyplot as plt
import math
import seaborn as sns
from mpl_toolkits.mplot3d import Axes3D
from scipy.ndimage import filters

def createWindowAroundCentroid(array, xDistance, yDistance):
    xCen, yCen = BeamModule.findCentroid(beamFromPNGImage.beamArray)
    pixelsInXDir = xDistance/(2.0 * beamFromPNGImage.beamPixelSize[0])
    pixelsInYDir = yDistance/(2.0 * beamFromPNGImage.beamPixelSize[1])
    xMin = xCen - pixelsInXDir
    xMax = xCen + pixelsInXDir
    yMin = yCen - pixelsInYDir
    yMax = yCen + pixelsInYDir
    return xMin, xMax, yMin, yMax

#beamFromAperture = BeamModule.Beam.initialiseBeamFromApMeas("20141216/20141216_Beam_profile_x_sma_ap.dat","20141216/20141216_Beam_profile_y_sma_ap.dat","threshold", 10, 2, "beamtest.pgm",
#                                                            "nogaussconv",True,"smoothweiner")
beamFromPNGImage = BeamModule.Beam.initialiseBeamFromPNG("beam_z20_g1_e70_s100_60_t100.png",0.3,0.6,0.1,"anotherbeamtest.pgm")
#beamFromPGM = BeamModule.Beam.initialiseBeamFromPGM("jonny.pgm")
#raddose = RaddoseBeamRun.RunRaddose(beamFromAperture,"testInput.txt","Max")

# Function to create a specified window around the centroid of the image
# xDistance = beamFromPNGImage.collimation[0]
# yDistance = beamFromPNGImage.collimation[1]
# window = createWindowAroundCentroid(beamFromPNGImage, xDistance, yDistance)

# xMin = np.floor(window[0])
# xMax = np.ceil(window[1])
# yMin = np.floor(window[2])
# yMax = np.ceil(window[3])

# beamArray = beamFromPNGImage.beamArray
# backgroundList = []
# for i in xrange(0,beamArray.shape[0]):
#     for j in xrange(0,beamArray.shape[1]):
#         if j < xMin or j > xMax or i < yMin or i > yMax:
#             backgroundList.append(beamArray[i,j])
# averageBackground = np.mean(backgroundList)
# beamMinusBackground = beamArray - averageBackground
# beamMinusBackground[beamMinusBackground < 0] = 0
# beamMinusBackground = np.around(beamMinusBackground)
# beamMinusBackground = beamMinusBackground.astype(int)
# xmax,ymax = beamFromAperture.beamArray.shape
# x = np.linspace(0,ymax,ymax)
# y = np.linspace(0,xmax,xmax)
# X,Y = np.meshgrid(x,y)

# fig = plt.figure(figsize=(25,10))
# ax = fig.add_subplot(1, 2, 1, projection='3d')
# p = ax.plot_surface(X,Y,beamFromAperture.beamArray, rstride=1, cstride=1, cmap='seismic', linewidth=0, antialiased=False)
# cb = fig.colorbar(p, shrink=0.5)

# ax = fig.add_subplot(1, 2, 2, projection='3d')
# p = ax.plot_surface(X, Y, a, rstride=4, cstride=4, cmap='seismic', linewidth=0, antialiased=False)
# cb = fig.colorbar(p, shrink=0.5)
