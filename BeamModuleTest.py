import BeamModule
import RaddoseBeamRun
import numpy as np
import matplotlib.pyplot as plt
import math
import seaborn as sns
from mpl_toolkits.mplot3d import Axes3D
from scipy.ndimage import filters

def findCentroid(array):
        h, w = array.shape
        ygrid, xgrid  = np.mgrid[0:h:1, 0:w:1]
        xcen, ycen = xgrid[array == 255].mean(), ygrid[array == 255].mean()
        return xcen, ycen

beamFromAperture = BeamModule.Beam.initialiseBeamFromApMeas("20141216/20141216_Beam_profile_x_sma_ap.dat","20141216/20141216_Beam_profile_y_sma_ap.dat","threshold", 10, 2, "beamtest.pgm")
beamFromPNGImage = BeamModule.Beam.initialiseBeamFromPNG("beam_z20_g1_e70_s100_60_t100.png",0.3,0.6,0.1,"anotherbeamtest.pgm")
#beamFromPGM = BeamModule.Beam.initialiseBeamFromPGM("jonny.pgm")
raddose = RaddoseBeamRun.RunRaddose(beamFromAperture,"testInput.txt")

# xCen,yCen = findCentroid(beamFromPNGImage.beamArray)

# plt.figure(1)
# plt.imshow(beamFromPNGImage.beamArray,cmap='gray')
# plt.plot(xCen,yCen,'ro')
# plt.show()

# arrayHeight, arrayWidth = beamFromPNGImage.beamArray.shape

# if arrayHeight/2.0 < yCen:
#     leftIndex = math.floor(yCen - arrayHeight/2.0)
# else:
#     leftIndex = 0
#     rightIndex = math.ceil(2 * yCen)

# if arrayWidth/2.0 < xCen:
#     topIndex = math.floor(xCen - arrayWidth/2.0)
# else:
#     topIndex = 0
#     bottomIndex = math.ceil(2 * xCen)

# if leftIndex != 0 and topIndex != 0:
#     b = beamFromPNGImage.beamArray[leftIndex:,topIndex:]
# elif leftIndex != 0:
#     b = beamFromPNGImage.beamArray[leftIndex:,topIndex:bottomIndex]
# elif topIndex != 0:
#     b = beamFromPNGImage.beamArray[leftIndex:rightIndex,topIndex:]
# else:
#     b = beamFromPNGImage.beamArray[leftIndex:rightIndex,topIndex:bottomIndex]

# plt.figure(2)
# plt.imshow(b,cmap='gray')

a = filters.gaussian_filter(beamFromAperture.beamArray,1)

xmax,ymax = beamFromAperture.beamArray.shape
x = np.linspace(0,ymax,ymax)
y = np.linspace(0,xmax,xmax)
X,Y = np.meshgrid(x,y)

fig = plt.figure(figsize=(25,10))
ax = fig.add_subplot(1, 2, 1, projection='3d')
p = ax.plot_surface(X,Y,beamFromAperture.beamArray, rstride=1, cstride=1, cmap='seismic', linewidth=0, antialiased=False)
cb = fig.colorbar(p, shrink=0.5)

ax = fig.add_subplot(1, 2, 2, projection='3d')
p = ax.plot_surface(X, Y, a, rstride=4, cstride=4, cmap='seismic', linewidth=0, antialiased=False)
cb = fig.colorbar(p, shrink=0.5)
