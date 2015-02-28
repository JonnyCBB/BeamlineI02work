import BeamModule
import numpy as np
import matplotlib.pyplot as plt

def findCentroid(array):
        h, w = array.shape
        ygrid, xgrid  = np.mgrid[0:h:1, 0:w:1]
        xcen, ycen = xgrid[array == 255].mean(), ygrid[array == 255].mean()
        return xcen, ycen

#beamFromAperture = BeamModule.Beam.initialiseBeamFromApMeas("20141216/20141216_Beam_profile_x_sma_ap.dat","20141216/20141216_Beam_profile_y_sma_ap.dat","threshold", 10, 2, "beamtest.pgm")
beamFromPNGImage = BeamModule.Beam.initialiseBeamFromPNG("beam_z20_g1_e70_s100_60_t100.png",0.3,0.6,0.1,"anotherbeamtest.pgm")
#beamFromPGM = BeamModule.Beam.initialiseBeamFromPGM("jonny.pgm")

a = findCentroid(beamFromPNGImage.beamArray)

plt.figure(1)
plt.imshow(beamFromPNGImage.beamArray,cmap='gray')
plt.plot(a[0],a[1],'ro')
plt.show()
