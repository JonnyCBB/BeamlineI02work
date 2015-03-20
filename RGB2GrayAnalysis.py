#Script to analyse the effect of colour weighted when converting RGB images
#to grayscale images.

############################################################
# Sort imports
############################################################
import BeamModule
import RaddoseBeamRun
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

############################################################
# Choose a selection of color weights between 0 and 255
############################################################
#Create empty lists
redWeights = []
greenWeights = []
blueWeights = []

#Append the first three weighting regimes
redWeights.append(1)
greenWeights.append(0)
blueWeights.append(0)
redWeights.append(0)
greenWeights.append(1)
blueWeights.append(0)
redWeights.append(0)
greenWeights.append(0)
blueWeights.append(1)

#Add random weights to the color weighting lists
for i in xrange(1,1):
    redWeights.append(float(np.random.uniform(0, 1 ,1)))
    greenWeights.append(float(np.random.uniform(0, 1 - redWeights[-1], 1)))
    blueWeights.append(1 - (redWeights[-1] + greenWeights[-1]))


############################################################
# Create beams
############################################################
beamList = []
for red in redWeights:
    for green in greenWeights:
        for blue in blueWeights:
            beamList.append(BeamModule.Beam.initialiseBeamFromPNG("beam_z20_g1_e70_s100_60_t100.png",red,green,blue,"rgb2grayanalysis.pgm"))
