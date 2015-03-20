#Script to analyse the effect of colour weighted when converting RGB images
#to grayscale images.

############################################################
# Sort imports
############################################################
import BeamModule
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
numberOfAdditionalBeams = 50
for i in xrange(0,numberOfAdditionalBeams):
    redWeights.append(float(np.random.uniform(0, 1 ,1)))
    greenWeights.append(float(np.random.uniform(0, 1 - redWeights[-1], 1)))
    blueWeights.append(1 - (redWeights[-1] + greenWeights[-1]))

############################################################
# Create beams
############################################################
beamList = []
for i in xrange(0,len(redWeights)):
    beamList.append(BeamModule.Beam.initialiseBeamFromPNG("beam_z20_g1_e70_s100_60_t100.png",
                                                          redWeights[i], greenWeights[i], blueWeights[i],
                                                          "rgb2grayanalysis.pgm"))

##############################################################
# Get dose Type
##############################################################
b = beamList[0]
if "DWD" in b.doseType:
    doseType = "Average Diffraction Weighted Dose"
elif "Max" in b.doseType:
    doseType = "Max Dose"
elif "Av" in b.doseType:
    doseType = "Average Dose (Whole Crystal)"

############################################################
# Plot dose for each beam
############################################################
sns.set_context("talk")
fig = plt.figure(1)
ax = fig.add_subplot(1, 1, 1)
labels = []
doseList = []
for i in xrange(0,len(beamList)):
    beam = beamList[i]
    labels.append('beam' + str(i))
    doseList.append(beam.dose)
    plt.plot(i, beam.dose, marker='o', markersize=10, linestyle="None")
    #If it's the first iteration then set the min and max dose
    if i == 0:
        minDose = beam.dose
        maxDose = beam.dose

    #Check if the dose value
    if beam.dose < minDose:
        minDose = beam.dose
    elif beam.dose > maxDose:
        maxDose = beam.dose
avgDose = np.mean(doseList)
plt.xlim((-1,len(beamList)))
plt.xticks(np.arange(len(beamList)), labels, rotation=45)
plt.minorticks_off()
plt.title("Effect of RGB weights on the calculated dose profile")
plt.ylabel(doseType + " (MGy)" )
plt.show()

error = (maxDose - minDose)/avgDose * 100
print "Maximum Dose is: " + str(round(maxDose, 2)) + " MGy"
print "Minimum Dose is: " + str(round(minDose, 2)) + " MGy"
print "Average Dose is: " + str(round(avgDose, 2)) + " MGy"
print "Error is: " + str(round(error, 2)) + "%"
