import BeamModule
import RaddoseBeamRun
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

"""
This script analyses beams generated from the aperture measurements taken at
Diamond Light Source beamline I0apertureStep.
Generating beams from aperture measurements with different processing parameters
Will give different dose values for any given crystal and strategy. This script
analyses the sizes of the errors caused by processing the measurements differently.
"""

#Give location of the aperture measurement data:
xApertureMeas = "20141216/20141216_Beam_profile_x_sma_ap.dat"
yApertureMeas = "20141216/20141216_Beam_profile_y_sma_ap.dat"
apertureDiameter = 10
apertureStep = 2

################################################################
#First generate beams with different processing parameters.
################################################################
noDeconvGaussThresh = BeamModule.Beam.initialiseBeamFromApMeas(xApertureMeas,yApertureMeas,
                                                            "threshold", apertureDiameter, apertureStep,
                                                               "beam.pgm","gaussfitconv",False,"smoothweiner")

noDeconvAvDataThresh = BeamModule.Beam.initialiseBeamFromApMeas(xApertureMeas,yApertureMeas,
                                                            "threshold", apertureDiameter, apertureStep,
                                                                "beam.pgm","AvData",False,"smoothweiner")

deconvAvDataWeinerSmoothThresh = BeamModule.Beam.initialiseBeamFromApMeas(xApertureMeas,yApertureMeas,
                                                            "threshold", apertureDiameter, apertureStep,
                                                                          "beam.pgm","AvData",True,"smoothweiner")

deconvAvDataGaussSmoothThresh = BeamModule.Beam.initialiseBeamFromApMeas(xApertureMeas,yApertureMeas,
                                                            "threshold", apertureDiameter, apertureStep,
                                                                         "beam.pgm","AvData",True,"smoothgauss")

deconvAvDataNoSmoothThresh = BeamModule.Beam.initialiseBeamFromApMeas(xApertureMeas,yApertureMeas,
                                                            "threshold", apertureDiameter, apertureStep,
                                                                      "beam.pgm","AvData",True,"noSmooth")

deconvGaussWeinerSmoothThresh = BeamModule.Beam.initialiseBeamFromApMeas(xApertureMeas,yApertureMeas,
                                                            "threshold", apertureDiameter, apertureStep,
                                                                         "beam.pgm","gaussfitconv",True,"smoothweiner")

deconvGaussGaussSmoothThresh = BeamModule.Beam.initialiseBeamFromApMeas(xApertureMeas,yApertureMeas,
                                                            "threshold", apertureDiameter, apertureStep,
                                                                        "beam.pgm","gaussfitconv",True,"smoothgauss")

deconvGaussNoSmoothThresh = BeamModule.Beam.initialiseBeamFromApMeas(xApertureMeas,yApertureMeas,
                                                            "threshold", apertureDiameter, apertureStep,
                                                                     "beam.pgm","gaussfitconv",True,"noSmooth")

noDeconvGaussShift = BeamModule.Beam.initialiseBeamFromApMeas(xApertureMeas,yApertureMeas,
                                                            "shift", apertureDiameter, apertureStep,
                                                              "beam.pgm","gaussfitconv",False,"smoothweiner")

noDeconvAvDataShift = BeamModule.Beam.initialiseBeamFromApMeas(xApertureMeas,yApertureMeas,
                                                            "shift", apertureDiameter, apertureStep,
                                                               "beam.pgm","AvData",False,"smoothweiner")

deconvAvDataWeinerSmoothShift = BeamModule.Beam.initialiseBeamFromApMeas(xApertureMeas,yApertureMeas,
                                                            "shift", apertureDiameter, apertureStep,
                                                                        "beam.pgm","AvData",True,"smoothweiner")

deconvAvDataGaussSmoothShift = BeamModule.Beam.initialiseBeamFromApMeas(xApertureMeas,yApertureMeas,
                                                            "shift", apertureDiameter, apertureStep,
                                                                        "beam.pgm","AvData",True,"smoothgauss")

deconvAvDataNoSmoothShift = BeamModule.Beam.initialiseBeamFromApMeas(xApertureMeas,yApertureMeas,
                                                            "shift", apertureDiameter, apertureStep,
                                                                     "beam.pgm","AvData",True,"noSmooth")

deconvGaussWeinerSmoothShift = BeamModule.Beam.initialiseBeamFromApMeas(xApertureMeas,yApertureMeas,
                                                            "shift", apertureDiameter, apertureStep,
                                                                        "beam.pgm","gaussfitconv",True,"smoothweiner")

deconvGaussGaussSmoothShift = BeamModule.Beam.initialiseBeamFromApMeas(xApertureMeas,yApertureMeas,
                                                            "shift", apertureDiameter, apertureStep,
                                                                       "beam.pgm","gaussfitconv",True,"smoothgauss")

deconvGaussNoSmoothShift = BeamModule.Beam.initialiseBeamFromApMeas(xApertureMeas,yApertureMeas,
                                                            "shift", apertureDiameter, apertureStep,
                                                                    "beam.pgm","gaussfitconv",True,"noSmooth")

##############################################################
# Get dose Type
##############################################################
if "DWD" in noDeconvGaussShift.doseType:
    doseType = "Average Diffraction Weighted Dose"
elif "Max" in noDeconvGaussShift.doseType:
    doseType = "Max Dose"
elif "Av" in noDeconvGaussShift.doseType:
    doseType = "Average Dose (Whole Crystal)"

##############################################################
# Run RADDOSE-3D for Ideal Gaussian Beam
##############################################################
#Write beam block
beamLine        = "Beam"
typeLine        = "Type Gaussian"
fwhmLine        = "FWHM {} {}".format(24.48, 81.96)
fluxLine        = "Flux {flux}".format(flux=noDeconvGaussShift.beamFlux)
energyLine      = "Energy {energy}".format(energy=noDeconvGaussShift.beamEnergy)
collimationLine = "Collimation Rectangular {vert} {horz}".format(horz=noDeconvGaussShift.collimation[0], vert=noDeconvGaussShift.collimation[1])
inputSequence = (beamLine,typeLine,fwhmLine,fluxLine
                ,energyLine,collimationLine)
newline = "\n"
beamBlock = newline.join(inputSequence)

raddoseGaussian = RaddoseBeamRun.RunRaddose(beamBlock, 'IdealGaussianInput.txt', noDeconvGaussShift.doseType)

##############################################################
# Find points close to ideal Gaussian
##############################################################
dataPoints = np.zeros((16,2))
dataPoints[0,0], dataPoints[0,1] = noDeconvGaussThresh.zeroBackgroundPercentage, noDeconvGaussThresh.dose
dataPoints[1,0], dataPoints[1,1] = noDeconvAvDataThresh.zeroBackgroundPercentage, noDeconvAvDataThresh.dose
dataPoints[2,0], dataPoints[2,1] = deconvAvDataWeinerSmoothThresh.zeroBackgroundPercentage, deconvAvDataWeinerSmoothThresh.dose
dataPoints[3,0], dataPoints[3,1] = deconvAvDataGaussSmoothThresh.zeroBackgroundPercentage, deconvAvDataGaussSmoothThresh.dose
dataPoints[4,0], dataPoints[4,1] = deconvAvDataNoSmoothThresh.zeroBackgroundPercentage, deconvAvDataNoSmoothThresh.dose
dataPoints[5,0], dataPoints[5,1] = deconvGaussWeinerSmoothThresh.zeroBackgroundPercentage, deconvGaussWeinerSmoothThresh.dose
dataPoints[6,0], dataPoints[6,1] = deconvGaussGaussSmoothThresh.zeroBackgroundPercentage, deconvGaussGaussSmoothThresh.dose
dataPoints[7,0], dataPoints[7,1] = deconvGaussNoSmoothThresh.zeroBackgroundPercentage, deconvGaussNoSmoothThresh.dose

dataPoints[8,0], dataPoints[8,1] = noDeconvGaussShift.zeroBackgroundPercentage, noDeconvGaussShift.dose
dataPoints[9,0], dataPoints[9,1] = noDeconvAvDataShift.zeroBackgroundPercentage, noDeconvAvDataShift.dose
dataPoints[10,0], dataPoints[10,1] = deconvAvDataWeinerSmoothShift.zeroBackgroundPercentage, deconvAvDataWeinerSmoothShift.dose
dataPoints[11,0], dataPoints[11,1] = deconvAvDataGaussSmoothShift.zeroBackgroundPercentage, deconvAvDataGaussSmoothShift.dose
dataPoints[12,0], dataPoints[12,1] = deconvAvDataNoSmoothShift.zeroBackgroundPercentage, deconvAvDataNoSmoothShift.dose
dataPoints[13,0], dataPoints[13,1] = deconvGaussWeinerSmoothShift.zeroBackgroundPercentage, deconvGaussWeinerSmoothShift.dose
dataPoints[14,0], dataPoints[14,1] = deconvGaussGaussSmoothShift.zeroBackgroundPercentage, deconvGaussGaussSmoothShift.dose
dataPoints[15,0], dataPoints[15,1] = deconvGaussNoSmoothShift.zeroBackgroundPercentage, deconvGaussNoSmoothShift.dose

allDoses = dataPoints[0:,1]
acceptableDoses = [dose for dose in allDoses if dose > 10]
maxAcceptableDose = max(acceptableDoses)
minAcceptableDose = min(acceptableDoses)

#Calculate error
error = (maxAcceptableDose - minAcceptableDose)/((maxAcceptableDose + minAcceptableDose)/2.0) * 100
##############################################################
# Create plot
##############################################################
sns.set_context("talk")
fig = plt.figure(1)
ax = fig.add_subplot(1, 1, 1)
plt.plot(noDeconvGaussThresh.zeroBackgroundPercentage, noDeconvGaussThresh.dose, marker='o', markersize=20, label='NoDecGauThr', linestyle="None")
plt.plot(noDeconvAvDataThresh.zeroBackgroundPercentage, noDeconvAvDataThresh.dose, marker='o', markersize=20, label='NoDecAvThr', linestyle="None")
plt.plot(deconvAvDataWeinerSmoothThresh.zeroBackgroundPercentage, deconvAvDataWeinerSmoothThresh.dose, marker='o', markersize=20, label='DecAvSmWeinThr', linestyle="None")
plt.plot(deconvAvDataGaussSmoothThresh.zeroBackgroundPercentage, deconvAvDataGaussSmoothThresh.dose, marker='o', markersize=20, label='DecAvSmGausThr',linestyle="None")
plt.plot(deconvAvDataNoSmoothThresh.zeroBackgroundPercentage, deconvAvDataNoSmoothThresh.dose, marker='o', markersize=20, label='DecAvSmNoneThr', linestyle="None")
plt.plot(deconvGaussWeinerSmoothThresh.zeroBackgroundPercentage, deconvGaussWeinerSmoothThresh.dose, marker='o', markersize=20, label='DecGauSmWeinThr', linestyle="None")
plt.plot(deconvGaussGaussSmoothThresh.zeroBackgroundPercentage, deconvGaussGaussSmoothThresh.dose, marker='o', markersize=20, label='DecGauSmGausThr', linestyle="None")
plt.plot(deconvGaussNoSmoothThresh.zeroBackgroundPercentage, deconvGaussNoSmoothThresh.dose, marker='o', markersize=20, label='DecGauSmNoneThr', linestyle="None")

plt.plot(noDeconvGaussShift.zeroBackgroundPercentage, noDeconvGaussShift.dose, marker='*', markersize=20, label='NoDecGauShi', linestyle="None")
plt.plot(noDeconvAvDataShift.zeroBackgroundPercentage, noDeconvAvDataShift.dose, marker='*', markersize=20, label='NoDecAvShi', linestyle="None")
plt.plot(deconvAvDataWeinerSmoothShift.zeroBackgroundPercentage, deconvAvDataWeinerSmoothShift.dose, marker='*', markersize=20, label='DecAvSmWeinShi', linestyle="None")
plt.plot(deconvAvDataGaussSmoothShift.zeroBackgroundPercentage, deconvAvDataGaussSmoothShift.dose, marker='*', markersize=20, label='DecAvSmGausShi', linestyle="None")
plt.plot(deconvAvDataNoSmoothShift.zeroBackgroundPercentage, deconvAvDataNoSmoothShift.dose, marker='*', markersize=20, label='DecAvSmNoneShi', linestyle="None")
plt.plot(deconvGaussWeinerSmoothShift.zeroBackgroundPercentage, deconvGaussWeinerSmoothShift.dose, marker='*', markersize=20, label='DecGauSmWeinShi', linestyle="None")
plt.plot(deconvGaussGaussSmoothShift.zeroBackgroundPercentage, deconvGaussGaussSmoothShift.dose, marker='*', markersize=20, label='DecGauSmGausShi', linestyle="None")
plt.plot(deconvGaussNoSmoothShift.zeroBackgroundPercentage, deconvGaussNoSmoothShift.dose, marker='*', markersize=20, label='DecGauSmNoneShi', linestyle="None")

# draw a red hline at dose corresponding to ideal Gaussian that spans the xrange
l = plt.axhline(linewidth=4, color='r', linestyle='--', y=raddoseGaussian.dose)
#Draw translucent bar corresponding to acceptable dose ranges for a given beam
p = plt.axhspan(minAcceptableDose, maxAcceptableDose, facecolor='g', alpha=0.5)

#--------------------------------------------------
#Annotate plot
ax.annotate('Error is ' + str(round(error, 2))+ '%', xy=(20, 11.7), xytext=(5, 9),
            arrowprops=dict(facecolor='black', shrink=0.05))

#--------------------------------------------------
# Change padding and margins, insert legend

fig.tight_layout() #tight margins
leg = plt.legend(loc='lower left', bbox_to_anchor=(1.05, 0))
plt.draw() #to know size of legend

padLeft   = ax.get_position().x0 * fig.get_size_inches()[0]
padBottom = ax.get_position().y0 * fig.get_size_inches()[1]
padTop    = ( 1 - ax.get_position().y0 - ax.get_position().height ) * fig.get_size_inches()[1]
padRight  = ( 1 - ax.get_position().x0 - ax.get_position().width ) * fig.get_size_inches()[0]
dpi       = fig.get_dpi()
padLegend = ax.get_legend().get_frame().get_width() / dpi

widthAx = 3 #inches
heightAx = 3 #inches
widthTot = widthAx+padLeft+padRight+padLegend
heightTot = heightAx+padTop+padBottom

fig.set_size_inches(widthTot,heightTot)
ax.set_position([padLeft/widthTot, padBottom/heightTot, widthAx/widthTot, heightAx/heightTot])

plt.xlim((-5,80))
plt.title("Effect of processing pipeline on dose calculation in RADDOSE-3D")
plt.ylabel(doseType + " (MGy)" )
plt.xlabel("Zero Background (%)")
plt.show()
