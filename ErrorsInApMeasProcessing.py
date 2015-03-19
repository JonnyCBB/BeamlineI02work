import BeamModule
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
# Find max and min dose values
##############################################################


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

# draw a red hline at y=0 that spans the xrange
l = plt.axhline(linewidth=4, color='r', linestyle='--', y=11.88)

p = plt.axhspan(11, 12.5, facecolor='0.5', alpha=0.5)
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
