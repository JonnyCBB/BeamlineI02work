## A script written by Carina Lobley February 2015
## This script enables the user to produce an input file for RADDOSE 3D immediately after a data collection 

import gda.epics
from time import sleep
import beamLineSpecificEnergy
import time
import os
import errno
import shutil
import java.io.File
from gda.data import NumTracker
from gda.data import PathConstructor 
import scisoftpy as dnp
from gda.analysis.functions import Gaussian
from gda.util import VisitPath
from goToPeak import *
from gda.px.bcm import BCMFinder
from gda.util import QuantityFactory, ElogEntry
finderNameMap = beamline_parameters.FinderNameMapping()


year = time.strftime("%Y")
#today = time.strftime("%Y%m%d")
today = "20150217"
comDir = "/dls/i02/data/2015/cm12149-1/"  # Manually change this as the run changes
                                          # Ideally change this to go to the visit directory and use the copy of the day's startup files


# Code to produce the input file for the crystal block
def crystalBlock():

    # Crystal type, dimensions and orientation

    crystalType = requestInput("Is the crystal Cuboid or Spherical?")
    if "s" in crystalType.lower():
        global crystalType; crystalType = "Spherical"
        global crystalDimensions; crystalDimensions = requestInput("What is the diameter of the crystal in  μm? (e.g. 100)")
    elif "c" in crystalType.lower():
        global crystalType; crystalType = "Cuboid"
        print "A cuboid has 3 dimensions." 
        print "These can be defined with the crystal aligned on the rotation axis and omega rotation of 0 degrees."
        print "The X axis is the vertical dimension on the OAV View - perpendicular to the rotation axis."
        print "The Y axis is horizontal dimension on the OAV View - parallel to the rotation axis."
        print "The Z axis is the dimension going in to the screen in the OAV View - parallel to the beam."
        global crystalDimensions; crystalDimensions = requestInput("What is are the dimensions of the crystal in um? (e.g. 100 100 100)")
        global orientation; orientation = requestInput("Is the crystal sitting parallel and perpendicular to the beam? (Y/N)")
        if "n" in orientation.lower():
            print "With the crystal aligned on the rotation axis and omega rotation of 0 degrees."
            print "AngleP is the angle between the crystal and the rotation axis."
            global angleP; angleP = requestInput("What is AngleP in degrees? (e.g. 10)")
            print "With the crystal aligned on the rotation axis and omega rotation of 90 degrees."
            print "AngleF is the angle between the crystal and the rotation axis."
            global angleF; angleF = requestInput("What is AngleF in degrees? (e.g. 10)")
        else:
            pass
    else:
        print "The crystal must be either Cuboid or Spherical"
        crystalBlock()

    # Absorption coefficient determination
    print "To calculate the absorption coefficient you need:"
    print "- a PDB 4 letted code to use as a model OR"
    print "- the unit cell parameters and content details."
    global absorptionCoeff1; absorptionCoeff1 = requestInput("Do you want to calculate the absorption coefficient?")
    if "n" in absorptionCoeff1.lower():
        global absCoeff1; absCoeff1 = "Average"
    else:
        global absorptionCoeff2; absorptionCoeff2 = requestInput("Do you want to use a published PDB?")
        if "y" in absorptionCoeff2.lower():
            global absCoeff1; absCoeff1 = "EXP"
            global pdbCode; pdbCode = requestInput("Please give the four letter code of the PDB entry?")
            global solHAC; solHAC = requestInput("Please define the heavy atoms in the solvent (e.g. Na 1000 Cl 1000 adds 1M NaCl to the solvent)")
        else:
            global absCoeff1; absCoeff1 = "RDV3"
            print "Please complete a number of details about your crystal contents:"
            global unitCell; unitCell = requestInput("Please enter the unit cell dimensions and angles (e.g. A B C a b c)")
            global nMon; nMon = requestInput("Please enter the integar number of monomers in the unit cell")
            global nRes; nRes = requestInput("Please enter the integar number of amino acids per monomer")
            global nRNA; nRNA = requestInput("Please enter the integar number of RNA nucleotides per monomer")
            global nDNA; nDNA = requestInput("Please enter the integar number of DNA nucleotides per monomer")
            global nPHA; nPHA = requestInput("Please define the heavy atoms in the protein (e.g. S 10 Se 2 adds 10 sulphur and 2 selenium atoms to the monomer)")
            global solHAC; solHAC = requestInput("Please define the heavy atoms in the solvent (e.g. Na 1000 Cl 1000 adds 1M NaCl to the solvent)")
            global solF; solF = requestInput("Please proivde the fraction of the unit cell occupied by solvent")
            print "That completed the crystal definition."


# Code to produce the input file for the beam block
def beamBlock():
   
    # FWHM from morning setup - vertical then horizontal
    xvData=dnp.io.load(comDir+"setup"+"/"+today+"/"+today+"_Beam_profile_y_sma_ap.dat")['miniap_y']
    yvData=dnp.io.load (comDir+"setup"+"/"+today+"/"+today+"_Beam_profile_y_sma_ap.dat")['i_pin']

    xhData=dnp.io.load(comDir+"setup"+"/"+today+"/"+today+"_Beam_profile_x_sma_ap.dat")['miniap_x']
    yhData=dnp.io.load (comDir+"setup"+"/"+today+"/"+today+"_Beam_profile_x_sma_ap.dat")['i_pin']


    def profile_scan_fit(xData, yData):
        dnp.plot.line(xData)

        fr=dnp.fit.fit([ dnp.fit.function.gaussian, dnp.fit.function.offset], xData, yData, [xData.mean(), xData.max(), xData.ptp()*yData.ptp(),0.0],[(xData.min(),xData.max()), (0.0,xData.ptp()), (0.0,xData.ptp()*yData.ptp()),(-1.0,1.0)],ptol=1e-20,optimizer='global')

        return fr
    
    fitResV = profile_scan_fit(xvData, yvData)
    FWHMV=fitResV[1]
    global fwhmvMod; fwhmvMod=(round(FWHMV*1000,2))

    fitResH = profile_scan_fit(xhData, yhData)
    FWHMH=fitResH[1]
    global fwhmhMod; fwhmhMod=(round(FWHMH*1000,2))

    # Flux from the current beamline flux value 
    # This requires the aperture factor and pin diode factors in beamline parameters and the energy conversion lookupTable to be correct
    fluxUsed = flux()
    global fluxMod; fluxMod=(round(fluxUsed,2)) 

    # Energy from the current beamline energy
    energyUsed = BeamLineEnergy_eV()
    global energyMod; energyMod=(round(energyUsed/1000,4))

    # Rectangular collimation from the current beamline s4 slit positions
    s4y = s4_ysize()
    global s4yMod; s4yMod=(round(s4y,2))    
    s4x = s4_xsize()
    global s4xMod; s4xMod=(round(s4x,2))    

# Code to produce the input file for the wedge block
def wedgeBlock:
    global wedgeStart; wedgeStart = requestInput("Please enter the start angle for the data collection")
    global wedgeEnd; wedgeEnd = requestInput("Please enter the end angle for the data collection")

    global expTime; expTime = requestInput("Please enter the total exposure time for the data collection")

# Code to clear up the namespace for the next user

def tidyUp():
    del crystalType
    del crystalDimensions
    del orientation
    del angleP
    del angleF
    del absorptionCoeff1
    del absCoeff1
    del absorptionCoeff2
    del pdbCode
    del solHAC
    del unitCell
    del nMon
    del nRes
    del nRNA
    del nDNA
    del nPHA
    del solHAC
    del solF
    del fwhmvMod
    del fwhmhMod
    del fluxMod
    del energyMod
    del s4yMod
    del s4xMod

# Code to produce the input file
def printInput():
    
    print
    print
    print "Crystal"
    print "Type "+crystalType	
    print "Dimensions "+crystalDimensions
    if "n" in `orientation`.lower():
        print "ANGLEP "+str(angleP)
        print "ANGLEF "+str(angleF)
    else:
        pass 
    print "AbsCoefCalc "+str(absCoeff1)
    if "y" in absorptionCoeff1.lower():
        if "y" in absorptionCoeff2.lower():
            print "PDB "+str(pdbCode)
            print "SOLVENTHEAVYCONC "+str(solHAC)
        else:
            print "UNITCELL "+str(unitCell)
            print "NUMMONOMERS "+str(nMon)
            print "NUMRESIDUES "+str(nRes)
            print "NUMRNA "+str(nRNA)
            print "NUMDNA "+str(nDNA)
            print "PROTEINHEAVYATOMS "+str(nPHA)
            print "SOLVENTHEAVYCONC "+str(solHAC)
            print "SOLVENTFRACTION "+str(solF)
    else:
        pass    
    print "Beam"
    print "Type Gaussian" 			
    print "FWHM "+str(fwhmvMod)+" "+str(fwhmhMod)
    print "Flux "+str(fluxMod)
    print "Energy "+str(energyMod)
    print "Collimation Rectangular "+str(s4yMod)+" "+str(s4xMod)
    print "Wedge "+str(wedgeStart)+""+str(wedgeEnd)
    print "ExposureTime "+str(expTime)
    print
    print
    print "# Please see the User Guide for a number of additional parameters."

crystalBlock()
beamBlock()
wedgeBlock()
printInput()
tidyUp()


