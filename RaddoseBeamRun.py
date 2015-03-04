import sys
from subprocess import check_output


class RunRaddose():

    def __init__(self, beamBlock, raddoseInputFile, doseType):
        self.inputFilename = raddoseInputFile
        self.raddose3dpath = "raddose3d.jar"
        self.crystalBlock = self.createCrystalBlock()
        self.beamBlock = self.createBeamBlock(beamBlock)
        self.wedgeBlock = self.createWedgeBlock()
        self.writeRaddose3Dinputfile()
        raddoseLog = self.runRADDOSE3D()
        self.dose = self.parseRaddoseOutput(raddoseLog, doseType)

    def createCrystalBlock(self):
        """Define the crystal block for the RADDOSE-3D input file.
        """
        crystalLine = "Crystal"
        typeLine = "Type Cuboid"
        dimensionsLine = "Dimensions 100 80 60"
        pixPerMicLine = "PixelsPerMicron 0.5"
        coefCalcLine = "CoefCalc exp"
        pdbLine = "pdb 4WUK"
        solventHeavyConcLine = "SolventHeavyConc Na 0"
        inputSequence = (crystalLine,typeLine,dimensionsLine,pixPerMicLine
                        ,coefCalcLine,pdbLine,solventHeavyConcLine)
        newline = "\n"
        crystalBlock = newline.join(inputSequence)
        return crystalBlock

    def createWedgeBlock(self):
        """Define the wedge block for the RADDOSE-3D input file.
        """
        wedgeLine = "Wedge 0 360"
        exposureLine = "ExposureTime 360"
        angResLine = "AngularResolution 1"
        inputSequence = (wedgeLine,exposureLine,angResLine)
        newline = "\n"
        wedgeBlock = newline.join(inputSequence)
        return wedgeBlock

    def createBeamBlock(self, beamBlock):
        """Define the beam block for the RADDOSE-3D input file.
        Note that the block is defined in the Beam module since the
        beams are the only things that differ.
        """
        return beamBlock

    def writeRaddose3Dinputfile(self):
        """Using the Crystal, Beam and Wedge blocks this method writes
        a RADDOSE-3D input text file.
        """
        newline = "\n"
        with open(self.inputFilename,"w") as inputFile:
            inputFile.write(self.crystalBlock)
            inputFile.write(newline)
            inputFile.write(newline)
            inputFile.write(self.beamBlock)
            inputFile.write(newline)
            inputFile.write(newline)
            inputFile.write(self.wedgeBlock)

    def runRADDOSE3D(self):
        """Method that runs RADDOSE-3D and returns the output of the program
        as a byte string
        """
        raddoseCommand = "java -jar {raddose3dPath} -i {inputFile}".format(
            raddose3dPath=self.raddose3dpath,
            inputFile=self.inputFilename)
        raddoseOutput = check_output(raddoseCommand, shell=True)
        return raddoseOutput

    def parseRaddoseOutput(self,raddoseOutput,doseType):
        """Parses the RADDOSE-3D log and returns the specified dose value.
        """
        #Check which operating system is being used so that the correct
        #line endings are being used
        if "win" in sys.platform:
            stringSeparator = "\r\n"
        else:
            stringSeparator = "\n"

        outputLines = raddoseOutput.split(stringSeparator)

        if "DWD" in doseType:
            parseDose = "Average Diffraction Weighted Dose"
        elif "Max" in doseType:
            parseDose = "Max Dose"
        elif "Av" in doseType:
            parseDose = "Average Dose (Whole Crystal)"
        else:
            print "Dose type unrecognised!"
            print "Using Diffraction Weighted Dose (DWD)"
            parseDose = "Average Diffraction Weighted Dose"

        for line in outputLines:
            if parseDose in line and "MGy" in line:
                print line
                dose = self.parseRaddoseLine(line)

        return dose

    def parseRaddoseLine(self, raddoseLine):
        """Parses a line from the RADDDOSE-3D and returns the laast numerical
        value from the line.
        """
        splitLine = raddoseLine.split(" ")
        for element in splitLine:
            if self.isfloat(element):
                value = float(element)

        return value


    def isfloat(self, value):
        """Checks whether a value can be parsed as a float.
        """
        try:
            float(value)
            return True
        except ValueError:
            return False
