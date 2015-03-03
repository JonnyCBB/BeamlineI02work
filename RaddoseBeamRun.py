import BeamModule

class RunRaddose():

    def __init__(self, beam, raddoseInputFile):
        self.crystalBlock = self.createCrystalBlock()
        self.beamBlock = self.createBeamBlock(beam)
        self.wedgeBlock = self.createWedgeBlock()
        self.writeRaddose3Dinputfile(raddoseInputFile)

    def createCrystalBlock(self):
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
        wedgeLine = "Wedge 0 360"
        exposureLine = "ExposureTime 360"
        angResLine = "AngularResolution 1"
        inputSequence = (wedgeLine,exposureLine,angResLine)
        newline = "\n"
        wedgeBlock = newline.join(inputSequence)
        return wedgeBlock

    def createBeamBlock(self, beam):
        beamBlock = beam.beamBlock
        return beamBlock

    def writeRaddose3Dinputfile(self,raddoseFilename):
        newline = "\n"
        with open(raddoseFilename,"w") as inputFile:
            inputFile.write(self.crystalBlock)
            inputFile.write(newline)
            inputFile.write(newline)
            inputFile.write(self.beamBlock)
            inputFile.write(newline)
            inputFile.write(newline)
            inputFile.write(self.wedgeBlock)
