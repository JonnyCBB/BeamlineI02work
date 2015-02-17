import numpy as np

def generateBeamArrayFromPGM(pgmImageFile):
    """Function that generates a beam array from a PGM file

    INPUTS:
        pgmImageFile     -A string giving the location of the pgm file containing the beam image

    OUTPUTS:
        beamArray        -2D numpy array of ints representing a spatial relative intensity distribution.
    """

    #Print an informative string to the console.
    print "Generating a 2D array from file: \"" + pgmImageFile + "\""

    #Local function variables
    maxValueLine = False     #boolean variable to determine if the line corresponding to the max pixel value has been reached
    pixelCounter = 0     #variable to count the number of pixel values
    beamListCounter = 0     #variable to count the number pixel elements in the pixel value list

    pgmfile = open(pgmImageFile,'r') #Open pgm file for reading

    #Loop through each line in the pgm file
    for line in pgmfile:
        #If the line begins with a 'P' (magic number line) or a '#' (comment line) then do nothing. Else...
        if (line[0] != 'P' and line[0] != '#'):
            #If the line has two values then these represent the pixel width and height
            if len(line.split()) == 2:
                beamArrayWidth = int(line.split()[0])     #extract the width
                beamArrayHeight = int(line.split()[1])     #extract the height
                beamList = np.zeros(beamArrayWidth * beamArrayHeight)     #calculate the total number of pixel values and preallocate an array to store them
                maxValueLine = True     #The next line is the max value line which we can ignore
            elif maxValueLine == True:     #If the current line is the max value line...
                maxValueLine = False    #...then ignore it but set the max value line checker to false for the subsequent lines
            else:
                #The rest of the lines represent pixel values so we want to store these
                beamList[pixelCounter] = int(line)
                pixelCounter += 1     #increment pixel counter

    pgmfile.close()     #close the file

    #preallocate the beam array
    beamArray = np.zeros((beamArrayHeight,beamArrayWidth), dtype=np.int)

    #Loop through each value in the beam array and insert the corresponding pixel value.
    for row in xrange(0,beamArrayHeight):
        for col in xrange(0,beamArrayWidth):
            beamArray[row,col] = beamList[beamListCounter]
            beamListCounter += 1

    print "Successful: Beam has been generated."

    return beamArray     #Return the beam array.

def writePGMFile(beamArray,pgmFileName):
    """Function that writes a plain text PGM file to the specified location.
    Note: The magic number at the header of the file gives the type of PGM.
    In this case the magic number is P2.

    INPUTS:
        beamArray     -A 2D numpy array of integer values as a spatially resolved representation
                        of relative intensities. The values should be between 0 and 255 to be
                        compatible with the .pgm file format.
        pgmFileName   -Name (and location if the full file path is given) of the PGM to be created
                        from the beam array.

    OUTPUTS:
        N/A           -No explicit return value is given but a pgm file is output at the location
                        specified by the user.
    """
    #Create header lines for the pgm file
    magicNumberLine = "P2\n"
    commentLine = "# CREATOR: Diamond Light Source Beamline I02\n"
    arrayHeight, arrayWidth = beamArray.shape     #Get array dimensions
    maxPixelValue = 255     #Get max value of array

    #Print informative string
    print "Writing a PGM file: \"" + pgmFileName + "\" to current directory..."

    #Open a pgm file for writing
    with open(pgmFileName,"w") as pgmFile:
        #Write header lines to pgm file
        pgmFile.write(magicNumberLine)
        pgmFile.write(commentLine)
        pgmFile.write(str(arrayWidth) + " " + str(arrayHeight) + "\n")
        pgmFile.write(str(maxPixelValue) + "\n")
        #Loop through the beam array and insert each element value into the pgm file.
        for row in xrange(0,arrayHeight):
            for col in xrange(0,arrayWidth):
                pgmFile.write(str(beamArray[row,col]) + "\n")

    print "Finished writing PGM file."

############## Main script
beamArray = generateBeamArrayFromPGM('jonny.pgm')
pgmFileName = "jonny_test.pgm"
writePGMFile(beamArray,pgmFileName)
