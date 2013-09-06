"""
This script extracts data value names form MIDAS file.
Takes two arguments: input and output file paths.
Input file must be in MIDAS format.
Output file can ve a text file.
Each name will be printed in a new line in output file.

Gungor Budak @gungorbudak
Department of Bioinformatics, BiGCaT, Maastricht University
2013 Summer, Maastricht, Netherlands
"""

import sys
args = sys.argv[1:]

# File definitions
fIn = open(args[0], "r") # input file (MIDAS/.csv file)
fOut = open(args[1], "w") # output file (Text/.txt file)

# Read the first line
line = fIn.readline().strip()

# Split line & make a list of DV names
elements = line.split(",")
dvElements = [ x for x in elements if "DV" in x ]

# Initiate DV Names dictionary
dvNames = {}

# Split protein names
for dvElement in dvElements:
	dvElement = dvElement.split(":")[1]
	dvElement = dvElement.split("_")[0]
	dvNames[dvElement] = 1

# Output string
dvNamesString = ""

# Loop through names
for dvName in dvNames.keys():
	dvNamesString += dvName + "\n"

# Write output string
fOut.write(dvNamesString)


