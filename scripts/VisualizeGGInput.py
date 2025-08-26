from GOATpy import pl as plotter
from GOATpy import dh as dh
import os
import numpy as np
import sys

# Description
#------------
# Script to visualize grid generator input. It is assumed that all necessary
# files have been printed out and are up to date. The data directory in
# which all the files reside has to be defined in the directory above
# 'datadir' (datadir should be the output directory)

# Inputs for the grid generator are the rzpsi.dat file (or in the 
# future the .equ file with the equilibrium data) and the structure.dat
# file with vessel structures. 

# Data directory
#---------------
# Simply the current directory + output directory
datadir = dh.GetDataDirectory()

# Print
print('VisualizeGGInput: reading from directory: ' + os.getcwd())

# Set default names
structurename = 'structure.dat'

# Check if names were given as input (first argument is python script name, 
# second is assumed to be gridname, third is topomesh name)
narg = len(sys.argv)
print("total number of arguments passed: ", narg)
for i in range(1, narg):
    # Check
    if (i == 1):
        structurename = sys.argv[i]

print('VisualizeGGInput: reading structure from file: ' + structurename)

# Inputs
#-------
# Magnetic field
try:
    # Load
    [R, Z, Psi] = dh.ReadRZPsiFile('./rzpsi.dat')

    # Visualize
    resc = 50
    mylevels = np.arange(0, resc, 1)*(np.max(Psi) - np.min(Psi))/resc + np.min(Psi)
    # mylevels = [-6., -5, -4., -3., -2., -1., -0.5, -0.1, 0.0, 0.1, 0.5, 1., 2., 3., 4., 5., 6.]
    plotter.PlotStructured2DContour(R, Z, Psi, 1, levels=mylevels)

except:
    print("Could not load the magnetic field from rzpsi.dat")

# Structure
try:
    # Load 
    structure = dh.ReadStructureFile(structurename)

    # Visualize
    plotter.PlotStructure(structure, 1, linewidth=2)
except: 
    print("Could not load the structure")



# Show figures
#-------------
plotter.ShowFigures()
