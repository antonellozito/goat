from src import Plotter as plotter
from src import Datahandler as dh
import os
import numpy as np

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
# Simply the current directory

# Print
print('VisualizeGGInput: reading from directory: ' + os.getcwd())

# Inputs
#-------
# Magnetic field
try:
    # Load
    [R, Z, Psi] = dh.ReadRZPsiFile('./rzpsi.dat')

    # Visualize
    resc = 50
    mylevels = np.arange(0, resc, 1)*(np.max(Psi) - np.min(Psi))/resc + np.min(Psi)
    mylevels = [-6., -5, -4., -3., -2., -1., -0.5, -0.1, 0.0, 0.1, 0.5, 1., 2., 3., 4., 5., 6.]
    plotter.PlotStructured2DContourf(R, Z, Psi, 1, levels=mylevels)

except:
    print("Could not load the magnetic field")

# Structure
try:
    # Load 
    structure = dh.ReadStructureFile('structure.dat')

    # Visualize
    plotter.PlotStructure(structure, 1, color='k', linewidth=2)
except: 
    print("Could not load the structure")



# Show figures
#-------------
plotter.ShowFigures()
