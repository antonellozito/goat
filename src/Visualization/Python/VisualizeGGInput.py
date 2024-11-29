import Plotter as plotter
import Datahandler as dh
import os

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
    plotter.PlotStructured2DContourf(R, Z, Psi, 1)

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
