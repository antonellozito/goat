from src import Plotter as plotter
from src import Datahandler as dh
import os

# Description
#------------
# Script to visualize grid generator output. It is assumed that all necessary
# files have been printed out and are up to date. The data directory in
# which all the files reside has to be defined in the 'datadir' directory

# Outputs of the grid generator are only the grid, but we also plot the
# topological mesh, face labels, cell regions, etc

# Data directory
#---------------
# Check if solps is present
datadir = dh.GetDataDirectory()

# Print
print('VisualizeGGOutput: reading output from directory: ' + datadir)
print('VisualizeGGOutput: reading grid from directory: ' + os.getcwd())

# Outputs
#-------
# Topological mesh
try:
    # Load
    topomesh = dh.ReadTopomeshFile(datadir + '/topomesh.dat')

    # Visualize
    plotter.PlotTopologicalMesh(topomesh, 1)

except:
    print("Could not load the topological mesh")

# Grid (no labels)
try:
    # Load (normally in folder above)
    simgrid = dh.ReadTraduitOutB2us('traduit.out.b2us')

    # Visualize
    plotter.PlotGridCells(simgrid, 2)

except: 
    print("Could not load the grid")

# Face labels
try:
    # At this point, the simulation grid has to be loaded or it will not 
    # work anyway

    # Visualize
    plotter.PlotGridCells(simgrid, 3)
    plotter.PlotGridFaceLabels(simgrid, 3)

except: 
    print("Could not load the grid")

# Cell regions
try:
    # At this point, the simulation grid has to be loaded or it will not 
    # work anyway

    # Visualize
    plotter.PlotGridCells(simgrid, 4)
    plotter.PlotCellBasedQuantity2D(simgrid, simgrid.cell.region, 4)

except: 
    print("Could not load the grid")



# Show figures
#-------------
plotter.ShowFigures()
