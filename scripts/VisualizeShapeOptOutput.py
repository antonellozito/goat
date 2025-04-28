from GOATpy import pl as plotter
from GOATpy import dh as dh

# Description
#------------
# Script to visualize optimization output. It is assumed that all necessary
# files have been printed out and are up to date. The data directory in
# which all the files reside has to be defined in 'datadir'

# Data directory
#---------------
# Check if solps is present
datadir = dh.GetDataDirectory()

# Print
print('VisualizeShapeOptOutput: reading from directory: ' + datadir)

# Design
#-------
# Plot the grid
try:
    plotter.PlotGridCells(datadir, 1)
except: 
    print('Could not plot grid cells')

try: 
    plotter.PlotVesselPolygon(datadir, -1)
except: 
    print('Could not plot vessel polygon')

plotter.PlotVesselDisplacement(datadir, 3)
try:
    plotter.PlotVesselDisplacement(datadir, 3)
except:
    print('Could not plot vessel displacement')

# Optimization history
#---------------------
try:
    plotter.PlotShapeOptimizationHistory(datadir, 2)
except: 
    print('Could not plot shape optimization history')

# Show figures
#-------------
plotter.ShowFigures()
