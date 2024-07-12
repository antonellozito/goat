import Plotter as plotter
import Datahandler as dh

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
plotter.PlotGridCells(datadir, 1)
plotter.PlotGridCellsIterate(datadir, 0)
plotter.PlotVesselPolygon(datadir, -1)

# Optimization history
#---------------------
plotter.PlotShapeOptimizationHistory(datadir, 2)

# Show figures
#-------------
plotter.ShowFigures()
