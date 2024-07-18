import Plotter as plotter
import Datahandler as dh

# Description
#------------
# Script to visualize optimization input. It is assumed that all necessary
# files have been printed out and are up to date. The data directory in
# which all the files reside has to be defined in 'datadir'

# Data directory
#---------------
# Check if solps is present
datadir = dh.GetDataDirectory()

# Print
print('VisualizeShapeOptInput: reading from directory: ' + datadir)

# Design
#-------
# Plot the grid
plotter.PlotGridCells(datadir, 1)
plotter.PlotGridCellsIterate(datadir, 0)
plotter.PlotVesselPolygon(datadir, -1)

# Constraints
#------------
# Vertices constrained by fixed points
try:
    plotter.PlotFixedVesselPointsConstraintVertices(datadir, 2)
except:
    # don't do anything
    print("could not plot fixed vessel points constraint vertices")

# Vertices constrained by fixed flux constraints
try:
    plotter.PlotFixedVesselFluxConstraintVertices(datadir, 3)
except:
    print("could not plot fixed vessel flux constraint vertices")

# Show figures
#-------------
plotter.ShowFigures()
