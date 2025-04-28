from GOATpy import pl as plotter
from GOATpy import dh as dh

# Description
#------------
# Script to visualize optimization input. It is assumed that all necessary
# files have been printed out and are up to date. The data directory in
# which all the files reside has to be defined in 'datadir'

# Data directory
#---------------
# Check if solps is present
datadir = dh.GetDataDirectory()
#datadir = './goatf/Runs/adjoint_shape_case/compass_adjoint_shape/run/output'

# Print
print('VisualizeShapeOptInput: reading from directory: ' + datadir)

# Design
#-------
# Plot the grid
try:
    plotter.PlotGridCells(datadir, 1)
except:
    print('Could not plot grid cells')
try:
    plotter.PlotInitialVesselPolygon(datadir, -1)
except: 
    print('Could not plot initial vessel')

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
