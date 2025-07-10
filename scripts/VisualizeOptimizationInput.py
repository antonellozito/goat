import Plotter as plotter

# Description
#------------
# Script to visualize optimization input. It is assumed that all necessary
# files have been printed out and are up to date. The data directory in
# which all the files reside has to be defined in 'datadir'

# Data directory
datadir = './goatf/Examples/DEMO'

# Design
#-------
# Plot the grid
plotter.PlotGridCells(datadir, 1)
plotter.PlotGridCellsIterate(datadir, 0)
plotter.PlotVesselPolygon(datadir, -1)

# Constraints
#------------
# Vertices constrained by boundary constraints
try:
    plotter.PlotBoundaryConstraintVertices(datadir, 2)
except:
    # don't do anything
    print("could not plot boundary constraint vertices")

# Flux function
try:
    plotter.PlotFluxfunctionConstraintVertices(datadir, 3)
except:
    print("could not plot flux function constraint vertices")


# X-points
try: 
    plotter.PlotXPointConstraintVertices(datadir, 4)
except:
    print("could not plot xpoint constraint vertices")

# Orthogonality
try: 
    plotter.PlotOrthogonalityConstraintEdges(datadir, 5)
except: 
    print("could not plot orthogonality constraint vertices")

# Edge lengths
try: 
    plotter.PlotEdgelengthsConstraintEdges(datadir, 6)
except: 
    print("could not plot edge length constraint vertices")

# Show figures
#-------------
plotter.ShowFigures()
