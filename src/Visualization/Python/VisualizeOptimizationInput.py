import Plotter as plotter

# Description
#------------
# Script to visualize optimization input. It is assumed that all necessary
# files have been printed out and are up to date. The data directory in
# which all the files reside has to be defined in 'datadir'

# Data directory
datadir = '../Plotdata'

# Design
#-------
# Plot the grid
plotter.PlotGridCells(datadir, 1)
plotter.PlotGridCellsIterate(datadir, 0)

# Constraints
#------------
# Vertices constrained by boundary constraints
plotter.PlotBoundaryConstraintVertices(datadir, 2)

# Flux function
plotter.PlotFluxfunctionConstraintVertices(datadir, 3)

# X-points
plotter.PlotXPointConstraintVertices(datadir, 4)

# Orthogonality
plotter.PlotOrthogonalityConstraintEdges(datadir, 5)

# Show figures
#-------------
plotter.ShowFigures()
