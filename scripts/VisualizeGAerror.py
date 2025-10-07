from GOATpy import pl as plotter
from GOATpy import dh as dh

# Description
#-------------
# Script to visualize grid adaptation output. 
# It is assumed that all necessary
# files have been printed out and are up to date.

# Data directory
#---------------
# Check if solps is present
datadir = dh.GetDataDirectory()

# Print
print('VisualizeGAerror: reading from directory: ' + datadir)

# Grid where error occurred
fullgridfile_in = 'grid_error.dat'

# Reading input grid
#-------------------
filepath = datadir + '/' + fullgridfile_in
grid_in = dh.ReadGAGridDataFile(filepath)

# Read vertices of cell
filepath = datadir + '/' + 'vertices_error.dat'
verts = dh.ReadGAIntegerArrayFile(filepath)

# Design
#-------
# Plot the grid
plotter.PlotGridCellVertArray(grid_in, verts, 0)
plotter.PlotGridCellsAlignedFaces(grid_in,1)
plotter.PlotGridCellsBoundaryFaces(grid_in,2)

# Show figures
#-------------
plotter.ShowFigures()


