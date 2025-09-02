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
print('VisualizeGAoutput: reading from directory: ' + datadir)
fullgridfile = 'grid_after_GA.dat'


# Reading
#--------
filepath = datadir + '/' + fullgridfile
grid = dh.ReadGAGridDataFile(filepath)

print('Total number of cells: ' + str(grid.cell.ntot))
print('Total number of faces: ' + str(grid.face.ntot))
print('Total number of vertices: ' + str(grid.vert.ntot))

# Design
#-------
# Plot the grid
plotter.PlotGridCells(grid,0)
plotter.PlotGridCellsAlignedFaces(grid,1)
plotter.PlotGridCellsBoundaryFaces(grid,2)


# Show figures
#-------------
plotter.ShowFigures()


