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

# Grid at input
#fullgridfil_input = 'grid_before_GA.dat'
#filepath = datadir + '/' + fullgridfil_input
#grid_input = dh.ReadGAGridDataFile(filepath)

# Grid where error occurred
fullgridfile_in = 'grid_error.dat'

# Reading input grid
#-------------------
gagrid_plot = 1
tria_plot = 0
try : 
    filepath = datadir + '/' + fullgridfile_in
    grid_in = dh.ReadGAGridDataFile(filepath)
except :
    print('No GAGrid found, trying triangular grid')
    gagrid_plot = 0
    filepath = datadir + '/' + 'tria_error.dat'
    try : 
        grid_in = dh.ReadTriaGrid(filepath) 
        tria_plot = 1
    except :
        error('Also no triangular grid found')

# Read vertices of cell
filepath = datadir + '/' + 'vertices_error.dat'
verts = dh.ReadGAIntegerArrayFile(filepath)

# Design
#-------
# Plot the grid
if gagrid_plot == 1:
    #plotter.PlotGridCells(grid_input, 0)
    plotter.PlotGridCellVertArray(grid_in, verts, 1)
    plotter.PlotGridCellsAlignedFaces(grid_in, 2)
    #plotter.PlotGridCellsBoundaryFaces(grid_in,3)
elif tria_plot == 1:
    plotter.PlotTriaCellVertArray(grid_in, verts, 1)



# Legend
#-------
print('Figure 0: grid at input')
print('Figure 1: grid at error')
print('Figure 2: aligned faces')
print('Figure 3: boundary faces')

# Show figures
#-------------
plotter.ShowFigures()


