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

# Grid before GA
fullgridfile_in = 'grid_before_GA.dat'

# Grid after GA
fullgridfile_out = 'grid_after_GA.dat'

# Reading input grid
#-------------------
filepath = datadir + '/' + fullgridfile_in
grid1 = dh.ReadGAGridDataFile(filepath)

print('Before GA')
print('---------')
print('Total number of cells : ' + str(grid1.cell.ntot))
print('Total number of faces: ' + str(grid1.face.ntot))
print('Total number of vertices: ' + str(grid1.vert.ntot))

# Design
#-------
# Plot the grid
plotter.PlotGridCells(grid1,0)

# Reading output grid
#--------------------
filepath = datadir + '/' + fullgridfile_out
grid2 = dh.ReadGAGridDataFile(filepath)
print('\n')
print('After GA')
print('---------')
print('Total number of cells: ' + str(grid2.cell.ntot))
print('Total number of faces: ' + str(grid2.face.ntot))
print('Total number of vertices: ' + str(grid2.vert.ntot))

# Design
#-------
# Plot the grid
plotter.PlotGridCells(grid2,1)
#plotter.PlotGridCellsAlignedFaces(grid,1)
#plotter.PlotGridCellsBoundaryFaces(grid,2)

# Reading cutcell arrays
#-----------------------
#filepath = datadir + '/' + 'cctria.dat'
#cctria = dh.ReadGAArrayFile(filepath)
#filepath = datadir + '/' + 'cctrapsP1.dat'
#cctrapsP1 = dh.ReadGAArrayFile(filepath)
#filepath = datadir + '/' + 'cctrapsP2.dat'
#cctrapsP2 = dh.ReadGAArrayFile(filepath)
#filepath = datadir + '/' + 'cctraps.dat'
#cctraps = dh.ReadGAArrayFile(filepath)

#print(str(cctria[0]))
#plotter.PlotGridCellCutcells(grid2,cctria,cctrapsP1,cctrapsP2,cctraps, 2)



# Show figures
#-------------
plotter.ShowFigures()


