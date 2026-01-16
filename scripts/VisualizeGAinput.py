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
print('VisualizeGAintput: reading from directory: ' + datadir)

# Grid before GA
fullgridfile_in = 'grid_before_GA.dat'
fullgridfile_fs = 'grid_fluxsurface.dat'

# Reading input grid
#-------------------
filepath = datadir + '/' + fullgridfile_in
grid_in = dh.ReadGAGridDataFile(filepath)

print('Before GA')
print('---------')
print('Total number of cells : ' + str(grid_in.cell.ntot))
print('Total number of faces: ' + str(grid_in.face.ntot))
print('Total number of vertices: ' + str(grid_in.vert.ntot))

# Design
#-------
# Plot the grid
plotter.PlotGridCells(grid_in,0)

# Reading output grid
#--------------------
filepath = datadir + '/' + fullgridfile_fs
grid_fs = dh.ReadGAGridDataFile(filepath)
filepath_fsfc = datadir + '/' + 'fluxsurfacefaces.dat'
filepath_fsfcP1 = datadir + '/' + 'fluxsurfacefacesP1.dat'
filepath_fsfcP2 = datadir + '/' + 'fluxsurfacefacesP2.dat'
fsfc = dh.ReadGAIntegerArrayFile(filepath_fsfc)
fsfcP1 = dh.ReadGAIntegerArrayFile(filepath_fsfcP1)
fsfcP2 = dh.ReadGAIntegerArrayFile(filepath_fsfcP2)

# Design
#-------
# Plot the grid
plotter.PlotGridCellsFsFc(grid_fs, fsfc, fsfcP1, fsfcP2, 1)

plotter.PlotGridCellsAlignedFaces(grid_in,2)
#plotter.PlotGridCellsBoundaryFaces(grid,3)

# Reading cutcell arrays
#-----------------------
#filepath = datadir + '/' + 'cctria.dat'
#cctria = dh.ReadGAIntegerArrayFile(filepath)
#filepath = datadir + '/' + 'cctrapsP1.dat'
#cctrapsP1 = dh.ReadGAIntegerArrayFile(filepath)
#filepath = datadir + '/' + 'cctrapsP2.dat'
#cctrapsP2 = dh.ReadGAIntegerArrayFile(filepath)
#filepath = datadir + '/' + 'cctraps.dat'
#cctraps = dh.ReadGAIntegerArrayFile(filepath)

#print(str(cctria[0]))
#plotter.PlotGridCellCutcells(grid2,cctria,cctrapsP1,cctrapsP2,cctraps, 2)

# Reading farSOL interpolation value
#-----------------------------------
#filepath = datadir + '/' + 'farSOLint.dat'
#farSOLint = dh.ReadGARealArrayFile(filepath)
#plotter.PlotGridCellValue(grid1,farSOLint, 0.95, 2)

# Legend
#-------
print('Figure 0: grid at input')
print('Figure 1: flux surfaces')
print('Figure 2: aligned faces')
#print('Figure 3: boundary faces')

# Show figures
#-------------
plotter.ShowFigures()


