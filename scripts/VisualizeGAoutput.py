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

# Grid after rem_trias
gridfile_inter = 'grid_after_rem_trias.dat'

# Grid after GA
fullgridfile_out = 'grid_after_GA.dat'

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
filepath = datadir + '/' + fullgridfile_out
grid_out = dh.ReadGAGridDataFile(filepath)
print('\n')
print('After GA')
print('---------')
print('Total number of cells: ' + str(grid_out.cell.ntot))
print('Total number of faces: ' + str(grid_out.face.ntot))
print('Total number of vertices: ' + str(grid_out.vert.ntot))

# Design
#-------
# Plot the grid
plotter.PlotGridCells(grid_out,1)

# Intermediate
#--------------
filepath = datadir + '/' + gridfile_inter
grid_int = dh.ReadGAGridDataFile(filepath)
plotter.PlotGridCells(grid_int,2)

# outershell
try: 
    filepath = datadir + '/' + 'outershell.dat'
    tube_rem = dh.ReadGAIntegerArrayFile(filepath)
    plotter.PlotGridCellArray(grid_int, tube_rem, 2)
except:
    print('No outershell data found')

#filepath = datadir + '/' + 'grid_after_53.dat'
#grid3 = dh.ReadGAGridDataFile(filepath)
#plotter.PlotGridCells(grid3,2)

#plotter.PlotGridCellsAlignedFaces(grid_int,3)
#plotter.PlotGridCellsBoundaryFaces(grid,2)

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

# Reading flux tube cells
#------------------------
filepath = datadir + '/' + 'grid_fluxtube.dat'
grid_ft = dh.ReadGAGridDataFile(filepath)
filepath_ftcv = datadir + '/' + 'fluxtubecells.dat'
filepath_ftcv1 = datadir + '/' + 'fluxtubecellsP1.dat'
filepath_ftcvP2 = datadir + '/' + 'fluxtubecellsP2.dat'
ftcv = dh.ReadGAIntegerArrayFile(filepath_ftcv)
ftcvP1 = dh.ReadGAIntegerArrayFile(filepath_ftcv1)
ftcvP2 = dh.ReadGAIntegerArrayFile(filepath_ftcvP2)

# Plot the grid
plotter.PlotGridCellsFtCv(grid_ft, ftcv, ftcvP1, ftcvP2, 3)

# Show figures
#-------------
plotter.ShowFigures()


