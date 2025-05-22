from GOATpy import pl as plotter
from GOATpy import dh as dh
import os
import sys 
import numpy as np

# Description
#------------
# Script to visualize the output of b2ag (i.e. the generated b2fgmtry
# file with metrics etc). It is assumed that the b2fgmtry is 
# located in the current folder. 

# Data directory
#---------------
# Check if solps is present
datadir = dh.GetDataDirectory()

# Print
dir = os.getcwd()
print('VisualizeB2agOutput: reading b2fgmtry from directory: ' + dir)

# Set default names
gridname = 'b2fgmtry'

# Check if names were given as input 
narg = len(sys.argv)
print("total number of arguments passed: ", narg)
for i in range(1, narg):
    # Check
    if (i == 1):
        gridname = sys.argv[i]

griddir = dir + '/' + gridname
print('VisualizeB2agOutput: reading b2fgmtry from file: ' + griddir)

# Outputs
#-------
# Initialize
fignum = 0

# Load (normally in folder above)
simgrid = dh.ReadGridFromB2fgmtryus(griddir)

# Plot face labels
plotter.PlotGridCells(simgrid, fignum)
plotter.PlotGridFaceLabels(simgrid, fignum)

fignum = fignum + 1

# Plot cell regions
plotter.PlotGridCells(simgrid, fignum)
plotter.PlotCellBasedQuantity2D(simgrid, simgrid.cell.region[0:simgrid.cell.ntot], fignum)
fignum = fignum + 1

# Plot magnetic field quantities
plotter.PlotGridCells(simgrid, fignum)
plotter.PlotCellBasedQuantity2D(simgrid, simgrid.cell.bb[:, 0], fignum)
ind = np.argmax(simgrid.cell.bb[:, 0])
plotter.PlotPoints2D(simgrid.cell.x[ind], simgrid.cell.y[ind], fignum, marker='o', color='r')
fignum = fignum + 1

plotter.PlotGridCells(simgrid, fignum)
plotter.PlotCellBasedQuantity2D(simgrid, simgrid.cell.bb[:, 1], fignum)
ind = np.argmax(simgrid.cell.bb[:, 1])
plotter.PlotPoints2D(simgrid.cell.x[ind], simgrid.cell.y[ind], fignum, marker='o', color='r')
fignum = fignum + 1

plotter.PlotGridCells(simgrid, fignum)
plotter.PlotCellBasedQuantity2D(simgrid, simgrid.cell.bb[:, 2], fignum)
ind = np.argmax(simgrid.cell.bb[:, 2])
plotter.PlotPoints2D(simgrid.cell.x[ind], simgrid.cell.y[ind], fignum, marker='o', color='r')
fignum = fignum + 1

plotter.PlotGridCells(simgrid, fignum)
plotter.PlotCellBasedQuantity2D(simgrid, simgrid.cell.bb[:, 3], fignum)
ind = np.argmax(simgrid.cell.bb[:, 3])
plotter.PlotPoints2D(simgrid.cell.x[ind], simgrid.cell.y[ind], fignum, marker='o', color='r')
fignum = fignum + 1




# Show figures
#-------------
plotter.ShowFigures()
