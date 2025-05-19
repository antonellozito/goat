from GOATpy import pl as plotter
from GOATpy import dh as dh
import numpy as np
import os 

# Description
#------------
# Script to visualize grid deformation cost function and constraint 
# data per iterate

# Initialize
#-----------
# figure data
fignum = 0 # counter
maxfignum = 8 # total number of figures to be plotted

# Data directory
#---------------
# Check if solps is present
datadir = dh.GetDataDirectory()
datadir = '/mnt/c/Users/u0110555/Desktop/code_werk/goatf/goatf/Runs/DEMO_P2/output'

# Print
print('VisualizeGDInput: reading from directory: ' + datadir)

# Design
#-------
# Plot the grid
try: 
    for i in np.arange(fignum, maxfignum):
        plotter.PlotGridCellsIterate(datadir, i)

except:
    print("could not plot initial grid")
#plotter.PlotGridCellsFromFile(datadir, 1)
#plotter.PlotGridCellsIterate(datadir, 0)
#plotter.PlotVesselPolygon(datadir, -1)

# Cost function
#--------------
# Length ratio 
try: 
    plotter.PlotLRCostfunctionValueAtVertices(datadir, fignum)
except: 
    print("could not plot length ratio cost function value")
fignum = fignum + 1


# Face angle difference
try: 
    plotter.PlotFADCostfunctionValueAtVertices(datadir, fignum)
except: 
    print("could not plot face angle difference cost function value")
fignum = fignum + 1


# Face angle
try: 
    plotter.PlotFACostfunctionValueAtVertices(datadir, fignum)
except: 
    print("could not plot face angle cost function value")
fignum = fignum + 1


# Length ratio, radial
try: 
    plotter.PlotLRradCostfunctionValueAtVertices(datadir, fignum)
except: 
    print("could not plot radial length ratio cost function value")
fignum = fignum + 1

# Psi ratio, psi based
try: 
    plotter.PlotPRPBCostfunctionValueAtVertices(datadir, fignum)
except: 
    print("could not plot psi ratio cost function value")
fignum = fignum + 1


# Constraints
#------------
# Flux function
try: 
    plotter.PlotFluxFunctionConstraintValueAtVertices(datadir, fignum)
except: 
    print("could not plot flux function constraint value")
fignum = fignum + 1

# Boundary
try: 
    plotter.PlotBoundaryConstraintValueAtVertices(datadir, fignum)
except: 
    print("could not plot boundary constraint value")
fignum = fignum + 1



# orthogonality
try: 
    plotter.PlotOrthogonalityConstraintValueAtVertices(datadir, fignum)
except: 
    print("could not plot orthogonality constraint value")
fignum = fignum + 1

# Show figures
#-------------
plotter.ShowFigures()
