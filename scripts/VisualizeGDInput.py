from GOATpy import pl as plotter
from GOATpy import dh as dh
import numpy as np
import os 

# Description
#------------
# Script to visualize optimization input. It is assumed that all necessary
# files have been printed out and are up to date. The data directory in
# which all the files reside has to be defined in 'datadir'

# Initialize
#-----------
# figure data
fignum = 0 # counter
maxfignum = 14 # total number of figures to be plotted

# Data directory
#---------------
# Check if solps is present
datadir = dh.GetDataDirectory()
inputdir = datadir + '/..'

# Print
print('VisualizeGDInput: reading from directory: ' + datadir)


# Design
#-------
# Plot the grid
try: 
    for i in np.arange(fignum-1, maxfignum):
        plotter.PlotGridCellsFromFile(datadir, i)

except:
    print("could not plot initial grid")
#plotter.PlotGridCellsFromFile(datadir, 1)
#plotter.PlotGridCellsIterate(datadir, 0)
#plotter.PlotVesselPolygon(datadir, -1)

# Cost function
#--------------
# Length ratio 
try: 
    plotter.PlotLRCostfunctionVertexPairs(datadir, fignum)
except: 
    print("could not plot length ratio cost function vertex pairs")
fignum = fignum + 1
try: 
    plotter.PlotLRCostfunctionDesiredBias(datadir, fignum)
except: 
    print("could not plot length ratio desired bias")
fignum = fignum + 1


# Face angle difference
plotter.PlotFADCostfunctionVertexPairs(datadir, fignum )
try: 
    plotter.PlotFADCostfunctionVertexPairs(datadir, fignum )
except: 
    print("could not plot face angle difference cost function pairs of vertex pairs")
fignum = fignum + 1
try: 
    plotter.PlotFADCostfunctionWeights(datadir, fignum)
except: 
    print("could not plot face angle difference cost function weights")
fignum = fignum + 1


# Face angle
try: 
    plotter.PlotFACostfunctionVertexPairs(datadir, fignum)
except: 
    print("could not plot face angle cost function vertex pairs")
fignum = fignum + 1


# Length ratio, radial
try: 
    plotter.PlotLRradCostfunctionVertexPairs(datadir, fignum)
except: 
    print("could not plot radial length ratio cost function vertex pairs")
fignum = fignum + 1

# Psi ratio
try: 
    plotter.PlotPRPBCostfunctionVertices(datadir, fignum)
except: 
    print("could not plot psi ratio cost function vertices")
fignum = fignum + 1



# Constraints
#------------
# Vertices constrained by boundary constraints
try:
    plotter.PlotBoundaryConstraintVertices(datadir, fignum)
except:
    # don't do anything
    print("could not plot boundary constraint vertices")
fignum = fignum + 1

# Flux function
try:
    plotter.PlotFluxfunctionConstraintVertices(datadir, fignum)
except:
    print("could not plot flux function constraint vertices")
fignum = fignum + 1

# Fixed flux function value
try:
    plotter.PlotFixedFluxFunctionConstraintVertices(datadir, fignum)
except:
    print("could not plot fixed flux function constraint vertices")
fignum = fignum + 1

# X-points
try: 
    plotter.PlotXPointConstraintVertices(datadir, fignum)
except:
    print("could not plot xpoint constraint vertices")
fignum = fignum + 1

# Orthogonality
try: 
    plotter.PlotOrthogonalityConstraintEdges(datadir, fignum)
except: 
    print("could not plot orthogonality constraint vertices")
fignum = fignum + 1

# Edge lengths
try: 
    plotter.PlotEdgelengthsConstraintEdges(datadir, fignum)
except: 
    print("could not plot edge length constraint vertices")
fignum = fignum + 1

# Line folding
try:
    plotter.PlotLinefoldingConstraintEdges(datadir, fignum)
except: 
    print("could not plot linefolding constraint vertices")
fignum = fignum + 1

# Show figures
#-------------
plotter.ShowFigures()
