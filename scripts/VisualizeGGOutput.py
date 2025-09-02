from GOATpy import pl as plotter
from GOATpy import dh as dh
import os
import sys 
import numpy as np
from matplotlib import pyplot as plt

# Description
#------------
# Script to visualize grid generator output. It is assumed that all necessary
# files have been printed out and are up to date. The data directory in
# which all the files reside has to be defined in the 'datadir' directory

# Outputs of the grid generator are only the grid, but we also plot the
# topological mesh, face labels, cell regions, etc

# Data directory
#---------------
# Check if solps is present
datadir = dh.GetDataDirectory()

# Print
print('VisualizeGGOutput: reading output from directory: ' + datadir)
print('VisualizeGGOutput: reading grid from directory: ' + os.getcwd())

# Set default names
topomeshname = 'topomesh.dat'
gridname = 'traduit.out.b2us'
structurename = 'structure.dat'

# Check if names were given as input (first argument is python script name, 
# second is assumed to be gridname, third is topomesh name)
narg = len(sys.argv)
print("total number of arguments passed: ", narg)
for i in range(1, narg):
    # Check
    if (i == 1):
        gridname = sys.argv[i]
    elif (i == 2): 
        topomeshname = sys.argv[i]
    elif (i == 3): 
        structurename = sys.argv[i]

print('VisualizeGGOutput: reading grid from file: ' + gridname)
print('VisualizeGGOutput: reading topomesh from file: ' + topomeshname)
print('VisualizeGGOutput: reading structure from file: ' + structurename)


# Outputs
#-------
try: 
    l1 = dh.GetPolygonCoordinates(datadir + '/' + 'lpgraph_edges.dat')
    plotter.PlotPolygons2D(l1[:, 0], l1[:, 1], 3, color='r')
except:
    pass 
    
# Topological mesh
try:
    # Load
    topomesh = dh.ReadTopomeshFile(datadir + '/' + topomeshname)

    # Visualize
    plotter.PlotTopologicalMesh(topomesh, 1)
    thisaxes = plt.gca()
    thisaxes.set_title('Topological mesh')
    thisaxes.set_xlabel('x [m]')
    thisaxes.set_ylabel('y [m]')

except:
    print("Could not load the topological mesh")

# Grid (no labels)
try:
    # Load (normally in folder above)
    simgrid = dh.ReadTraduitOutB2us(gridname)

    # Visualize
    plotter.PlotGridCells(simgrid, 2)
    thisaxes = plt.gca()
    thisaxes.set_title('Grid cells')
    thisaxes.set_xlabel('x [m]')
    thisaxes.set_ylabel('y [m]')

except: 
    print("Could not load the grid")

# Structure
try: 
    # Load 
    structure = dh.ReadStructureFile(structurename)

    # Visualize
    plotter.PlotStructure(structure, 2, linewidth=2)
    thisaxes = plt.gca()
    thisaxes.set_title('Vessel structures')
    thisaxes.set_xlabel('x [m]')
    thisaxes.set_ylabel('y [m]')
except:
    print("Could not plot structure.dat")

# Face labels
try:
    # At this point, the simulation grid has to be loaded or it will not 
    # work anyway

    # Visualize
    plotter.PlotGridCells(simgrid, 3)
    plotter.PlotGridFaceLabels(simgrid, 3)
    thisaxes = plt.gca()
    thisaxes.set_title('Face labels')
    thisaxes.set_xlabel('x [m]')
    thisaxes.set_ylabel('y [m]')

except: 
    print("Could not load the grid")

# Cell regions
try:
    # At this point, the simulation grid has to be loaded or it will not 
    # work anyway

    # Visualize
    plotter.PlotGridCells(simgrid, 4)
    plotter.PlotCellBasedQuantity2D(simgrid, simgrid.cell.ft, 4)
    thisaxes = plt.gca()
    thisaxes.set_title('Cell flux tube numbers')
    thisaxes.set_xlabel('x [m]')
    thisaxes.set_ylabel('y [m]')

except: 
    print("Could not load the grid")

try: 
    plotter.PlotGridCells(simgrid, 5)
    vals = dh.GetGeneral2DSurfaceData(datadir + '/' + 'gg_vd_radialdensityfunction.dat')
    plotter.PlotGeneral2DContourf(vals[:, 0], vals[:, 1], vals[:, 2], 5)
    thisaxes = plt.gca()
    thisaxes.set_title('GG radial density function [# m^-1]')
    thisaxes.set_xlabel('x [m]')
    thisaxes.set_ylabel('y [m]')
except:
    print("Could not plot radial density distribution function")

try: 
    plotter.PlotGridCells(simgrid, 6)
    vals = dh.GetGeneral2DSurfaceData(datadir + '/' + 'Lmaxpolref.dat')
    plotter.PlotGeneral2DContourf(vals[:, 0], vals[:, 1], vals[:, 2], 6)
    thisaxes = plt.gca()
    thisaxes.set_title('GG maximal poloidal length [m] (without hard bounds)')
    thisaxes.set_xlabel('x [m]')
    thisaxes.set_ylabel('y [m]')
except:
    print("Could not maximal poloidal length distribution function")

try: 
    plotter.PlotGridCells(simgrid, 7)
    vals = dh.GetGeneral2DSurfaceData(datadir + '/' + 'Lminpolref.dat')
    plotter.PlotGeneral2DContourf(vals[:, 0], vals[:, 1], vals[:, 2], 7)
    thisaxes = plt.gca()
    thisaxes.set_title('GG minimal poloidal length [m] (without hard bounds)')
    thisaxes.set_xlabel('x [m]')
    thisaxes.set_ylabel('y [m]')
except:
    print("Could not minimal poloidal length distribution function")

    
# Faces
try: 
    plotter.PlotGridFaces(simgrid, 8)
except:
    print("Could not plot grid faces")

# Surface area
try: 
    plotter.PlotGridCellArea(simgrid, 9, doinverse=True)
    plotter.PlotGridCells(simgrid, 9)
except:
    print("Could not plot grid cell area")

# Face intersections (if present)
try:
    coord = dh.GetVertexCoordinates(datadir + '/' + 'gg_faceintersections.dat')
except:
    print("Could not load intersection data")
try: 
    plotter.PlotGridCells(simgrid, 10)
    plotter.PlotPoints2D(coord[:, 0], coord[:, 1], 10, color='r', marker='o')
except:
    print("Could not plot intersection data")


# Face regions
try:
    # At this point, the simulation grid has to be loaded or it will not 
    # work anyway

    # Visualize
    plotter.PlotGridCells(simgrid, 11)
    plotter.PlotGridFaceRegions(simgrid, 11)
    thisaxes = plt.gca()
    thisaxes.set_title('Face regions')
    thisaxes.set_xlabel('x [m]')
    thisaxes.set_ylabel('y [m]')

except: 
    print("Could not plot face regions")

# Cell regions
try:
    # At this point, the simulation grid has to be loaded or it will not 
    # work anyway

    # Visualize
    plotter.PlotGridCells(simgrid, 12)
    plotter.PlotCellBasedQuantity2D(simgrid, simgrid.cell.region, 12)
    thisaxes = plt.gca()
    thisaxes.set_title('Cell regions')
    thisaxes.set_xlabel('x [m]')
    thisaxes.set_ylabel('y [m]')

except: 
    print("Could not plot cell regions")

# Show figures
#-------------
plotter.ShowFigures()
