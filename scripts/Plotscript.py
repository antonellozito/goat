# Description
#------------
# Simple plot script that uses the other python functionality
import sys
import os
import numpy as np
from matplotlib import pyplot as plt

script_dir = os.path.dirname(os.path.abspath(sys.argv[0]))  # Get the script's directory
sys.path.append(script_dir)  # Add Visualization/ to sys.path


import GOATpy as gp

# Set base_dir, should not be changed
bd = './'

# Get all subdirectories in base_dir
subfolders = [f for f in os.listdir(bd) if os.path.isdir(os.path.join(bd, f))]

# Check if there are any subfolders
if not subfolders:
    print("No subfolders found in Runs/. Exiting.")
    exit()

# Print the subfolders with numbers
print("Select a folder:")
for idx, folder in enumerate(subfolders, start=1):
    print(f"{idx}) {folder}")

# User input selection
while True:
    try:
        choice = int(input("Enter the number of the folder: "))
        if 1 <= choice <= len(subfolders):
            sd = subfolders[choice - 1]
            break
        else:
            print("Invalid number. Please enter a valid option.")
    except ValueError:
        print("Invalid input. Please enter a number.")

# Print and store the selected folder name
print(f"You selected: {sd}")

# Set base_dir, should not be changed
base_dir = sd

# Get all subdirectories in base_dir
subfolders = [f for f in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, f))]

# Check if there are any subfolders
if not subfolders:
    print("No subfolders found in Runs/. Exiting.")
    exit()

# Print the subfolders with numbers
print("Select a folder:")
for idx, folder in enumerate(subfolders, start=1):
    print(f"{idx}) {folder}")

# User input selection
while True:
    try:
        choice = int(input("Enter the number of the folder: "))
        if 1 <= choice <= len(subfolders):
            s_dir = subfolders[choice - 1]
            break
        else:
            print("Invalid number. Please enter a valid option.")
    except ValueError:
        print("Invalid input. Please enter a number.")

# Print and store the selected folder name
print(f"You selected: {s_dir}")

simdir = './Runs/' + s_dir
outputdir = simdir + '/output'
structuredir = './goatf/src/Visualization'
writedir = './goatf/src/Visualization'

try: 
    l1 = gp.dh.GetPolygonCoordinates(outputdir + '/l1.dat')
    l2 = gp.dh.GetPolygonCoordinates(outputdir + '/l2.dat')
    l3 = gp.dh.GetPolygonCoordinates(outputdir + '/l3.dat')
    l4 = gp.dh.GetPolygonCoordinates(outputdir + '/l4.dat')
    k1 = 1
    k2 = 1
    v11 = k1 - 1
    v12 = k1 
    v13 = k1+1
    v21 = k2 - 1
    v22 = k2 
    v23 = k2+1
    gp.pl.PlotPolygons2D(l1[:, 0], l1[:, 1], 1, color='r', marker='o')
    gp.pl.PlotPolygons2D(l2[:, 0], l2[:, 1], 1, color='g', marker='o')
    gp.pl.PlotPolygons2D(l3[:, 0], l3[:, 1], 1, color='b', marker='o')
    gp.pl.PlotPolygons2D(l4[:, 0], l4[:, 1], 1, color='m', marker='o')


    #x1 = [4.4388371941530949, 4.4879558708267488]
    #y1 = [-3.7138585019166284, -3.6936360969799154]
    #x2 = [4.3916720661976889, 4.4426284461966903]
    #y2 = [-3.7435296275342731, -3.7230237727210027]
    x1 = l1[[v11, v12, v13], 0]
    y1 = l1[[v11, v12, v13], 1]
    x2 = l2[[v21, v22, v23], 0]
    y2 = l2[[v21, v22, v23], 1]

    #simgrid = dh.ReadTraduitOutB2us(simdir + '/traduit.out.b2us')
    #gp.pl.PlotGridCells(simgrid, 0)

    gp.pl.PlotPolygons2D(x1, y1, 1, color='r', marker='o')
    gp.pl.PlotPolygons2D(x2, y2, 1, color='g', marker='o')
    gp.pl.PlotPolygons2D([x1[0], x2[0]], [y1[0], y2[0]], 1, color='k') # previous face
    gp.pl.PlotPolygons2D([x1[0], x2[1]], [y1[0], y2[1]], 1, color='b') # triangle 1
    gp.pl.PlotPolygons2D([x1[1], x2[0]], [y1[1], y2[0]], 1, color='m') # triangle 2
    gp.pl.PlotPolygons2D([x1[1], x2[1]], [y1[1], y2[1]], 1, color='k', linestyle='--') # quad
    gp.pl.PlotPolygons2D([x1[1], x2[2]], [y1[1], y2[2]], 1, color='b', linestyle='--') # next tria 1
    gp.pl.PlotPolygons2D([x1[2], x2[1]], [y1[2], y2[1]], 1, color='m', linestyle='--') # next tria 2
    gp.pl.PlotPoints2DWithID(x1, y1, [11, 12], 1)
    gp.pl.PlotPoints2DWithID(x2, y2, [21, 22], 1)
    gp.pl.SetAxesLimits2D(plt.gca(), l1[:, 0], l1[:, 1])
    #gp.pl.SetAxesLimits2D(plt.gca(), x1, y1)
    #gp.pl.ShowFigures()
except:
    print("Could not plot vessel vertices.")

# Plot stuff
#-----------
#[R, Z, Psi] = dh.ReadRZPsiFile(simdir + '/rzpsi.dat')
structure = gp.dh.ReadStructureFile(simdir + '/structure.dat')
#dh.WriteRZPsiFile(writedir, R, Z, Psi)
#gp.pl.PlotStructured2DContourf(R, Z, Psi, 3)
gp.pl.PlotStructure(structure, 3)
#gp.pl.PlotPolygonData(outputdir + '/temppol.dat', 5, color = 'r')
#structure = dh.ReadStructureFile(structuredir + '/structure.dat')
#gp.pl.PlotStructure(structure, 2)

try:
    simgrid = gp.dh.ReadTraduitOutB2us(simdir + '/traduit.out.b2us')
    #tf = 31
    #tfind = tf - 1
    #gp.pl.PlotPoints2D(simgrid.vert.x[[simgrid.face.v1[tfind]-1, simgrid.face.v2[tfind]-1]], 
    #    simgrid.vert.y[[simgrid.face.v1[tfind]-1, simgrid.face.v2[tfind]-1]], 0, 
    #    color='r', marker='o', linewidth=3)
    
    gp.pl.PlotGridCells(simgrid, 0)
    gp.pl.PlotGridFaceLabels(simgrid, 0)
    # gp.pl.PlotGridVertFieldlineID(simgrid, 0)
    
    #gp.pl.PlotGridCells(simgrid, 1)
    gp.pl.PlotGridCells(simgrid, 3)
    gp.pl.PlotGridCells(simgrid, -1)
    gp.pl.PlotGridFaceRegions(simgrid, -1)
    #gp.pl.PlotCellBasedQuantity2D(simgrid, simgrid.cell.ft, 2)
    #gp.pl.PlotCellBasedQuantity2D(simgrid, simgrid.cell.cflags, 3)
    gp.pl.PlotCellBasedQuantity2D(simgrid, simgrid.cell.bt, 14)
    gp.pl.PlotCellBasedQuantity2D(simgrid, simgrid.cell.bp, 15)
    #gp.pl.PlotCellBasedQuantity2D(simgrid, simgrid.cell.region, 4)
except: 
    print("Could not plot simulation grid.")

try:
    topomesh = gp.dh.ReadTopomeshFile(outputdir + '/topomesh.dat')
    gp.pl.PlotTopologicalMesh(topomesh, 5)
    
    gp.pl.PlotTopologicalMesh(topomesh, 0)
    gp.pl.PlotTopologicalMesh(topomesh, 1)
    gp.pl.PlotTopologicalMesh(topomesh, 2)
    gp.pl.PlotTopologicalMesh(topomesh, 3)
    gp.pl.PlotTMCellBasedQuantity(topomesh, np.arange(0, topomesh.cell.ntot), 7)
except:
    print("Could not plot final topological mesh.")

try: 
    topomesh_base = gp.dh.ReadTopomeshFile(outputdir + '/topomesh_base.dat')
    gp.pl.PlotTopologicalMesh(topomesh_base, 6)
except:
    print("Could not plot base topological mesh.")

try:
    print('Ciao')
    ggtmdata = gp.dh.ReadGGTMDataFile(outputdir + '/ggtmdata_after_vertexdistribution.dat')
    grid = gp.dh.ReadGGGridDataFile(outputdir + '/grid_after_cellconstruction.dat')

    gp.pl.PlotGGTMDataFaceVertexDistribution(ggtmdata, topomesh, -2)
    gp.pl.PlotGGTMDataCellVertexDistribution(ggtmdata, topomesh, -2)
    gp.pl.PlotGridVertices(grid, 9)
    gp.pl.PlotGridFaces(grid, 9)
    gp.pl.PlotTopologicalMesh(topomesh, 9)
    print('Fine ciao')
except:
    print("Could not plot intermediate data.")

# Diagnostics
try: 
    # Intersecting faces
    faceIDs, vals = gp.dh.GetVertexCoordinatesWithID(outputdir + '/gg_faceintersections.dat')
    gp.pl.PlotPoints2D(vals[:, 0], vals[:, 1], 9, color='g', marker='o')
except:
    print("Could not plot diagnostics.")


#gp.pl.PlotGridVertices(grid, 2)
#gp.pl.PlotGridCells(grid, 3)
#gp.pl.PlotPolygonData(outputdir + '/vesselpolygon.dat', 3, color='r', linewidth=0.5)
#gp.pl.PlotPolygonData(outputdir + '/vesselpolygon.dat', 1, color='r', linewidth=0.5)

#gp.pl.PlotGridVertices(grid, 1)
#gp.pl.PlotTopologicalMeshCells(topomesh, 1)
#gp.pl.PlotPolygonData(outputdir + '/extrema_fx_lines.dat', 2, color='r')
#gp.pl.PlotPolygonData(outputdir + '/extrema_fy_lines.dat', 2, color='g')
#gp.pl.PlotPolygonData(outputdir + '/mfcontours.dat', 5, color='k')

#mylevels = [-1, -0.1, -0.01, -0.001, -0.0001, -0.000001, 0, 0.000001, 0.0001, 0.001, 0.01, 0.1, 1]
#mylevelsbias = [0.95, 0.96, 0.97, 0.98, 0.99, 1, 1.01, 1.02, 1.03, 1.04, 1.05]
#gp.pl.Plot2DSurfaceDataContourf(outputdir + '/generalplf.dat', 0, levels=mylevels)
#gp.pl.Plot2DSurfaceDataContourf(outputdir + '/vesselplf.dat', 5, levels=mylevels)
#gp.pl.PlotPolygonData(outputdir + '/vesselpolygon.dat', 5, color='m')
# Show figures
#-------------


# gp.pl.PlotPolygonData(outputdir + '/testpolyg.dat', 2)
# gp.pl.Plot2DSurfaceDataContourf(outputdir + '/costfunctionLR_vesselcontours.dat', 2)
#gp.pl.Plot2DSurfaceDataContourf(outputdir + '/constraints_boundary_plf.dat', 2, levels = mylevels)
#gp.pl.PlotBoundaryConstraintVertices(outputdir, 3)

gp.pl.ShowFigures()
