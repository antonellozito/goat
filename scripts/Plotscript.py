# Description
#------------
# Simple plot script that uses the other python functionality
from src import Plotter as pl 
from src import Datahandler as dh 
from src import goat_types as gt
import numpy as np
from matplotlib import pyplot as plt

# Set paths
#----------
# Simulation directory
#thisdir = './goatf/Runs/TCV_Anthony'
#Rpath = thisdir + '/R.csv'
#Zpath = thisdir + '/Z.csv'
#Psipath = thisdir + '/psi.csv'
#separator = ';'
#[R, Z, PsimWb, nZ, nR] = dh.ReadRZPsiFromCSV(Rpath, Zpath, Psipath, separator)
# [R, Z, PsiWb] = dh.ReadRZPsiFromEqdskFile('./goatf/Runs/DEMO_P2/Equil_DEMO_PROCESS_SOB_COCOS11.eqdsk')
#Psi = PsiWb
#Psi = PsimWb/1000.0 # rescale
#pl.PlotStructured2DContourf(R,Z, Psi, 10, levels=20)
#dh.WriteRZPsiFile('./goatf/Runs/TCV_Anthony', R, Z, Psi)

simdir = dh.GetDataDirectory()
#simdir = './goatf/Examples/JT60SA/5degreesRP'
simdir = './goatf/Runs/JT60SA/Giulio'
#simdir = './goatf/Runs/DEMO_P2'
outputdir = simdir + '/output'
structuredir = './goatf/src/Visualization'
writedir = './goatf/src/Visualization'
griddir = 'traduit.out.b2us'
structurefile = 'structure.dat.dat'
simgrid = dh.ReadTraduitOutB2us(simdir + '/' + griddir)

#mylevels = [0.001, 0.005, 0.01, 0.02, 0.05, 0.1]
#mylevels = [0.01, 0.02, 0.05, 0.1, 0.2, 0.5]
#mylevelsmf = [-10, -8, -6, -4, -2, -1, 0, 1, 2, 4, 8, 10]
mylevelsdf = 20
#mylevelsdf = 2*mylevelsdf
#pl.Plot2DSurfaceDataContourf(outputdir + '/Lmaxradref.dat', 0, levels=mylevels)
#pl.Plot2DSurfaceDataContourf(outputdir + '/magneticfield_visualization.dat', 11, levels=mylevelsmf)
#simgrid = dh.ReadTraduitOutB2us(simdir + '/' + griddir)
#pl.PlotGridTopologicalData(simgrid, 17)

#thisinterp = gt.GridInterpolant2D()
#thisinterp.Construct(simgrid, simgrid.cell.region, 'cartesian')

#pl.VisualizeGridInterpolant2D(thisinterp, -5)

try: 
    l1 = dh.GetPolygonCoordinates(outputdir + '/l1.dat')
    l2 = dh.GetPolygonCoordinates(outputdir + '/l2.dat')
    l3 = dh.GetPolygonCoordinates(outputdir + '/l3.dat')
    l4 = dh.GetPolygonCoordinates(outputdir + '/l4.dat')
    k1 = 1
    k2 = 1
    v11 = k1 - 1
    v12 = k1 
    v13 = k1+1
    v21 = k2 - 1
    v22 = k2 
    v23 = k2+1
    pl.PlotPolygons2D(l1[:, 0], l1[:, 1], 1, color='r', marker='o')
    pl.PlotPolygons2D(l2[:, 0], l2[:, 1], 1, color='g', marker='o')
    pl.PlotPolygons2D(l3[:, 0], l3[:, 1], 1, color='b', marker='o')
    pl.PlotPolygons2D(l4[:, 0], l4[:, 1], 1, color='m', marker='o')


    #x1 = [4.4388371941530949, 4.4879558708267488]
    #y1 = [-3.7138585019166284, -3.6936360969799154]
    #x2 = [4.3916720661976889, 4.4426284461966903]
    #y2 = [-3.7435296275342731, -3.7230237727210027]
    x1 = l1[[v11, v12, v13], 0]
    y1 = l1[[v11, v12, v13], 1]
    x2 = l2[[v21, v22, v23], 0]
    y2 = l2[[v21, v22, v23], 1]

    #simgrid = dh.ReadTraduitOutB2us(simdir + '/traduit.out.b2us')
    #pl.PlotGridCells(simgrid, 0)

    pl.PlotPolygons2D(x1, y1, 1, color='r', marker='o')
    pl.PlotPolygons2D(x2, y2, 1, color='g', marker='o')
    pl.PlotPolygons2D([x1[0], x2[0]], [y1[0], y2[0]], 1, color='k') # previous face
    pl.PlotPolygons2D([x1[0], x2[1]], [y1[0], y2[1]], 1, color='b') # triangle 1
    pl.PlotPolygons2D([x1[1], x2[0]], [y1[1], y2[0]], 1, color='m') # triangle 2
    pl.PlotPolygons2D([x1[1], x2[1]], [y1[1], y2[1]], 1, color='k', linestyle='--') # quad
    pl.PlotPolygons2D([x1[1], x2[2]], [y1[1], y2[2]], 1, color='b', linestyle='--') # next tria 1
    pl.PlotPolygons2D([x1[2], x2[1]], [y1[2], y2[1]], 1, color='m', linestyle='--') # next tria 2
    pl.PlotPoints2DWithID(x1, y1, [11, 12], 1)
    pl.PlotPoints2DWithID(x2, y2, [21, 22], 1)
    pl.SetAxesLimits2D(plt.gca(), l1[:, 0], l1[:, 1])
    #pl.SetAxesLimits2D(plt.gca(), x1, y1)
    #pl.ShowFigures()
except:
    print("whatever")

# Plot stuff
#-----------
#[R, Z, Psi] = dh.ReadRZPsiFile(simdir + '/rzpsi.dat')
structure = dh.ReadStructureFile(simdir + '/' + structurefile)
#dh.WriteRZPsiFile(writedir, R, Z, Psi)
#pl.PlotStructured2DContourf(R, Z, Psi, 3)
pl.PlotStructure(structure, 5)
pl.PlotStructure(structure, 12)
pl.PlotStructure(structure, 10)
pl.PlotStructure(structure, 16)
pl.PlotStructure(structure, 18)

#pl.PlotPolygonData(outputdir + '/temppol.dat', 5, color = 'r')
#structure = dh.ReadStructureFile(structuredir + '/' + structurefile)
#pl.PlotStructure(structure, 2)

try: 
    simgrid = dh.ReadTraduitOutB2us(simdir + '/' + griddir)
    #tf = 31
    #tfind = tf - 1
    #pl.PlotPoints2D(simgrid.vert.x[[simgrid.face.v1[tfind]-1, simgrid.face.v2[tfind]-1]], 
    #    simgrid.vert.y[[simgrid.face.v1[tfind]-1, simgrid.face.v2[tfind]-1]], 0, 
    #    color='r', marker='o', linewidth=3)
    
    pl.PlotGridCells(simgrid, 0)
    pl.PlotGridFaceLabels(simgrid, 0)
    pl.PlotGridFaceLabels(simgrid, 18)
    #pl.PlotGridVertFieldlineID(simgrid, 0)
    
    #pl.PlotGridCells(simgrid, 1)
    pl.PlotGridCells(simgrid, 3)
    pl.PlotGridCells(simgrid, -1)
    pl.SetAxesLimits2D(plt.gca(), np.array([1.8, 3.0]), np.array([-3.0, -1.8]))
    #pl.PlotGridFaceRegions(simgrid, -1)
    #pl.PlotCellBasedQuantity2D(simgrid, simgrid.cell.ft, 2)
    #pl.PlotCellBasedQuantity2D(simgrid, simgrid.cell.cflags, 3)
    pl.PlotCellBasedQuantity2D(simgrid, simgrid.cell.bt, 14)
    pl.PlotCellBasedQuantity2D(simgrid, simgrid.cell.bp, 13)
    pl.PlotGridTopologicalData(simgrid, 17)
    #pl.PlotCellBasedQuantity2D(simgrid, simgrid.cell.region, 4)
except: 
    print("could not plot simulation grid")

try: 
    topomesh = dh.ReadTopomeshFile(outputdir + '/topomesh.dat')
    pl.PlotTopologicalMesh(topomesh, 5)
    
    pl.PlotTopologicalMesh(topomesh, 0)
    pl.PlotTopologicalMesh(topomesh, 1)
    pl.PlotTopologicalMesh(topomesh, 2)
    pl.PlotTopologicalMesh(topomesh, 3)
    pl.PlotTMCellBasedQuantity(topomesh, np.arange(0, topomesh.cell.ntot), 7)
except:
    print("could not plot final topological mesh")

try: 
    topomesh_base = dh.ReadTopomeshFile(outputdir + '/topomesh_base.dat')
    pl.PlotTopologicalMesh(topomesh_base, 6)
except:
    print("could not plot base topological mesh")
try:
    ggtmdata2 = dh.ReadGGTMDataFile(outputdir + '/ggtmdata_before_pp.dat')
    ggtmdata = dh.ReadGGTMDataFile(outputdir + '/ggtmdata_after_vertexdistribution.dat')
    grid = dh.ReadGGGridDataFile(outputdir + '/grid_after_cellconstruction.dat')

    pl.PlotGGTMDataFaceVertexDistribution(ggtmdata, topomesh, -2)
    pl.PlotGGTMDataCellVertexDistribution(ggtmdata, topomesh, -2)
    pl.PlotGGTMDataFaceVertexDistribution(ggtmdata2, topomesh, -4)
    pl.PlotGGTMDataCellVertexDistribution(ggtmdata2, topomesh, -4)
    pl.PlotGridVertices(grid, 9)
    pl.PlotGridFaces(grid, 9)
    pl.PlotTopologicalMesh(topomesh, 9)
    pl.Plot2DSurfaceDataContourf(outputdir + '/gg_vd_radialdensityfunction.dat', 12, levels=mylevelsdf)
    pl.Plot2DSurfaceDataContourf(outputdir + '/gg_vd_poloidaldensityfunction.dat', 11, levels=mylevelsdf)
    pl.Plot2DSurfaceDataContourf(outputdir + '/Lmaxpolref.dat', 16, levels=mylevelsdf)

except:
    print("could not plot intermediate data")

# Diagnostics
try: 
    # Intersecting faces
    faceIDs, vals = dh.GetVertexCoordinatesWithID(outputdir + '/gg_faceintersections.dat')
    pl.PlotPoints2D(vals[:, 0], vals[:, 1], 9, color='g', marker='o')
except:
    print("could not plot diagnostics")


#pl.PlotGridVertices(grid, 2)
#pl.PlotGridCells(grid, 3)
#pl.PlotPolygonData(outputdir + '/vesselpolygon.dat', 3, color='r', linewidth=0.5)
#pl.PlotPolygonData(outputdir + '/vesselpolygon.dat', 1, color='r', linewidth=0.5)

#pl.PlotGridVertices(grid, 1)
#pl.PlotTopologicalMeshCells(topomesh, 1)
#pl.PlotPolygonData(outputdir + '/extrema_fx_lines.dat', 2, color='r')
#pl.PlotPolygonData(outputdir + '/extrema_fy_lines.dat', 2, color='g')
#pl.PlotPolygonData(outputdir + '/mfcontours.dat', 5, color='k')

#mylevels = [-1, -0.1, -0.01, -0.001, -0.0001, -0.000001, 0, 0.000001, 0.0001, 0.001, 0.01, 0.1, 1]
#mylevelsbias = [0.95, 0.96, 0.97, 0.98, 0.99, 1, 1.01, 1.02, 1.03, 1.04, 1.05]
#pl.Plot2DSurfaceDataContourf(outputdir + '/generalplf.dat', 0, levels=mylevels)
#pl.Plot2DSurfaceDataContourf(outputdir + '/vesselplf.dat', 5, levels=mylevels)
#pl.PlotPolygonData(outputdir + '/vesselpolygon.dat', 5, color='m')
# Show figures
#-------------


# pl.PlotPolygonData(outputdir + '/testpolyg.dat', 2)
# pl.Plot2DSurfaceDataContourf(outputdir + '/costfunctionLR_vesselcontours.dat', 2)
#pl.Plot2DSurfaceDataContourf(outputdir + '/constraints_boundary_plf.dat', 2, levels = mylevels)
#pl.PlotBoundaryConstraintVertices(outputdir, 3)


pl.ShowFigures()
