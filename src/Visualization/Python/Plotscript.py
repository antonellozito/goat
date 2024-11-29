# Description
#------------
# Simple plot script that uses the other python functionality
import Plotter as pl 
import Datahandler as dh 

# Set paths
#----------
# Simulation directory
simdir = dh.GetDataDirectory()
simdir = './goatf/Examples/DEMO_WIM'
outputdir = simdir + '/output'
structuredir = './goatf/src/Visualization'
writedir = './goatf/src/Visualization'

# Plot stuff
#-----------
#[R, Z, Psi] = dh.ReadRZPsiFile(simdir + '/rzpsi.dat')
structure = dh.ReadStructureFile(simdir + '/structure.dat')
#dh.WriteRZPsiFile(writedir, R, Z, Psi)
#pl.PlotStructured2DContourf(R, Z, Psi, 3)
pl.PlotStructure(structure, 3)
#pl.PlotPolygonData(outputdir + '/temppol.dat', 5, color = 'r')
#structure = dh.ReadStructureFile(structuredir + '/structure.dat')
#pl.PlotStructure(structure, 2)

#simgrid = dh.ReadTraduitOutB2us(simdir + '/traduit.out.b2us')
#pl.PlotGridCells(simgrid, 0)
#pl.PlotGridFaceLabels(simgrid, 0)
#pl.PlotGridVertFieldlineID(simgrid, 0)
#pl.PlotGridCells(simgrid, 1)
#pl.PlotGridFaceRegions(simgrid, -1)
#pl.PlotCellBasedQuantity2D(simgrid, simgrid.cell.ft, 2)
#pl.PlotCellBasedQuantity2D(simgrid, simgrid.cell.cflags, 3)
#pl.PlotCellBasedQuantity2D(simgrid, simgrid.cell.region, 4)

topomesh = dh.ReadTopomeshFile(outputdir + '/topomesh.dat')
#ggtmdata = dh.ReadGGTMDataFile(outputdir + '/ggtmdata_after_vertexdistribution.dat')
#grid = dh.ReadGGGridDataFile(outputdir + '/grid_after_cellconstruction.dat')

#pl.PlotGGTMDataFaceVertexDistribution(ggtmdata, topomesh, 2)
#pl.PlotGGTMDataCellVertexDistribution(ggtmdata, topomesh, 2)
#pl.PlotGridVertices(grid, 1)
#pl.PlotGridFaces(grid, 1)

#pl.PlotGridVertices(grid, 2)
#pl.PlotGridCells(grid, 3)
#pl.PlotPolygonData(outputdir + '/vesselpolygon.dat', 3, color='r', linewidth=0.5)
#pl.PlotPolygonData(outputdir + '/vesselpolygon.dat', 1, color='r', linewidth=0.5)
pl.PlotTopologicalMesh(topomesh, 5)
#pl.PlotGridVertices(grid, 1)
#pl.PlotTopologicalMeshCells(topomesh, 1)
#pl.PlotPolygonData(outputdir + '/extrema_fx_lines.dat', 2, color='r')
#pl.PlotPolygonData(outputdir + '/extrema_fy_lines.dat', 2, color='g')
#pl.PlotPolygonData(outputdir + '/mfcontours.dat', 5, color='k')
# pl.PlotPolygonData(outputdir + '/topocontours_all_beforeinsertion.dat', 5, color='g')
mylevels = [-0.01, -0.001, -0.0001, -0.000001, 0, 0.000001, 0.0001, 0.001, 0.01]
#mylevelsbias = [0.95, 0.96, 0.97, 0.98, 0.99, 1, 1.01, 1.02, 1.03, 1.04, 1.05]
#pl.Plot2DSurfaceDataContourf(outputdir + '/generalplf.dat', 0, levels=mylevels)
#pl.Plot2DSurfaceDataContourf(outputdir + '/vesselplf.dat', 6, levels=mylevels)
pl.PlotPolygonData(outputdir + '/vesselpolygon.dat', 5, color='m')
# Show figures
#-------------


#pl.Plot2DSurfaceDataContourf(outputdir + '/closedapproximationplf.dat', 2, levels=mylevels)
# pl.PlotPolygonData(outputdir + '/testpolyg.dat', 2)
# pl.Plot2DSurfaceDataContourf(outputdir + '/costfunctionLR_vesselcontours.dat', 2)
#pl.Plot2DSurfaceDataContourf(outputdir + '/constraints_boundary_plf.dat', 2, levels = mylevels)
#pl.PlotBoundaryConstraintVertices(outputdir, 3)


pl.ShowFigures()
