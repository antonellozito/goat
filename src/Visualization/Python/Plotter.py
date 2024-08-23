# Main plotter class
import matplotlib as mpl
mpl.use('TkAgg') # set the gui backend for the cluster...
from matplotlib import pyplot as plt
import numpy as np
import Datahandler as dh
import time
import goat_types as gt

#==========================================================================#
#                                                                          #
#                              GLOBAL VARIABLES                            #
#                                                                          #
#==========================================================================#

# Special characters
#-------------------
filesep = '/' # file separator

# Enable gui events
#------------------

# Paths & files
#--------------
# file where grid vertices are stored in [ID, x, y] format
gridverticesfile = 'vertices_init.dat' # all grid vertices (initial coordinates)
bndconverticesfile = 'con_bnd_vertices.dat' # vertices constrained by boundary constraints
ffconverticesfile = 'con_ff_vertices.dat' # vertices constrained by flux function constraints
ffconverticesfilestem = 'con_ff_vertices_'
xpconverticesfile = 'con_xp_vertices.dat' # vertices constrained as x-points
elconvertexpairsfile = 'con_el_vertices.dat' # vertex pairs constrained for edge lengths
orthconvertexpairsfile = 'con_orth_vertices.dat' # vertex pairs constrained for orthogonality
vesselpolygonfile = 'vesselpolygon.dat'

# file where grid cells are stored in polygon format ([x, y] with blank
# lines between polygons)
gridcellsfile = 'cells_init.dat'
gridcellsiteratefile = 'cells_iterate.dat'

# Files with optimization history
goathistoryfile = 'goat_optimization_history.dat'
shapeopthistoryfile = 'so_optimization_history.dat'

# Shape optimization paths
fvpfile = 'so_con_fvp_vertices.dat' # fixed vessel points file
fvffile = 'so_con_fvf_vertices.dat' # fixed vessel flux file
origvesselpolygonfile = 'vesselpolygon_orig_so.dat' # original/initial vessel file
currentvesselpolygonfile = 'vesselpolygon_iterate_so.dat' # current/new vessel file

# Grid generation paths
tm_beforecellsfile = 'topomesh_beforecells.dat'
tmfile = 'topomesh.dat'

#==========================================================================#
#                                                                          #
#                          SPECIFIC PLOTTING ROUTINES                      #
#                                                                          #
#==========================================================================#

#--------------------------------------------------------------------------#
#                                   Grid                                   #
#--------------------------------------------------------------------------#


def PlotGridCells(dirpath, fignum):
    # Description
    #------------
    # Plot the grid cells once by reading in the grid cell polygon data as
    # written out in the cells.dat file

    # Get filepath
    filepath = dirpath + filesep + gridcellsfile
    vesselfilepath = dirpath + filesep + vesselpolygonfile

    # Get data
    vals = dh.GetPolygonCoordinates(filepath)
    vesselvals = dh.GetPolygonCoordinates(vesselfilepath)

    # Plot
    PlotPolygons2D(vals[:, 0], vals[:, 1], fignum, color='r', marker='',
        label='Grid cells')
    PlotPolygons2D(vesselvals[:, 0], vesselvals[:, 1], fignum, color='b', marker='',
        label='Vessel')

    # Set axes
    SetAxesLimits2D(plt.gca(), vesselvals[:, 0], vesselvals[:, 1])

    # Set title and other descriptors
    thisaxes = plt.gca()
    thisaxes.set_title('Grid cells')
    thisaxes.set_xlabel('x [m]')
    thisaxes.set_ylabel('y [m]')
    thisaxes.legend(loc='upper right')

def PlotGridCellsIterate(dirpath, fignum):
    # Description
    #------------
    # Plot the grid cells once by reading in the grid cell polygon data as
    # written out in the cells.dat file

    # Get filepath
    filepath = dirpath + filesep + gridcellsiteratefile

    # Get data
    vals = dh.GetPolygonCoordinates(filepath)

    # Plot
    PlotPolygons2D(vals[:, 0], vals[:, 1], fignum, color='r', marker='',
        linewidth=0.1, label='Grid cells')

    # Set axes
    SetAxesLimits2D(plt.gca(), vals[:, 0], vals[:, 1])

    # Set title and other descriptors
    thisaxes = plt.gca()
    thisaxes.set_title('Grid cells')
    thisaxes.set_xlabel('x [m]')
    thisaxes.set_ylabel('y [m]')
    thisaxes.legend(loc='upper right')

def MonitorGrid(datadir, num, pausetime, maxruntime):
    # Description
    #------------
    # This routine makes a plot that is continuously updated.

    # Time interval to replot [s]
    starttime = time.time()

    # Prepare for gui event loop
    plt.ion()

    # Plot the data
    PlotGridCellsIterate(datadir, 1)
    thisfig = plt.gcf()

    # Loop until time has passed
    while (time.time() - starttime <= maxruntime):
        try: 
            # Plot the data
            PlotGridCellsIterate(datadir, 1)

            # Draw
            thisfig.canvas.draw()
            thisfig.canvas.flush_events()

            # Pause
            time.sleep(pausetime)

            # Clear figure
            ClearCurrentAxes()
        except: 
            time.sleep(pausetime)

#--------------------------------------------------------------------------#
#                              Grid Optimization                           #
#--------------------------------------------------------------------------#

def PlotFluxfunctionConstraintVertices(dirpath, fignum):
    # Description
    #------------
    # Plot the vertices that have their flux values constrained. Vertices
    # should be stored in ID, x, y format.

    # Set filepaths
    # Set the filepaths
    vertfilepath = dirpath + filesep + gridverticesfile
    fffilepath = dirpath + filesep + ffconverticesfilestem

    # Get the grid data
    vals = dh.GetVertexCoordinates(vertfilepath)
    PlotPoints2D(vals[:, 0], vals[:, 1], fignum, color='r', marker='o',
        facecolors='none', label='Vertices')

    # Special points
    try: 
        valscon = dh.GetVertexCoordinates(fffilepath+'sp.dat')
        PlotPoints2D(valscon[:, 0], valscon[:, 1], fignum, color='m',
                 marker='*', label='Special vertices')
        valscon = dh.GetVertexCoordinates(fffilepath+'spfs.dat')
        PlotPoints2D(valscon[:, 0], valscon[:, 1], fignum, color='m',
                 marker='+', label='Special flux surface vertices')
    except:
        print('could not print special points of flux function constraints')

    # Tangency points
    try: 
        valscon = dh.GetVertexCoordinates(fffilepath+'tp.dat')
        PlotPoints2D(valscon[:, 0], valscon[:, 1], fignum, color='m',
                 marker='s', label='Tangency points')
    except:
        print('could not print tangency points of flux function constraints')

    # Fixed points
    try: 
        valscon = dh.GetVertexCoordinates(fffilepath+'fp.dat')
        PlotPoints2D(valscon[:, 0], valscon[:, 1], fignum, color='g',
                 marker='*', label='Fixed points')
    except:
        print('could not print fixed points of flux function constraints')

    # Flux surfaces
    try: 
        valscon = dh.GetVertexCoordinates(fffilepath+'fs.dat')
        PlotPoints2D(valscon[:, 0], valscon[:, 1], fignum, color='b',
                 marker='+', label='Flux surfaces')
    except:
        print('could not print flux surface points of flux function constraints')
    

    # Set axes
    SetAxesLimits2D(plt.gca(), vals[:, 0], vals[:, 1])

    # Set title and other descriptors
    thisaxes = plt.gca()
    thisaxes.set_title('Flux function constraint vertices')
    thisaxes.set_xlabel('x [m]')
    thisaxes.set_ylabel('y [m]')
    thisaxes.legend(loc='upper right')

def PlotBoundaryConstraintVertices(dirpath, fignum):
    # Description
    #------------
    # This routine plots the vertices that are constrained by boundary
    # constraints. The data of all vertex coordinates should be stored in
    # a vertices.dat file, together with the boundary vertex data in a
    # con_bnd_vertices.dat file. Both should be present in the folder
    # specified in 'dirpath'.

    # Set the filepaths
    vertfilepath = dirpath + filesep + gridverticesfile
    bndfilepath = dirpath + filesep + bndconverticesfile

    # Get the data
    vals = dh.GetVertexCoordinates(vertfilepath)
    valscon = dh.GetVertexCoordinates(bndfilepath)

    # Plot the data
    PlotPoints2D(vals[:, 0], vals[:, 1], fignum, color='r', marker='o',
        facecolors='none', label='Vertices')
    PlotPoints2D(valscon[:, 0], valscon[:, 1], fignum, color='b',
        marker='+', label='Constrained vertices')

    # Set axes
    SetAxesLimits2D(plt.gca(), vals[:, 0], vals[:, 1])

    # Set title and other descriptors
    thisaxes = plt.gca()
    thisaxes.set_title('Vertices constrained with boundary constraints')
    thisaxes.set_xlabel('x [m]')
    thisaxes.set_ylabel('y [m]')
    thisaxes.legend(loc='upper right')

def PlotXPointConstraintVertices(dirpath, fignum):
    # Description
    #------------
    # Plot the vertices that are x-points and are constrained as such

    # Set filepaths
    # Set the filepaths
    vertfilepath = dirpath + filesep + gridverticesfile
    confilepath = dirpath + filesep + xpconverticesfile

    # Get the data
    vals = dh.GetVertexCoordinates(vertfilepath)
    valscon = dh.GetVertexCoordinates(confilepath)

    # Plot the data
    PlotPoints2D(vals[:, 0], vals[:, 1], fignum, color='r', marker='o',
                 facecolors='none', label='Vertices')
    PlotPoints2D(valscon[:, 0], valscon[:, 1], fignum, color='b',
                 marker='+', label='Constrained vertices')

    # Set axes
    SetAxesLimits2D(plt.gca(), vals[:, 0], vals[:, 1])

    # Set title and other descriptors
    thisaxes = plt.gca()
    thisaxes.set_title('X-point constraint vertices')
    thisaxes.set_xlabel('x [m]')
    thisaxes.set_ylabel('y [m]')
    thisaxes.legend(loc='upper right')

def PlotEdgelengthsConstraintEdges(dirpath, fignum):
    # Description
    #------------
    # Plot the vertex pairs that belong to the edgelengths.

    # Set filepaths
    # Set the filepaths
    cellfilepath = dirpath + filesep + gridcellsfile
    confilepath = dirpath + filesep + elconvertexpairsfile

    # Get the data
    vals = dh.GetPolygonCoordinates(cellfilepath)
    valscon = dh.GetVertexPairCoordinates(confilepath)

    # Plot the data
    PlotPolygons2D(vals[:, 0], vals[:, 1], fignum, color='r',
        label='Grid faces')
    PlotPoints2D(0.5*valscon[:, 0] + 0.5*valscon[:, 2],
        0.5*valscon[:, 1] + 0.5*valscon[:, 3], fignum, color='b',
        marker='+', label='Constrained edges')

    # Set axes
    SetAxesLimits2D(plt.gca(), vals[:, 0], vals[:, 1])

    # Set title and other descriptors
    thisaxes = plt.gca()
    thisaxes.set_title('Edge length constrained edges')
    thisaxes.set_xlabel('x [m]')
    thisaxes.set_ylabel('y [m]')
    thisaxes.legend(loc='upper right')

def PlotOrthogonalityConstraintEdges(dirpath, fignum):
    # Description
    #------------
    # Plot the vertex pairs that belong to the edgelengths.

    # Set filepaths
    # Set the filepaths
    cellfilepath = dirpath + filesep + gridcellsfile
    confilepath = dirpath + filesep + orthconvertexpairsfile

    # Get the data
    vals = dh.GetPolygonCoordinates(cellfilepath)
    valscon = dh.GetVertexPairCoordinates(confilepath)

    # Plot the data
    PlotPolygons2D(vals[:, 0], vals[:, 1], fignum, color='r',
        label='Grid faces')
    PlotPoints2D(0.5*valscon[:, 0] + 0.5*valscon[:, 2],
        0.5*valscon[:, 1] + 0.5*valscon[:, 3], fignum, color='b',
        marker='o', label='Constrained edges')

    # Set axes
    SetAxesLimits2D(plt.gca(), vals[:, 0], vals[:, 1])

    # Set title and other descriptors
    thisaxes = plt.gca()
    thisaxes.set_title('Orthogonality constrained edges')
    thisaxes.set_xlabel('x [m]')
    thisaxes.set_ylabel('y [m]')
    thisaxes.legend(loc='upper right')

def PlotLinefoldingConstraintEdges(dirpath, fignum):
    # Description
    #------------
    # Plot the vertex pairs that belong to the linefolding constraints.
    # This includes poloidal, radial, and vessel line folding 
    # constraints. 

    # General
    #--------
    # Underlying grid
    cellfilepath = dirpath + filesep + gridcellsfile
    vals = dh.GetPolygonCoordinates(cellfilepath)
    PlotPolygons2D(vals[:, 0], vals[:, 1], fignum, color='r',
        label='Grid faces')

    # Poloidal
    #---------
    # Set filepaths
    confilepath = dirpath + filesep + 'con_lf_vpairspol.dat'

    # Get the data
    valscon = dh.GetVertexPairCoordinates(confilepath)

    # Plot the data
    PlotPoints2D(0.5*valscon[:, 0] + 0.5*valscon[:, 2],
        0.5*valscon[:, 1] + 0.5*valscon[:, 3], fignum, color='b',
        marker='o', label='Constrained edges')
    
    # Radial
    #-------
    confilepath = dirpath + filesep + 'con_lf_vpairsrad.dat'

    # Get the data
    valscon = dh.GetVertexPairCoordinates(confilepath)

    # Plot the data
    PlotPoints2D(0.5*valscon[:, 0] + 0.5*valscon[:, 2],
        0.5*valscon[:, 1] + 0.5*valscon[:, 3], fignum, color='g',
        marker='x', label='Constrained edges')
    
    # Vessel
    #-------
    confilepath = dirpath + filesep + 'con_lf_vpairsves.dat'

    # Get the data
    valscon = dh.GetVertexPairCoordinates(confilepath)

    # Plot the data
    PlotPoints2D(0.5*valscon[:, 0] + 0.5*valscon[:, 2],
        0.5*valscon[:, 1] + 0.5*valscon[:, 3], fignum, color='k',
        marker='+', label='Constrained edges')

    # Set axes
    SetAxesLimits2D(plt.gca(), vals[:, 0], vals[:, 1])

    # Set title and other descriptors
    thisaxes = plt.gca()
    thisaxes.set_title('Linefolding constrained edges')
    thisaxes.set_xlabel('x [m]')
    thisaxes.set_ylabel('y [m]')
    thisaxes.legend(loc='upper right')

def PlotGoatOptimizationHistory(dirpath, fignum):
    # Description
    #------------
    # Plot the goat optimization history on a log(value)-iterate plot.
    # This includes the cost function, (max) value of 
    # constraints, the convergence history (infinity norm of 
    # lagrangian gradient), and the linesearch step length 

    # Read data
    #----------
    historypath = dirpath + filesep + goathistoryfile
    [valnames, vals] = dh.ReadGeneralColumnwiseFloatData(historypath)

    # Hedge for negative numbers
    for j in range(len(vals[1,:])):
        for i in range(len(vals[:, j])):
            if (vals[i, j] <= 0):
                vals[i, j] = np.nan

    # Set figure
    #-----------
    plt.figure(fignum)

    # Plot data
    #----------
    # First entry should be iteration counter, then convnorm, dphi, L, 
    # J, max(G), max(H), rxf, step, tol, ...

    # Convnorm
    plt.plot(vals[:, 0], vals[:, 1], 'rx-', label='max(abs(grad L))')
    
    # Cost function
    plt.plot(vals[:, 0], vals[:, 4], 'bx-', label='J')

    # Equality constraints
    plt.plot(vals[:, 0], vals[:, 5], 'gx-', label='max(G)')

    # Inequality constraints
    plt.plot(vals[:, 0], vals[:, 6], 'mx-', label='max(H)')

    # Line search step length
    plt.plot(vals[:, 0], vals[:, 8], 'kx-', label='alpha_ls')

    # Set figure data
    #----------------
    # Set axes
    SetAxesLimitsLogplot(plt.gca(), vals[:, 0], vals[:, [1, 4, 5, 6, 8]])

    # Set title and other descriptors
    thisaxes = plt.gca()
    thisaxes.set_title('Goat grid deformation convergence history')
    thisaxes.set_xlabel('iteration number')
    thisaxes.set_ylabel('value')
    thisaxes.set_yscale('log')
    thisaxes.legend(loc='upper right')

#--------------------------------------------------------------------------#
#                             Grid Generation                              #
#--------------------------------------------------------------------------#

# Topological mesh plotting from file
def PlotTopologicalMeshFromFile(dirpath, fignum):
    filepath = dirpath + filesep + tmfile
    topomesh = dh.ReadTopomeshFile(dirpath)
    PlotTopologicalMesh(topomesh, fignum)

# Topological mesh plotting
def PlotTopologicalMesh(topomesh, fignum):
    # Description
    #------------
    # Plot the topological mesh 

    # Plot the vertices
    #------------------
    # Plot per vertex type
    regular = np.where(topomesh.vert.type == 0) 
    split = np.where(topomesh.vert.type == gt.TMvertexsplitID) 
    minima = np.where(topomesh.vert.type == gt.TMvertexminID)
    saddle = np.where(topomesh.vert.type == gt.TMvertexsaddleID)
    maxima = np.where(topomesh.vert.type == gt.TMvertexmaxID)
    tp1 = np.where(topomesh.vert.type == gt.TMvertextp1ID)
    tp2 = np.where(topomesh.vert.type == gt.TMvertextp2ID)
    PlotPoints2D(topomesh.vert.x[maxima], topomesh.vert.y[maxima], fignum, color='b', 
        marker='o', label='field maxima')
    PlotPoints2D(topomesh.vert.x[saddle], topomesh.vert.y[saddle], fignum, color='b', 
        marker='x', label='field saddle')
    PlotPoints2D(topomesh.vert.x[minima], topomesh.vert.y[minima], fignum, color='b', 
        marker='s', label='field minima')
    PlotPoints2D(topomesh.vert.x[tp1], topomesh.vert.y[tp1], fignum, color='r', 
        marker='o', label='field tangency point type 1')
    PlotPoints2D(topomesh.vert.x[tp2], topomesh.vert.y[tp2], fignum, color='r', 
        marker='x', label='field tangency point type 2')
    PlotPoints2D(topomesh.vert.x[regular], topomesh.vert.y[regular], fignum, color='k', 
        marker='.', label='regular vertex')
    PlotPoints2D(topomesh.vert.x[split], topomesh.vert.y[split], fignum, color='g', 
        marker='d', label='splitted face vertex')
    
    # Store bounds
    xb = [np.min(topomesh.vert.x), np.max(topomesh.vert.x)]
    yb = [np.min(topomesh.vert.y), np.max(topomesh.vert.y)]
    
    # Plot the faces
    #---------------
    # Plot per face type
    pf = np.where(topomesh.face.type == gt.TMfacepolID)
    rf = np.where(topomesh.face.type == gt.TMfaceradID)
    bndf = np.where(topomesh.face.type == gt.TMfacebndID)
    sepf = np.where(topomesh.face.type == gt.TMfacesepID)

    for i in pf[0]:
        PlotPolygons2D(topomesh.face.data[i].x, topomesh.face.data[i].y, 
            fignum, color='r')
    for i in rf[0]:
        PlotPolygons2D(topomesh.face.data[i].x, topomesh.face.data[i].y, 
            fignum, color='g')
    for i in bndf[0]:
        PlotPolygons2D(topomesh.face.data[i].x, topomesh.face.data[i].y, 
            fignum, color='k')
    for i in sepf[0]:
        PlotPolygons2D(topomesh.face.data[i].x, topomesh.face.data[i].y, 
            fignum, color='m')
        
    # Check bounds
    for i in np.arange(0, topomesh.face.ntot, 1):
        xb[0] = np.min([xb[0], np.min(topomesh.face.data[i].x)])
        xb[1] = np.max([xb[1], np.max(topomesh.face.data[i].x)])
        yb[0] = np.min([yb[0], np.min(topomesh.face.data[i].y)])
        yb[1] = np.max([yb[1], np.max(topomesh.face.data[i].y)])

    # Set axes
    SetAxesLimits2D(plt.gca(), xb, yb)

    # Set title and other descriptors
    thisaxes = plt.gca()
    thisaxes.set_title('topomesh')
    thisaxes.set_xlabel('x [m]')
    thisaxes.set_ylabel('y [m]')
    thisaxes.legend(loc='upper right')


#--------------------------------------------------------------------------#
#                             Shape Optimization                           #
#--------------------------------------------------------------------------#

def PlotFixedVesselPointsConstraintVertices(dirpath, fignum):
    # Description
    #------------
    # Plot the vessel vertices that have their position constrained. 

    # Set filepaths
    # Set the filepaths
    vertfilepath = dirpath + filesep + origvesselpolygonfile
    filepath = dirpath + filesep + fvpfile 

    # Get the grid data
    vals = dh.GetPolygonCoordinates(vertfilepath)
    PlotPoints2D(vals[:, 0], vals[:, 1], fignum, color='r', marker='o',
        facecolors='none', label='Vertices')

    # Special points
    try: 
        valscon = dh.GetVertexCoordinates(filepath)
        PlotPoints2D(valscon[:, 0], valscon[:, 1], fignum, color='m',
                 marker='*', label='Fixed vessel points')
    except:
        print('could not print fixed vessel point vertex constraints')

    # Set axes
    SetAxesLimits2D(plt.gca(), vals[:, 0], vals[:, 1])

    # Set title and other descriptors
    thisaxes = plt.gca()
    thisaxes.set_title('Fixed vessel points constraint vertices')
    thisaxes.set_xlabel('x [m]')
    thisaxes.set_ylabel('y [m]')
    thisaxes.legend(loc='upper right')

def PlotFixedVesselFluxConstraintVertices(dirpath, fignum):
    # Description
    #------------
    # Plot the vessel vertices that have their flux value constrained. 

    # Set filepaths
    # Set the filepaths
    vertfilepath = dirpath + filesep + origvesselpolygonfile
    filepath = dirpath + filesep + fvffile

    # Get the grid data
    vals = dh.GetPolygonCoordinates(vertfilepath)
    PlotPoints2D(vals[:, 0], vals[:, 1], fignum, color='r', marker='o',
        facecolors='none', label='Vertices')

    # Special points
    try: 
        valscon = dh.GetVertexCoordinates(filepath)
        PlotPoints2D(valscon[:, 0], valscon[:, 1], fignum, color='m',
                 marker='*', label='Fixed vessel points')
    except:
        print('could not print fixed vessel point vertex constraints')

    # Set axes
    SetAxesLimits2D(plt.gca(), vals[:, 0], vals[:, 1])

    # Set title and other descriptors
    thisaxes = plt.gca()
    thisaxes.set_title('Fixed vessel flux constraint vertices')
    thisaxes.set_xlabel('x [m]')
    thisaxes.set_ylabel('y [m]')
    thisaxes.legend(loc='upper right')

def PlotShapeOptimizationHistory(dirpath, fignum):
        # Description
    #------------
    # Plot the goat optimization history on a log(value)-iterate plot.
    # This includes the cost function, (max) value of 
    # constraints, the convergence history (infinity norm of 
    # lagrangian gradient), and the linesearch step length 

    # Read data
    #----------
    historypath = dirpath + filesep + shapeopthistoryfile
    [valnames, vals] = dh.ReadGeneralColumnwiseFloatData(historypath)

    # Hedge for negative numbers
    for j in range(len(vals[1,:])):
        for i in range(len(vals[:, j])):
            if (vals[i, j] <= 0):
                vals[i, j] = np.nan

    # Set figure
    #-----------
    plt.figure(fignum)

    # Plot data
    #----------
    # First entry should be iteration counter, then convnorm, dphi, L, 
    # J, max(G), max(H), rxf, step, tol, ...

    # Convnorm
    plt.plot(vals[:, 0], vals[:, 1], 'rx-', label='max(abs(grad L))')
    
    # Cost function
    plt.plot(vals[:, 0], vals[:, 4], 'bx-', label='J')

    # Equality constraints
    plt.plot(vals[:, 0], vals[:, 5], 'gx-', label='max(G)')

    # Inequality constraints
    plt.plot(vals[:, 0], vals[:, 6], 'mx-', label='max(H)')

    # Line search step length
    plt.plot(vals[:, 0], vals[:, 8], 'kx-', label='alpha_ls')

    # Set figure data
    #----------------
    # Set axes
    SetAxesLimitsLogplot(plt.gca(), vals[:, 0], vals[:, [1, 4, 5, 6, 8]])

    # Set title and other descriptors
    thisaxes = plt.gca()
    thisaxes.set_title('Shape optimization convergence history')
    thisaxes.set_xlabel('iteration number')
    thisaxes.set_ylabel('value')
    thisaxes.set_yscale('log')
    thisaxes.legend(loc='upper right')

def PlotVesselDisplacement(dirpath, fignum):
    # Description
    #------------
    # Plot the vessel displacement between original and new vessel

    # Get filepath
    filepathorig = dirpath + filesep + origvesselpolygonfile 
    filepathnew = dirpath + filesep + vesselpolygonfile 

    # Get data
    valsorig = dh.GetPolygonCoordinates(filepathorig)
    valsnew = dh.GetPolygonCoordinates(filepathnew)

    # Plot both polygons
    PlotPolygons2D(valsorig[:, 0], valsorig[:, 1], fignum, color='b', marker='x',
        label='Original vessel polygon')
    PlotPolygons2D(valsnew[:, 0], valsnew[:, 1], fignum, color='r', marker='o',
        label='New vessel polygon')

    # Plot displacement vector
    d = valsnew - valsorig 
    plt.quiver(valsorig[:, 0], valsorig[:, 1], d[:, 0], d[:, 1], color='g', angles='xy', scale_units='xy', scale=1)

    # Set axes
    SetAxesLimits2D(plt.gca(), valsnew[:, 0], valsnew[:, 1])
    
#--------------------------------------------------------------------------#
#                                 Vessel                                   #
#--------------------------------------------------------------------------#

def PlotInitialVesselPolygon(dirpath, fignum):
    # Description
    # ------------
    # Plot the initial vessel polygon

    # Get filepath
    filepath = dirpath + filesep + origvesselpolygonfile 

    # Get data
    vals = dh.GetPolygonCoordinates(filepath)

    # Plot
    PlotPolygons2D(vals[:, 0], vals[:, 1], fignum, color='b', marker='x',
                   label='Vessel polygon')

    # Set axes
    SetAxesLimits2D(plt.gca(), vals[:, 0], vals[:, 1])

    # Set title and other descriptors
    thisaxes = plt.gca()
    thisaxes.set_title('Vessel polygon')
    thisaxes.set_xlabel('x [m]')
    thisaxes.set_ylabel('y [m]')
    thisaxes.legend(loc='upper right')

def PlotVesselPolygon(dirpath, fignum):
    # Description
    # ------------
    # Plot the grid cells once by reading in the grid cell polygon data as
    # written out in the cells.dat file

    # Get filepath
    filepath = dirpath + filesep + vesselpolygonfile

    # Get data
    vals = dh.GetPolygonCoordinates(filepath)

    # Plot
    PlotPolygons2D(vals[:, 0], vals[:, 1], fignum, color='r', marker='o',
                   label='Vessel polygon')

    # Set axes
    SetAxesLimits2D(plt.gca(), vals[:, 0], vals[:, 1])

    # Set title and other descriptors
    thisaxes = plt.gca()
    thisaxes.set_title('Vessel polygon')
    thisaxes.set_xlabel('x [m]')
    thisaxes.set_ylabel('y [m]')
    thisaxes.legend(loc='upper right')

def PlotVesselIterate(dirpath, fignum):
    # Description
    #------------
    # Plot the vessel polygon obtained during a shape iteration
    # Get filepath
    filepath = dirpath + filesep + currentvesselpolygonfile

    # Get data
    vals = dh.GetPolygonCoordinates(filepath)

    # Plot
    PlotPolygons2D(vals[:, 0], vals[:, 1], fignum, color='r', marker='o',
                   label='Vessel polygon')

    # Set axes
    SetAxesLimits2D(plt.gca(), vals[:, 0], vals[:, 1])

    # Set title and other descriptors
    thisaxes = plt.gca()
    thisaxes.set_title('Vessel polygon')
    thisaxes.set_xlabel('x [m]')
    thisaxes.set_ylabel('y [m]')
    thisaxes.legend(loc='upper right')

def MonitorGridAndVessel(datadir, num, pausetime, maxruntime):
    # Description
    #------------
    # This routine makes a plot that is continuously updated of the
    # current grid and vessel

    # Time interval to replot [s]
    starttime = time.time()

    # Prepare for gui event loop
    plt.ion()

    # Plot the data
    PlotGridCellsIterate(datadir, 1)
    PlotVesselIterate(datadir, 1)
    thisfig = plt.gcf()

    # Loop until time has passed
    while (time.time() - starttime <= maxruntime):
        try: 
            # Plot the data
            PlotGridCellsIterate(datadir, 1)
            PlotVesselIterate(datadir, 1)

            # Draw
            thisfig.canvas.draw()
            thisfig.canvas.flush_events()

            # Pause
            time.sleep(pausetime)

            # Clear figure
            ClearCurrentAxes()
        except: 
            time.sleep(pausetime)

#--------------------------------------------------------------------------#
#                               2D surface plots                           #
#--------------------------------------------------------------------------#
    
# Wrappers 
def Plot2DSurfaceData(filepath, fignum):
    # Description
    #------------
    # Read and plot 2D data as patchplot

    # Read
    data = dh.GetGeneral2DSurfacedata(filepath)

    # Plot
    PlotGeneral2DSurface(data[:, 0], data[:, 1], data[:, 2], fignum)

def Plot2DSurfaceDataContourf(filepath, fignum, **plotargs):
    # Description
    #------------
    # Read and plot 2D data with filled contours

    # Read
    data = dh.GetGeneral2DSurfaceData(filepath)

    # Plot
    PlotGeneral2DContourf(data[:, 0], data[:, 1], data[:, 2], fignum, **plotargs)

    # Add colorbar

def Plot2DSurfaceDataContour(filepath, fignum, **plotargs):
    # Description
    #------------
    # Read and plot 2D data with filled contours

    # Read
    data = dh.GetGeneral2DSurfaceData(filepath)

    # Plot
    PlotGeneral2DContour(data[:, 0], data[:, 1], data[:, 2], fignum, **plotargs)

    # Add colorbar



    

#==========================================================================#
#                                                                          #
#                           GENERAL PLOTTING ROUTINES                      #
#                                                                          #
#==========================================================================#
def PlotPolygons2D(x, y, fignum, **plotargs):
    # General polygon plotter. x and y should be np.arrays containing the
    # coordinates to plot. Fignum should contain the figure number on which
    # to plot the data.

    # Set the current figure
    plt.figure(fignum)

    # Plot the data as a polygon plot
    #fig, ax = plt.subplots()
    plt.plot(x, y, **plotargs)
    
    #myrange = range(len(x))
    #for i, txt in enumerate(myrange):
    #    ax.annotate(txt, (x[i], y[i]))
    #plt.draw()

def PlotPolygonData(filepath, fignum, **plotargs):
    # Read data
    data = dh.GetPolygonCoordinates(filepath)

    # Plot
    plt.figure(fignum)
    PlotPolygons2D(data[:, 0], data[:, 1], fignum, **plotargs)

def PlotVertexPairsData2D(filepath, fignum, **plotargs):
    # Read data
    data = dh.GetVertexPairCoordinates(filepath)

    # Compute faces
    xf = 0.5*(data[:, 0] + data[:, 2])
    yf = 0.5*(data[:, 1] + data[:, 3])

    # Plot
    PlotPoints2D(xf, yf, fignum, **plotargs)

def PlotPoints2D(x, y, fignum, **plotargs):
    # General point plotter.

    # Set the current figure
    plt.figure(fignum)
    plt.scatter(x, y, **plotargs)
    plt.draw()
    
def PlotGeneral2DSurface(x, y, z, fignum, **plotargs):
    # General z = f(x, y) surface plotter - may be expensive since
    # a triangulation is created under the hood. 

    # Set the current figure
    plt.figure(fignum)
    plt.tripcolor(x, y, z, **plotargs)
    plt.draw()

def PlotGeneral2DContourf(x, y, z, fignum, **plotargs):
    # Same as PlotGeneral2DSurface, but now we plot filled contours
    
    # Set the current figure
    fig = plt.figure(fignum)
    CS = plt.tricontourf(x, y, z, **plotargs)
    cbar = fig.colorbar(CS)
    plt.draw()

def PlotGeneral2DContour(x, y, z, fignum, **plotargs):
    # Same as PlotGeneral2DSurface, but now we plot filled contours
    
    # Set the current figure
    fig = plt.figure(fignum)
    CS = plt.tricontour(x, y, z, **plotargs)
    cbar = fig.colorbar(CS)
    plt.draw()

#==========================================================================#
#                                                                          #
#                                 AUXILIARY                                #
#                                                                          #
#==========================================================================#
def ShowFigures(*args, **kwargs):
    # Just a wrapper for plt.show()
    plt.show(*args, **kwargs)

def ClearCurrentAxes():
    # Just a wrapper for plt.clf()
    plt.clf()

def SetAxesLimits2D(thisaxes, xdata, ydata):
    # Automatically set the axes limits based on the x and y figure data.
    # Also applies true scaling

    # Compute data limits - ignore NaNs
    maxx = np.nanmax(xdata)
    minx = np.nanmin(xdata)
    maxy = np.nanmax(ydata)
    miny = np.nanmin(ydata)

    # Compute ranges and add some buffer to the sides
    xrange = maxx - minx
    yrange = maxy - miny
    xbuffer = 0.05*xrange
    ybuffer = 0.05*yrange

    # Compute limits
    xlim = [minx - xbuffer, maxx + xbuffer]
    ylim = [miny - ybuffer, maxy + ybuffer]

    # Set limits
    thisaxes.set_xlim(xlim)
    thisaxes.set_ylim(ylim)

    # Set proper scaling
    thisaxes.set_aspect(1)

def SetAxesLimitsLogplot(thisaxes, xdata, ydata):
    # Automatically set the axes limits based on the x and y figure data.
    # Also applies true scaling

    # Compute data limits - ignore NaNs
    maxx = np.nanmax(xdata)
    minx = np.nanmin(xdata)
    maxy = np.nanmax(ydata)
    miny = np.nanmin(ydata)


    # Compute limits
    xlim = [minx, maxx]
    ylim = [miny, maxy]

    # Set limits
    thisaxes.set_xlim(xlim)
    thisaxes.set_ylim(ylim)

    # Set proper scaling
    thisaxes.set_aspect(1)




