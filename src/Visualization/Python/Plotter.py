# Main plotter class
import matplotlib as mpl
mpl.use('TkAgg') # set the gui backend for the cluster...
from matplotlib import pyplot as plt
import numpy as np
import Datahandler as dh

#==========================================================================#
#                                                                          #
#                              GLOBAL VARIABLES                            #
#                                                                          #
#==========================================================================#

# Special characters
#-------------------
filesep = '/' # file separator

# Paths & files
#--------------
# file where grid vertices are stored in [ID, x, y] format
gridverticesfile = 'vertices_init.dat' # all grid vertices (initial coordinates)
bndconverticesfile = 'con_bnd_vertices.dat' # vertices constrained by boundary constraints
ffconverticesfile = 'con_ff_vertices.dat' # vertices constrained by flux function constraints
xpconverticesfile = 'con_xp_vertices.dat' # vertices constrained as x-points
elconvertexpairsfile = 'con_el_vertices.dat' # vertex pairs constrained for edge lengths
orthconvertexpairsfile = 'con_orth_vertices.dat' # vertex pairs constrained for orthogonality

# file where grid cells are stored in polygon format ([x, y] with blank
# lines between polygons)
gridcellsfile = 'cells.dat'


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

    # Get data
    vals = dh.GetPolygonCoordinates(filepath)

    # Plot
    PlotPolygons2D(vals[:, 0], vals[:, 1], fignum, color='r', marker='',
        label='Grid cells')

    # Set axes
    SetAxesLimits2D(plt.gca(), vals[:, 0], vals[:, 1])

    # Set title and other descriptors
    thisaxes = plt.gca()
    thisaxes.set_title('Grid cells')
    thisaxes.set_xlabel('x [m]')
    thisaxes.set_ylabel('y [m]')
    thisaxes.legend(loc='upper right')


#--------------------------------------------------------------------------#
#                                Optimization                              #
#--------------------------------------------------------------------------#

def PlotFluxfunctionConstraintVertices(dirpath, fignum):
    # Description
    #------------
    # Plot the vertices that have their flux values constrained. Vertices
    # should be stored in ID, x, y format.

    # Set filepaths
    # Set the filepaths
    vertfilepath = dirpath + filesep + gridverticesfile
    fffilepath = dirpath + filesep + ffconverticesfile

    # Get the data
    vals = dh.GetVertexCoordinates(vertfilepath)
    valscon = dh.GetVertexCoordinates(fffilepath)

    # Plot the data
    PlotPoints2D(vals[:, 0], vals[:, 1], fignum, color='r', marker='o',
                 facecolors='none', label='Vertices')
    PlotPoints2D(valscon[:, 0], valscon[:, 1], fignum, color='b',
                 marker='+', label='Constrained vertices')

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
    plt.plot(x, y, **plotargs)
    plt.draw()


def PlotPoints2D(x, y, fignum, **plotargs):
    # General point plotter.

    # Set the current figure
    plt.figure(fignum)
    plt.scatter(x, y, **plotargs)
    plt.draw()

#==========================================================================#
#                                                                          #
#                                 AUXILIARY                                #
#                                                                          #
#==========================================================================#


def ShowFigures():
    # Just a wrapper for plt.show()
    plt.show()

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



