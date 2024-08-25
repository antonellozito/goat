# This module handles plotting data. It reads in data and returns the data
# in proper structures.
import numpy as np
import os 
import goat_types as gt

def GetDataDirectory():
    # Description
    #------------
    # This routine determines based on the environment variables which
    # data directory has to be considered. For goat not coupled to 
    # solps, this is simply the current directory '.', when coupled to
    # solps, all output should be written to './output'.

    # Check
    if 'SOLPSTOP' in os.environ:
        datadir = './output'
    else:
        datadir = './output'

    return datadir

def GetVertexCoordinates(filepath):
    # Description
    #------------
    # This routine reads a formatted file where the vertex IDs and
    # coordinates (2D) are stored as [ID, x, y] columns. White/empty lines
    # are ignored.

    # Read in vertex coordinates from the vertices.dat file
    thisfile = open(filepath)

    alllines = thisfile.readlines()

    # Remove header
    del alllines[0]
    vals = np.zeros([len(alllines), 3])

    # Read in vertex data
    cc = 0

    for i in alllines:
        if i == '\n':  # empty string
            pass
            # Don't read in
        else:
            # Read
            vals[cc, 0:3] = np.fromstring(i, dtype=float, count=3, sep=' ')
            cc = cc + 1

    # Return values
    return vals[0:cc, 1:3]


def GetPolygonCoordinates(filepath):
    # Description
    #------------
    # This routine reads a formatted file where the polygon
    # coordinates (2D) are stored as [x, y] columns. White/empty lines
    # indicate the end of a polygon and are padded with NaNs in the
    # data.

    # Read in vertex coordinates from the vertices.dat file
    thisfile = open(filepath)

    alllines = thisfile.readlines()

    # Remove header
    del alllines[0]
    vals = np.zeros([len(alllines), 2])

    # Read in vertex data
    cc = 0

    for i in alllines:
        if i == '\n':  # empty string
            vals[cc, :] = np.NaN
        else:
            # Read
            vals[cc, 0:2] = np.fromstring(i, dtype=float, count=2, sep=' ')

        # Update counter
        cc = cc + 1

    # Return values
    return vals[0:cc, 0:2]


def GetVertexPairCoordinates(filepath):
    # Description
    #------------
    # Read in a vertex pair format file where the vertices are stored in
    # ID1, ID2, x1, x2, y1, y2 format. The values that are returned are
    # in [x1, y1, x2, y2] format. Empty lines are skipped

    # Read in vertex coordinates from the vertices.dat file
    thisfile = open(filepath)

    alllines = thisfile.readlines()

    # Retrieve sizes
    dim = np.fromstring(alllines[1], dtype=float, sep=' ')
    del alllines[0]
    npoints = np.floor(len(dim)/3).astype(int)
    vals = np.zeros([len(alllines), len(dim)])

    # Read in vertex data
    cc = 0

    for i in alllines:
        if i == '\n':  # empty string
            # don't read
            pass
        else:
            # Read
            vals[cc, 0:6] = np.fromstring(i, dtype=float, count=6, sep=' ')

        # Update counter
        cc = cc + 1

    # Return values
    returnvec = np.zeros([2*npoints], dtype=int)
    for i in np.arange(0, npoints, 1):
        returnvec[2*i] = npoints + i 
        returnvec[2*i+1] = 2*npoints + i 
    return vals[0:cc, returnvec]

def GetGeneral2DSurfaceData(filepath):
    # Description
    #------------
    # Read in x, y, z data where x and y are assumed to be spatial 
    # coordinates and z is some field. Blank lines are allowed but 
    # skipped. It is assumed that the first line is a header. 

    # Read in vertex coordinates from the file
    thisfile = open(filepath)

    alllines = thisfile.readlines()

    # Remove header
    del alllines[0]
    vals = np.zeros([len(alllines), 3])

    # Read in vertex data
    cc = 0

    for i in alllines:
        if i == '\n':  # empty string
            pass
            # Don't read in
        else:
            # Read
            vals[cc, 0:3] = np.fromstring(i, dtype=float, count=3, sep=' ')
            cc = cc + 1

    # Return values
    return vals[0:cc, 0:3]

def ReadGeneralColumnwiseFloatData(filepath):
    # Description
    #------------
    # Read in any columnwise organized float data. It is assumed that 
    # the first line is a header with the same amount of columns as the
    # data later on, and that the string of each header describes 
    # the data. 

    # Read data
    #----------
    # Open file
    thisfile = open(filepath)

    # Read lines
    alllines = thisfile.readlines()

    # Split the header
    valnames = alllines[0].split()

    # Delete the header
    del alllines[0]

    # Count the number of columns
    ncol = len(valnames)

    # Initialize
    vals = np.zeros([len(alllines), ncol])

    # Read data
    cc = 0
    for i in alllines:
        if i == '\n':  # empty string
            pass
            # Don't read in
        else:
            # Read
            vals[cc, 0:ncol] = np.fromstring(i, dtype=float, count=ncol, sep=' ')
            cc = cc + 1
        
    # Return 
    return valnames, vals

def ReadTopomeshFile(filepath):
    # Description
    #------------
    # Read in the data of a topological mesh by reading the topological
    # mesh file (format specific here)
    topomesh = gt.Topomesh()

    # Open file
    thisfile = open(filepath)

    # Read lines
    alllines = thisfile.readlines()

    # Read in vertex data
    #--------------------
    # Get the 'vertices' header
    i = 0
    while i < len(alllines): 
        if "vertices" in alllines[i]:
            break 
        else: 
            i = i + 1

    # Read total amount of vertices
    i = i + 1
    ntot = np.fromstring(alllines[i], dtype=int, count=1, sep=' ')
    topomesh.vert.ntot = ntot[0]
    
    # Initialize fields
    topomesh.vert.Initialize(topomesh.vert.ntot)

    # Skip header
    i = i + 2

    # Read vertex data
    for j in np.arange(0, topomesh.vert.ntot, 1):
        # Split the string
        values = alllines[i+j].split()

        # Add the values
        topomesh.vert.ID[j] = np.fromstring(values[0], dtype=int, count=1, sep=' ')
        topomesh.vert.x[j] = np.fromstring(values[1], dtype=float, count=1, sep=' ')
        topomesh.vert.y[j] = np.fromstring(values[2], dtype=float, count=1, sep=' ')
        topomesh.vert.type[j] = np.fromstring(values[3], dtype=int, count=1, sep=' ')
        topomesh.vert.fval[j] = np.fromstring(values[4], dtype=float, count=1, sep=' ')
        topomesh.vert.BV[j] = np.fromstring(values[5], dtype=int, count=1, sep=' ')

    # Read face data
    #---------------
    # Get the 'faces' header
    i = 0
    while i < len(alllines): 
        if "faces" in alllines[i]:
            break 
        else: 
            i = i + 1

    # Read total amount of faces
    i = i + 1
    ntot = np.fromstring(alllines[i], dtype=int, count=1, sep=' ')
    topomesh.face.ntot = ntot[0]
    
    # Initialize fields
    topomesh.face.Initialize(topomesh.face.ntot)

    # Skip header
    i = i + 2

    # Read vertex data
    for j in np.arange(0, topomesh.face.ntot, 1):
        # Split the string
        values = alllines[i+j].split()

        # Add the values
        topomesh.face.ID[j] = np.fromstring(values[0], dtype=int, count=1, sep=' ')
        topomesh.face.fsID[j] = np.fromstring(values[1], dtype=int, count=1, sep=' ')
        topomesh.face.type[j] = np.fromstring(values[2], dtype=int, count=1, sep=' ')
        topomesh.face.vert[j, 0] = np.fromstring(values[3], dtype=int, count=1, sep=' ')
        topomesh.face.vert[j, 1] = np.fromstring(values[4], dtype=int, count=1, sep=' ')
        topomesh.face.BF[j] = np.fromstring(values[5], dtype=int, count=1, sep=' ')
        topomesh.face.nc[j] = np.fromstring(values[6], dtype=int, count=1, sep=' ')

    # Update file position
    i = i + topomesh.face.ntot-1

    # Read face coordinates
    #----------------------
    # Get header position
    i = 0
    while i < len(alllines): 
        if "face coordinates" in alllines[i]:
            break 
        else: 
            i = i + 1

    # Update counter
    i = i + 1

    # Read coordinates
    allfound = False 
    counter = 0
    while (i < len(alllines)) and (not allfound):
        # Read until we find 'face'
        if "face" in alllines[i]:
            # Get face ID
            values = alllines[i].split()
            fID = np.fromstring(values[1], dtype=int, count=1, sep=' ')
            fID = fID[0] - 1 
            counter = counter + 1
            if (counter == topomesh.face.ntot):
                allfound = True 

            # Update position
            i = i + 1

            # Read nc coordinates
            topomesh.face.data[fID].Initialize(topomesh.face.nc[fID])
            for k in np.arange(0, topomesh.face.nc[fID], 1):
                values = alllines[i+k].split()
                this = np.fromstring(values[0], dtype=float, count=1)
                topomesh.face.data[fID].x[k] = np.fromstring(values[0], dtype=float, count=1, sep=' ')
                topomesh.face.data[fID].y[k] = np.fromstring(values[1], dtype=float, count=1, sep=' ')

            # Update position
            i = i + topomesh.face.nc[fID]
        else: 
            i = i + 1 

    # Read cell data
    #---------------
    # Get header position
    i = 0
    while i < len(alllines): 
        if "cells" in alllines[i]:
            break 
        else: 
            i = i + 1





    # Return 
    return topomesh 

    

