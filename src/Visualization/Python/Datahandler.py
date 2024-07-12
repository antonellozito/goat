# This module handles plotting data. It reads in data and returns the data
# in proper structures.
import numpy as np
import os 

def GetDataDirectory():
    # Description
    #------------
    # This routine determines based on the environment variables which
    # data directory has to be considered. For goat not coupled to 
    # solps, this is simply the current directory '.', when coupled to
    # solps, all output should be written to './output'.

    # Check
    if 'SOLPSTOP' in os.environ:
        datadir = '/output'
    else:
        datadir = '.'

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



