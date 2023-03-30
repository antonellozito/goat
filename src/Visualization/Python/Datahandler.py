# This module handles plotting data. It reads in data and returns the data
# in proper structures.
import numpy as np



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
    # ID1, ID2, x1, y1, x2, y2 format. The values that are returned are
    # in [x1, y1, x2, y2] format. Empty lines are skipped

    # Read in vertex coordinates from the vertices.dat file
    thisfile = open(filepath)

    alllines = thisfile.readlines()

    # Remove header
    del alllines[0]
    vals = np.zeros([len(alllines), 6])

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
    return vals[0:cc, 2:6]


