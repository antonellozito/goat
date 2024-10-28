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
    
    # Read total amount of cells
    i = i + 1
    ntot = np.fromstring(alllines[i], dtype=int, count=3, sep=' ')
    topomesh.cell.ntot = ntot[0]
    topomesh.cell.nvert = ntot[1]
    topomesh.cell.nface = ntot[2]
    
    # Initialize fields
    topomesh.cell.Initialize(topomesh.cell.ntot, topomesh.cell.nvert, topomesh.cell.nface)

    # Skip header
    i = i + 2

    # Read vertex data
    for j in np.arange(0, topomesh.cell.ntot, 1):
        # Split the string
        values = alllines[i+j].split()

        # Add the values
        topomesh.cell.ID[j] = np.fromstring(values[0], dtype=int, count=1, sep=' ')
        topomesh.cell.nc[j] = np.fromstring(values[1], dtype=int, count=1, sep=' ')

    # Update file position
    i = i + topomesh.cell.ntot-1

    # Get header position
    i = 0
    while i < len(alllines): 
        if "cell vertices" in alllines[i]:
            break 
        else: 
            i = i + 1

    # Update counter
    i = i + 1

    # Read vertex data
    for j in np.arange(0, topomesh.cell.nvert, 1):
        # Split the string
        values = alllines[i+j].split()

        # Add the values
        topomesh.cell.vert[j] = np.fromstring(values[0], dtype=int, count=1, sep=' ')

    # Get header position
    i = 0
    while i < len(alllines): 
        if "cell vertex pointer" in alllines[i]:
            break 
        else: 
            i = i + 1

    # Update counter
    i = i + 1

    # Read vertex data
    for j in np.arange(0, topomesh.cell.ntot, 1):
        # Split the string
        values = alllines[i+j].split()

        # Add the values
        tv = np.fromstring(values[0], dtype=int, count=2, sep=' ')
        topomesh.cell.vertP[j, 0:2] = tv

    # Read cell coordinates
    #----------------------
    # Get header position
    i = 0
    while i < len(alllines): 
        if "cell polygons" in alllines[i]:
            break 
        else: 
            i = i + 1

    # Update counter
    i = i + 1

    # Read coordinates
    allfound = False 
    if topomesh.cell.ntot == 0:
        allfound = True
    counter = 0
    while (i < len(alllines)) and (not allfound):
        # Read until we find 'cell'
        if "cell" in alllines[i]:
            # Get cell ID
            values = alllines[i].split()
            fID = np.fromstring(values[1], dtype=int, count=1, sep=' ')
            fID = fID[0] - 1 
            counter = counter + 1
            if (counter == topomesh.cell.ntot):
                allfound = True 

            # Update position
            i = i + 1

            # Read nc coordinates
            topomesh.cell.data[fID].Initialize(topomesh.cell.nc[fID])
            for k in np.arange(0, topomesh.cell.nc[fID], 1):
                values = alllines[i+k].split()
                this = np.fromstring(values[0], dtype=float, count=1)
                topomesh.cell.data[fID].x[k] = np.fromstring(values[0], dtype=float, count=1, sep=' ')
                topomesh.cell.data[fID].y[k] = np.fromstring(values[1], dtype=float, count=1, sep=' ')

            # Update position
            i = i + topomesh.cell.nc[fID]
        else: 
            i = i + 1 

    # Return 
    return topomesh 

def ReadGGTMDataFile(filepath):
    # Description
    #------------
    # Read in the grid generation data of a topological mesh by reading 
    # the ggtmdata file (format specific here)
    ggtmdata = gt.GGTMData()

    # Open file
    thisfile = open(filepath)

    # Read lines
    alllines = thisfile.readlines()

    # Read basic face and cell data
    i = 0
    while i < len(alllines):
        if "faces" in alllines[i]:
            break 
        else: 
            i = i + 1
    i = i + 1
    values = alllines[i].split()
    nftot = np.fromstring(values[0], dtype=int, count=1, sep =' ')
    nf = np.fromstring(values[1], dtype=int, count=1, sep =' ')
    nf = nf[0]

    i = 0
    while i < len(alllines):
        if "cells" in alllines[i]:
            break 
        else: 
            i = i + 1
    i = i + 1
    values = alllines[i].split()
    nctot = np.fromstring(values[0], dtype=int, count=1, sep =' ')
    nc = np.fromstring(values[1], dtype=int, count=1, sep =' ')
    nc = nc[0]

    # Initialize
    ggtmdata.Initialize(nftot[0], nctot[0])

    # Read in face data
    #------------------
    # Get the face header
    i = 0 # reset counter to start from beginning
    while i < len(alllines):
        if "faces" in alllines[i]:
            break 
        else: 
            i = i + 1

    # Read until ID header is found
    while i < len(alllines):
        if "ID, nc" in alllines[i]:
            break 
        else:
            i = i + 1

    # Skip header
    i = i + 1

    # Read face data 
    for j in np.arange(0, nf):
        # Split the string
        values = alllines[i+j].split()

        # Add the values
        ID = np.fromstring(values[0], dtype=int, count=1, sep=' ')
        ID = ID[0]
        nv = np.fromstring(values[1], dtype=int, count=1, sep=' ')
        nv = nv[0]
        ggtmdata.face[ID-1].nv = nv
        ggtmdata.face[ID-1].ID = ID

    # Update position
    i = i + nf - 1

    # Get the face coordinates
    while i < len(alllines): 
        if "face coordinates" in alllines[i]:
            break 
        else: 
            i = i + 1

    # Skip header
    i = i + 1

    # Read 
    counter = 0
    while (i < len(alllines)) and (counter < nf):
        # Read until we find 'face'
        if "face" in alllines[i]:
            # Get face ID
            values = alllines[i].split()
            fID = np.fromstring(values[1], dtype=int, count=1, sep=' ')
            fID = fID[0] - 1 
            counter = counter + 1

            # Update position
            i = i + 1

            # Read nc coordinates
            ggtmdata.face[fID].Initialize(ggtmdata.face[fID].nv, ggtmdata.face[fID].ID)
            for k in np.arange(0, ggtmdata.face[fID].nv, 1):
                values = alllines[i+k].split()
                ggtmdata.face[fID].vert[k] = np.fromstring(values[0], dtype=int, count=1, sep=' ')
                ggtmdata.face[fID].x[k] = np.fromstring(values[1], dtype=float, count=1, sep=' ')
                ggtmdata.face[fID].y[k] = np.fromstring(values[2], dtype=float, count=1, sep=' ')

            # Update position
            i = i + ggtmdata.face[fID].nv
        else: 
            i = i + 1 

    # Read cell data
    #---------------
    # Get the cell header
    i = 0 # reset counter to start from beginning
    while i < len(alllines):
        if "cells" in alllines[i]:
            break 
        else: 
            i = i + 1

    # Read until ID header is found
    while i < len(alllines):
        if "ID, srf, erf, cell line size" in alllines[i]:
            break 
        else:
            i = i + 1

    # Skip header
    i = i + 1

    # Read cell data 
    for j in np.arange(0, nc):
        # Split the string
        values = alllines[i+j].split()

        # Add the values
        ID = np.fromstring(values[0], dtype=int, count=1, sep=' ')
        ID = ID[0]
        ggtmdata.cell[ID-1].ID = ID
        srf = np.fromstring(values[1], dtype=int, count=1, sep=' ')
        erf = np.fromstring(values[2], dtype=int, count=1, sep=' ')
        nl = np.fromstring(values[3], dtype=int, count=1, sep=' ')
        ggtmdata.cell[ID-1].srf = srf[0]
        ggtmdata.cell[ID-1].erf = erf[0]
        ggtmdata.cell[ID-1].nl = nl[0]
        
    # Update file position
    i = i + nc-1

    # Get header position
    while i < len(alllines): 
        if "cell lines" in alllines[i]:
            break 
        else: 
            i = i + 1

    # Update counter
    i = i + 1

    # Read line data
    counter = 0
    while (i < len(alllines)) and (counter < nc):
        # Read until we find 'cell'
        if "cell" in alllines[i]:
            # Get cell ID
            values = alllines[i].split()
            fID = np.fromstring(values[1], dtype=int, count=1, sep=' ')
            fID = fID[0] - 1 
            counter = counter + 1

            # Read nl lines
            ggtmdata.cell[fID].Initialize(ggtmdata.cell[fID].nl-2, ggtmdata.cell[fID].ID)
            for k in np.arange(0, ggtmdata.cell[fID].nl+2, 1):
                # Update position
                i = i + 1

                # Read size of line
                values = alllines[i].split()

                # Initialize line
                nv = np.fromstring(values[0], dtype=int, count=1, sep=' ')
                nv = nv[0]
                line = gt.GGTMLine()
                line.Initialize(nv)

                # Read 
                for cc in np.arange(0, nv, 1): 
                    values = alllines[i+cc+1].split()
                    line.vert[cc] = np.fromstring(values[0], dtype=int, count=1, sep=' ')
                    line.x[cc] = np.fromstring(values[1], dtype=float, count=1, sep=' ')
                    line.y[cc] = np.fromstring(values[2], dtype=float, count=1, sep=' ')

                # Add line
                if k == 0:
                    # First line -> hfline
                    ggtmdata.cell[fID].hfline = line
                elif k == ggtmdata.cell[fID].nl+1:
                    # Last line -> lfline
                    ggtmdata.cell[fID].lfline = line
                else:
                    # intermediate line
                    ggtmdata.cell[fID].lines[k-1] = line

                # Update position
                i = i + nv

        else: 
            i = i + 1 

    # Return
    return ggtmdata

def ReadGGGridDataFile(filepath):
    # Description
    #------------
    # This routine reads in the grid data from an intermediate grid 
    # generator grid object. This is not the same as the full grid
    # that is used in the grid deformation module.

    # Initialize
    grid = gt.GGGrid()

    # Open file
    thisfile = open(filepath)

    # Read lines
    alllines = thisfile.readlines()

    # Read vertex data
    #-----------------
    i = 0
    while i < len(alllines):
        if "vertices" in alllines[i]:
            break 
        else: 
            i = i + 1

    # Skip header
    i = i + 1

    # Read number of vertices
    values = alllines[i].split()
    nv = np.fromstring(values[0], dtype=int, count=1, sep =' ')
    nv = nv[0]

    # Initialize vertex structure
    grid.vert.Initialize(nv)

    # Read vertex data
    while i < len(alllines):
        if "ID, x, y, fieldlineID" in alllines[i]:
            break 
        else: 
            i = i + 1

    # Skip header
    i = i + 1

    # Start reading
    for j in np.arange(0, nv):
        values = alllines[i+j].split()
        ID = np.fromstring(values[0], dtype=int, count=1, sep =' ')
        ID = ID[0]
        x = np.fromstring(values[1], dtype=float, count=1, sep =' ')
        y = np.fromstring(values[2], dtype=float, count=1, sep =' ')
        fID = np.fromstring(values[3], dtype=int, count=1, sep =' ')
        grid.vert.ID[ID-1] = ID
        grid.vert.x[ID-1] = x[0]
        grid.vert.y[ID-1] = y[0]
        grid.vert.fieldlineID[ID-1] = fID[0] 

    # Read in face data
    #------------------
    i = 0
    while i < len(alllines):
        if "faces" in alllines[i]:
            break 
        else: 
            i = i + 1

    # Skip header
    i = i + 1

    # Read number of faces
    values = alllines[i].split()
    nf = np.fromstring(values[0], dtype=int, count=1, sep =' ')
    nf = nf[0]

    # Initialize vertex structure
    grid.face.Initialize(nf)

    # Read vertex data
    while i < len(alllines):
        if "ID, v1, v2, label, region" in alllines[i]:
            break 
        else: 
            i = i + 1

    # Skip header
    i = i + 1

    # Start reading
    for j in np.arange(0, nf):
        values = alllines[i+j].split()
        ID = np.fromstring(values[0], dtype=int, count=1, sep =' ')
        ID = ID[0]
        v1 = np.fromstring(values[1], dtype=int, count=1, sep =' ')
        v2 = np.fromstring(values[2], dtype=int, count=1, sep =' ')
        label = np.fromstring(values[3], dtype=int, count=1, sep =' ')
        region = np.fromstring(values[4], dtype=int, count=1, sep =' ')
        grid.face.ID[ID-1] = ID
        grid.face.v1[ID-1] = v1[0]
        grid.face.v2[ID-1] = v2[0]
        grid.face.label[ID-1] = label[0] 
        grid.face.region[ID-1] = region[0]

    # Read in cell data
    #------------------
    i = 0
    while i < len(alllines):
        if "cells" in alllines[i]:
            break 
        else: 
            i = i + 1

    # Skip header
    i = i + 1

    # Read number of cells and cell vertices
    values = alllines[i].split()
    nc = np.fromstring(values[0], dtype=int, count=1, sep =' ')
    ncv = np.fromstring(values[1], dtype=int, count=1, sep =' ')
    nc = nc[0]
    ncv = ncv[0]

    # Initialize vertex structure
    grid.cell.Initialize(nc, ncv)

    # Read cell data
    while i < len(alllines):
        if "ID, vp1, vp2, region" in alllines[i]:
            break 
        else: 
            i = i + 1

    # Skip header
    i = i + 1

    # Start reading
    for j in np.arange(0, nc):
        values = alllines[i+j].split()
        ID = np.fromstring(values[0], dtype=int, count=1, sep =' ')
        ID = ID[0]
        v1 = np.fromstring(values[1], dtype=int, count=1, sep =' ')
        v2 = np.fromstring(values[2], dtype=int, count=1, sep =' ')
        region = np.fromstring(values[3], dtype=int, count=1, sep =' ')
        grid.cell.ID[ID-1] = ID
        grid.cell.vp1[ID-1] = v1[0]-1 # Need to account for 0-based indexing
        grid.cell.vp2[ID-1] = v2[0]
        grid.cell.region[ID-1] = region[0]

    # Read cell vertices
    while i < len(alllines):
        if "cell vertices" in alllines[i]:
            break 
        else: 
            i = i + 1

    # Skip header
    i = i + 1

    # Read vertices
    for j in np.arange(0, ncv):
        values = alllines[i+j].split()
        tv = np.fromstring(values[0], dtype=int, count=1, sep =' ')
        tv = tv[0]
        grid.cell.vert[j] = tv

    # Return
    return grid

def ReadTraduitOutB2us(filepath):
    # Description
    #------------
    # This routine reads in the grid data in unstructured traduit file
    # into a 'Grid' object. 

    # Initialize
    grid = gt.Grid()

    # Open file
    thisfile = open(filepath)

    # Read lines
    alllines = thisfile.readlines()

    # Read dimensions
    #----------------
    # Start at top of file
    i = 0
    while i < len(alllines):
        if "nCi,nFc,nVx,nCg,nFs,nFt" in alllines[i]:
            break 
        else: 
            i = i + 1

    # Skip header
    i = i + 1

    # Read dimensions
    values = alllines[i].split()
    nc = np.fromstring(values[0], dtype=int, count=1, sep =' '); nc = nc[0]
    nf = np.fromstring(values[1], dtype=int, count=1, sep =' '); nf = nf[0]
    nv = np.fromstring(values[2], dtype=int, count=1, sep =' '); nv = nv[0]
    nfs = np.fromstring(values[4], dtype=int, count=1, sep =' '); nfs = nfs[0]
    nft = np.fromstring(values[5], dtype=int, count=1, sep =' '); nft = nft[0]

    # Read secondary dimensions (continue from previous line)
    while i < len(alllines):
        if "nCmxVx,nCmxFc,nFmxCv,nVmxCv,nVmxFc" in alllines[i]:
            break 
        else: 
            i = i + 1

    # Skip header
    i = i + 1

    # Read dimensions
    values = alllines[i].split()
    ncv = np.fromstring(values[0], dtype=int, count=1, sep =' '); ncv = ncv[0]
    ncf = np.fromstring(values[1], dtype=int, count=1, sep =' '); ncf = ncf[0]

    # Initialize the grid (except flux surfaces etc - later on)
    grid.vert.Initialize(nv)
    grid.face.Initialize(nf)
    grid.cell.Initialize(nc, ncv, ncf)
    grid.ft.Initialize(nft)
    grid.fs.Initialize(nfs)

    # Initialize vertex structure
    grid.vert.Initialize(nv)

    # Read vertex data
    #-----------------
    # Continue from last line
    while i < len(alllines):
        if "*cf: Vx vxX vxY vxPsi vxBx vxBy vxFfbz" in alllines[i]:
            break 
        else: 
            i = i + 1

    # Skip header
    i = i + 1

    # Start reading
    for j in np.arange(0, nv):
        values = alllines[i+j].split()
        ID = np.fromstring(values[0], dtype=int, count=1, sep =' '); ID = ID[0]
        x = np.fromstring(values[1], dtype=float, count=1, sep =' '); x = x[0]
        y = np.fromstring(values[2], dtype=float, count=1, sep =' '); y = y[0]
        psi = np.fromstring(values[3], dtype=float, count=1, sep =' '); psi = psi[0]
        bx = np.fromstring(values[4], dtype=float, count=1, sep =' '); bx = bx[0]
        by = np.fromstring(values[5], dtype=float, count=1, sep =' '); by = by[0]
        ffbz = np.fromstring(values[6], dtype=int, count=1, sep =' '); ffbz = ffbz[0]
        grid.vert.ID[ID-1] = ID
        grid.vert.x[ID-1] = x
        grid.vert.y[ID-1] = y
        grid.vert.bx[ID-1] = bx
        grid.vert.by[ID-1] = by
        grid.vert.psi[ID-1] = psi
        grid.vert.ffbz[ID-1] = ffbz

    # Read cell data
    #---------------
    # Continue from last line
    while i < len(alllines):
        if "*cf: cv cvVxP(:,1) cvVxP(:,2) cvX cvY psi bp bt cflags(:) cvReg cvFt" in alllines[i]:
            break 
        else: 
            i = i + 1

    # Skip header
    i = i + 1

    # Start reading
    for j in np.arange(0, nc):
        values = alllines[i+j].split()
        ID = np.fromstring(values[0], dtype=int, count=1, sep =' '); ID = ID[0]
        vp1 = np.fromstring(values[1], dtype=int, count=1, sep =' '); vp1 = vp1[0]
        vp2 = np.fromstring(values[2], dtype=int, count=1, sep =' '); vp2 = vp2[0]
        x = np.fromstring(values[3], dtype=float, count=1, sep =' '); x = x[0]
        y = np.fromstring(values[4], dtype=float, count=1, sep =' '); y = y[0]
        psi = np.fromstring(values[5], dtype=float, count=1, sep =' '); psi = psi[0]
        bp = np.fromstring(values[6], dtype=float, count=1, sep =' '); bp = bp[0]
        bt = np.fromstring(values[7], dtype=float, count=1, sep =' '); bt = bt[0]
        cflags = np.fromstring(values[8], dtype=int, count=1, sep =' '); cflags = cflags[0]
        cvreg = np.fromstring(values[9], dtype=int, count=1, sep =' '); cvreg = cvreg[0]
        cvft = np.fromstring(values[10], dtype=int, count=1, sep =' '); cvft = cvft[0]
        
        grid.cell.ID[ID-1] = ID
        ID = ID-1
        grid.cell.vp1[ID] = vp1-1 # Need to account for 0-based indexing
        grid.cell.vp2[ID] = vp2
        grid.cell.x[ID] = x
        grid.cell.y[ID] = y
        grid.cell.psi[ID] = psi
        grid.cell.bp[ID] = bp
        grid.cell.bt[ID] = bt
        grid.cell.cflags[ID] = cflags
        grid.cell.region[ID] = cvreg
        grid.cell.ft[ID] = cvft

    # Update i 
    i = i + nc

    # Vertex and face pointers are the same
    grid.cell.fp1 = grid.cell.vp1
    grid.cell.fp2 = grid.cell.vp2

    # Skip next line for reading cell vertices
    i = i + 1

    # Read cell vertices
    k = 0
    while k < ncv:
        values = alllines[i].split()
        for j in np.arange(0, len(values)):
            ID = np.fromstring(values[j], dtype=int, count=1, sep =' '); ID = ID[0]
            grid.cell.vert[k] = ID
            k = k + 1
        i = i + 1

    # Skip header
    i = i + 1

    # Read cell faces
    k = 0
    while k < ncf:
        values = alllines[i].split()
        for j in np.arange(0, len(values)):
            ID = np.fromstring(values[j], dtype=int, count=1, sep =' '); ID = ID[0]
            grid.cell.face[k] = ID
            k = k + 1
        i = i + 1

    # Read face data
    #---------------
    # Continue from last line
    while i < len(alllines):
        if "*cf: fc fcVx(:,1) fcVx(:,2) fcLbl fcReg fcAligned" in alllines[i]:
            break 
        else: 
            i = i + 1

    # Start reading
    for j in np.arange(0, nf):
        values = alllines[i+j].split()
        ID = np.fromstring(values[0], dtype=int, count=1, sep =' '); ID = ID[0]
        vp1 = np.fromstring(values[1], dtype=int, count=1, sep =' '); vp1 = vp1[0]
        vp2 = np.fromstring(values[2], dtype=int, count=1, sep =' '); vp2 = vp2[0]
        label = np.fromstring(values[3], dtype=int, count=1, sep =' '); label = label[0]
        region = np.fromstring(values[4], dtype=int, count=1, sep =' '); region = region[0]
        aligned = np.fromstring(values[5], dtype=int, count=1, sep =' '); aligned = aligned[0]

        grid.face.ID[ID-1] = ID
        ID = ID-1
        grid.face.v1[ID] = vp1
        grid.face.v2[ID] = vp2
        grid.face.label[ID] = label
        grid.face.region[ID] = region
        grid.face.aligned[ID] = aligned

    # Update i 
    i = i + nf

    # Flux tube
    #----------
    # Continue from last line
    while i < len(alllines):
        if "*cf: ft ftCvP(:,1) ftCvP(:,2) ftFcP(:,1) ftFcP(:,2) ftReg" in alllines[i]:
            break 
        else: 
            i = i + 1

    # Skip header
    i = i + 1

    # Start reading
    for j in np.arange(0, nft):
        values = alllines[i+j].split()
        ID = np.fromstring(values[0], dtype=int, count=1, sep =' '); ID = ID[0]
        cp1 = np.fromstring(values[1], dtype=int, count=1, sep =' '); cp1 = cp1[0]
        cp2 = np.fromstring(values[2], dtype=int, count=1, sep =' '); cp2 = cp2[0]
        fp1 = np.fromstring(values[3], dtype=int, count=1, sep =' '); fp1 = fp1[0]
        fp2 = np.fromstring(values[4], dtype=int, count=1, sep =' '); fp2 = fp2[0]
        region = np.fromstring(values[5], dtype=int, count=1, sep =' '); region = region[0]

        grid.ft.ID[ID-1] = ID
        ID = ID-1
        grid.ft.cp1[ID] = cp1-1 # Need to account for 0-based indexing
        grid.ft.cp2[ID] = cp2
        grid.ft.fp1[ID] = fp1-1 # Need to account for 0-based indexing
        grid.ft.fp2[ID] = fp2
        grid.ft.region[ID] = region

    # Update
    i = i + nft 

    # Compute and initialize face and cell data
    nftf = grid.ft.fp1[nft-1]+grid.ft.fp2[nft-1]
    nftc = grid.ft.cp1[nft-1]+grid.ft.cp2[nft-1]
    grid.ft.InitializeCellData(nftc)
    grid.ft.InitializeFaceData(nftf)

    # Skip header and read
    i = i + 1
    k = 0
    while k < nftc:
        values = alllines[i].split()
        for j in np.arange(0, len(values)):
            ID = np.fromstring(values[j], dtype=int, count=1, sep =' '); ID = ID[0]
            grid.ft.cell[k] = ID
            k = k + 1
        i = i + 1
    

    # Update and read
    i = i + 1
    k = 0
    while k < nftf:
        values = alllines[i].split()
        for j in np.arange(0, len(values)):
            ID = np.fromstring(values[j], dtype=int, count=1, sep =' '); ID = ID[0]
            grid.ft.face[k] = ID
            k = k + 1
        i = i + 1
    

    # Flux surfaces
    #--------------
    # Continue from last line
    while i < len(alllines):
        if "*cf: fs fsFcP(:,1) fsFcP(:,2) fsPsi" in alllines[i]:
            break 
        else: 
            i = i + 1

    # Skip header
    i = i + 1

    # Start reading
    for j in np.arange(0, nfs):
        values = alllines[i+j].split()
        ID = np.fromstring(values[0], dtype=int, count=1, sep =' '); ID = ID[0]
        fp1 = np.fromstring(values[1], dtype=int, count=1, sep =' '); fp1 = fp1[0]
        fp2 = np.fromstring(values[2], dtype=int, count=1, sep =' '); fp2 = fp2[0]
        psi = np.fromstring(values[3], dtype=float, count=1, sep =' '); psi = psi[0]

        grid.fs.ID[ID-1] = ID
        ID = ID-1
        grid.fs.fp1[ID] = fp1-1 # Need to account for 0-based indexing
        grid.fs.fp2[ID] = fp2
        grid.fs.psi[ID] = psi

    # Update
    i = i + nfs 

    # Compute and initialize face data
    nfsf = grid.fs.fp1[nfs-1]+grid.fs.fp2[nfs-1]
    grid.fs.InitializeFaceData(nfsf)

    # Skip header and read
    i = i + 1
    k = 0
    while k < nfsf:
        values = alllines[i].split()
        for j in np.arange(0, len(values)):
            ID = np.fromstring(values[j], dtype=int, count=1, sep =' '); ID = ID[0]
            grid.fs.face[k] = ID
            k = k + 1
        i = i + 1

    # Return
    return grid





    

