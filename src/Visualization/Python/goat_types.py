# This module defines some convenient types for post-processing that
# mimick the types of goat 
import numpy as np 

#======================================================================#
#                                                                      #
#                       GRID GENERATOR TYPES                           #
#                                                                      #
#======================================================================#

#----------------------------------------------------------------------#
#                           DEFINITIONS                                #
#----------------------------------------------------------------------#
# Vessel parts definitions
targetID = 1
coreID = 2
outerboundaryID = 3
vesselID = 4
interiorID = 5

# Definitions for grid generation
# Topological mesh vertex IDs (field minimum, saddle point, maximum: 1, 2, 3,
# field & boundary tangency point types 1 and 2: 4, 5, boundary vertex
# that is no tangency point or whatever: 6, vertex constructed by 
# splitting face: -1)
TMvertexminID = 1
TMvertexsaddleID = 2
TMvertexmaxID = 3
TMvertextp1ID = 4
TMvertextp2ID = 5
TMvertexbndID = 6
TMvertexsplitID = -1
TMvertexregularID = 0

# Topological mesh boundary IDs (1: radial face, 2: poloidal face, &
# 3: boundary face)
TMfaceradID = 1
TMfacepolID = 2
TMfacebndID = 3
TMfacesepID = 4

#----------------------------------------------------------------------#
#                        TOPOLOGICAL MESH                              #
#----------------------------------------------------------------------#

# Topological mesh vertices
class TopomeshVert:
    # Definition
    def __init__(self):
        # total number of vertices
        self.ntot = 0

        # Coordinates
        self.x = np.zeros(0, dtype=float)
        self.y = np.zeros(0, dtype=float)

        # ID
        self.ID = np.zeros(0, dtype=int)
        self.type = np.zeros(0, dtype=int)

        # Field value
        self.fval = np.zeros(0, dtype=float)

        # Boundary vertex
        self.BV = np.zeros(0, dtype=int)

    # Initialization
    def Initialize(self, ntot):
        # total number of vertices
        self.ntot = ntot

        # Coordinates
        self.x = np.zeros(ntot, dtype=float)
        self.y = np.zeros(ntot, dtype=float)

        # ID
        self.ID = np.zeros(ntot, dtype=int)
        self.type = np.zeros(ntot, dtype=int)

        # Field value
        self.fval = np.zeros(ntot, dtype=float)

        # Boundary vertex
        self.BV = np.zeros(ntot, dtype=int)


# Topological mesh face data
class TopomeshFacedata: 
    # Definition
    def __init__(self):
        # coordinates
        self.x = np.zeros(0, dtype=float)
        self.y = np.zeros(0, dtype=float)

    # Initializer
    def Initialize(self, nc):
        # Initialize coordinates
        self.x = np.zeros(nc, dtype=float)
        self.y = np.zeros(nc, dtype=float)

# Topological mesh faces
class TopomeshFace:
    # Definition
    def __init__(self):
        # total number of faces
        self.ntot = 0

        # Coordinates
        self.data = TopomeshFacedata()

        # ID
        self.ID = np.zeros(0, dtype=int)

        # Field value
        self.fval = np.zeros(0, dtype=float)

        # Boundary face
        self.BF = np.zeros(0, dtype=int)

        # Flux surface ID
        self.fsID = np.zeros(0, dtype=int)
        
        # Type
        self.type = np.zeros(0, dtype=int)
        self.vert = np.zeros((0, 2), dtype=int)

        # Number of coordinates
        self.nc = np.zeros(0, dtype=int)

    # Initialization
    def Initialize(self, ntot):
        # total number of faces
        self.ntot = ntot

        # Coordinates
        self.data = [TopomeshFacedata() for i in range(ntot)]

        # ID
        self.ID = np.zeros(ntot, dtype=int)

        # Field value
        self.fval = np.zeros(ntot, dtype=float)

        # Boundary face
        self.BF = np.zeros(ntot, dtype=int)

        # Flux surface ID
        self.fsID = np.zeros(ntot, dtype=int)
        
        # Type
        self.type = np.zeros(ntot, dtype=int)
        self.vert = np.zeros((ntot, 2), dtype=int)

        # Number of coordinates
        self.nc = np.zeros(ntot, dtype=int)


    def AddFaceCoordinates(self, faceindex, xf, yf):
        # Add face coordinates
        self.data[faceindex].x = xf 
        self.data[faceindex].y = yf

# Topological mesh cells
class TopomeshCell:
    # Definition
    def __init__(self):
        # Total number of cells
        self.ntot = 0
        self.nvert = 0
        self.nface = 0

        # ID
        self.ID = np.zeros(0, dtype=int)

        # Number of coordinates
        self.nc = np.zeros(0, dtype=int)

        # Vertices
        self.vert = np.zeros(0, dtype=int)
        self.vertP = np.zeros((0, 2), dtype=int)

        # Faces
        self.face = np.zeros(0, dtype=int)
        self.faceP = np.zeros((0, 2), dtype=int)

        # Coordinates
        self.data = TopomeshFacedata() 

    # Initialization
    def Initialize(self, ntot, nvert, nface):
        # Total number of cells
        self.ntot = ntot
        self.nvert = nvert 
        self.nface = nface

        # ID
        self.ID = np.zeros(ntot, dtype=int)

        # Number of polygon coordinates
        self.nc = np.zeros(ntot, dtype=int)

        # Vertices
        self.vert = np.zeros(nvert, dtype=int)
        self.vertP = np.zeros((ntot, 2), dtype=int)

        # Faces
        self.face = np.zeros(nface, dtype=int)
        self.faceP = np.zeros((ntot, 2), dtype=int)

        # Coordinates
        self.data = [TopomeshFacedata() for i in range(ntot)]

    def AddCellCoordinates(self, cellindex, xc, yc):
        # Add face coordinates
        self.data[cellindex].x = xc 
        self.data[cellindex].y = yc
        
# Topological mesh
class Topomesh:
    def __init__(self):
        # Fields
        self.vert = TopomeshVert()
        self.face = TopomeshFace() 
        self.cell = TopomeshCell()
        
#----------------------------------------------------------------------#
#                         GRID GENERATOR                               #
#----------------------------------------------------------------------#

# GGTM lines
class GGTMLine:
    # Definition
    def __init__(self):
        # coordinates
        self.x = np.zeros(0, dtype=float)
        self.y = np.zeros(0, dtype=float)

        # Vertices
        self.vert = np.zeros(0, dtype=int)

    # Initializer
    def Initialize(self, nv):
        # Initialize coordinates
        self.x = np.zeros(nv, dtype=float)
        self.y = np.zeros(nv, dtype=float)

        # Vertices
        self.vert = np.zeros(nv, dtype=int)

    # Coordinate addition
    def AddCoordinates(self, xc, yc, vc):
        self.x = xc 
        self.y = yc 
        self.vert = vc

# GGTM faces
class GGTMFace:
    # Definition
    def __init__(self):
        # coordinates
        self.x = np.zeros(0, dtype=float)
        self.y = np.zeros(0, dtype=float)

        # Vertices
        self.vert = np.zeros(0, dtype=int)

        # Number
        self.nv = 0

        # ID
        self.ID = 0

    # Initializer
    def Initialize(self, nv, ID):
        # Initialize coordinates
        self.x = np.zeros(nv, dtype=float)
        self.y = np.zeros(nv, dtype=float)

        # Vertices
        self.vert = np.zeros(nv, dtype=int)

        # Number
        self.nv = nv

        # ID
        self.ID = ID

# GGTM cells
class GGTMCell:
    # Definition
    def __init__(self):

        # hfline, lfline
        self.hfline = GGTMLine()
        self.lfline = GGTMLine()

        # Start and end radial face IDs
        self.srf = 0
        self.erf = 0

        # Lines
        self.nl = 0
        self.lines = GGTMLine()

        # ID
        self.ID = 0

    # Initialization
    def Initialize(self, nl, ID):

        # Coordinates
        self.nl = nl
        self.lines = [GGTMLine() for i in range(nl)]

        # ID
        self.ID = ID

    # Adding line coordinates
    def AddLineCoordinates(self, lineindex, xl, yl, vl):
        # Add face coordinates
        self.lines[lineindex].AddCoordinates(xl, yl, vl)

# GGTM data class
class GGTMData:
    # Init
    def __init__(self):
        # Fields
        self.face = GGTMFace() 
        self.cell = GGTMCell()

    # Initializer
    def Initialize(self, nf, nc):
        # Fields
        self.face = [GGTMFace() for i in range(nf)]
        self.cell = [GGTMCell() for i in range(nc)]
        
# Grid vertices
class GGVert:
    # Definition
    def __init__(self):
        # Number
        self.ntot = 0

        # coordinates
        self.x = np.zeros(0, dtype=float)
        self.y = np.zeros(0, dtype=float)

        # ID
        self.ID =  np.zeros(0, dtype=int)

        # Fieldline ID
        self.fieldlineID = np.zeros(0, dtype=int)

    # Initializer
    def Initialize(self, nv):
        # Number
        self.ntot = nv 

        # Coordinates
        self.x = np.zeros(nv, dtype=float)
        self.y = np.zeros(nv, dtype=float)

        # ID
        self.ID =  np.zeros(nv, dtype=int)

        # Fieldline ID
        self.fieldlineID = np.zeros(nv, dtype=int)

# Grid faces
class GGFace:
    # Definition
    def __init__(self):
        # Number
        self.ntot = 0 

        # Vertices
        self.v1 = np.zeros(0, dtype=int)
        self.v2 = np.zeros(0, dtype=int)

        # ID
        self.ID =  np.zeros(0, dtype=int)

        # Label
        self.label = np.zeros(0, dtype=int)
        self.region = np.zeros(0, dtype=int)

    # Initializer
    def Initialize(self, nf):
        # Number
        self.ntot = nf

        # Vertices
        self.v1 = np.zeros(nf, dtype=int)
        self.v2 = np.zeros(nf, dtype=int)

        # ID
        self.ID =  np.zeros(nf, dtype=int)

        # Label
        self.label = np.zeros(nf, dtype=int)
        self.region = np.zeros(nf, dtype=int)

# Grid cells
class GGCell:
    # Definition
    def __init__(self):
        # Number
        self.ntot = 0
        self.nvert = 0

        # Vertex pointer
        self.vp1 = np.zeros(0, dtype=int)
        self.vp2 = np.zeros(0, dtype=int)
        
        # Vertices
        self.vert = np.zeros(0, dtype=int)

        # ID
        self.ID =  np.zeros(0, dtype=int)

        # Region
        self.region = np.zeros(0, dtype=int)

    # Initializer
    def Initialize(self, nc, ncv):
        # Number
        self.ntot = nc
        self.nvert = ncv

        # Vertex pointer
        self.vp1 = np.zeros(nc, dtype=int)
        self.vp2 = np.zeros(nc, dtype=int)
        
        # Vertices
        self.vert = np.zeros(ncv, dtype=int)

        # ID
        self.ID =  np.zeros(nc, dtype=int)

        # Region
        self.region = np.zeros(nc, dtype=int)

    # Vertex getter
    def GetVert(self, i):
        return self.vert[self.vp1[i]:self.vp1[i]+self.vp2[i]]

# Grid
class GGGrid:
    # Init
    def __init__(self):
        # Fields
        self.vert = GGVert()
        self.face = GGFace() 
        self.cell = GGCell()

#----------------------------------------------------------------------#
#                        SIMULATION GRID                               #
#----------------------------------------------------------------------#
        
# Grid vertices
class Vert:
    # Definition
    def __init__(self):
        # Number
        self.ntot = 0

        # Coordinates
        self.x = np.zeros(0, dtype=float)
        self.y = np.zeros(0, dtype=float)

        # Magnetic field
        self.psi = np.zeros(0, dtype=float)
        self.bxv = np.zeros(0, dtype=float)
        self.byv = np.zeros(0, dtype=float)
        self.ffbz = np.zeros(0, dtype=float)

        # ID
        self.ID =  np.zeros(0, dtype=int)

        # Fieldline ID
        self.fieldlineID = np.zeros(0, dtype=int)

    # Initializer
    def Initialize(self, nv):
        # Number
        self.ntot = nv 

        # Coordinates
        self.x = np.zeros(nv, dtype=float)
        self.y = np.zeros(nv, dtype=float)

        # Magnetic field
        self.psi = np.zeros(nv, dtype=float)
        self.bx = np.zeros(nv, dtype=float)
        self.by = np.zeros(nv, dtype=float)
        self.ffbz = np.zeros(nv, dtype=float)

        # ID
        self.ID =  np.zeros(nv, dtype=int)

        # Fieldline ID
        self.fieldlineID = np.zeros(nv, dtype=int)

# Grid faces
class Face:
    # Definition
    def __init__(self):
        # Number
        self.ntot = 0 

        # Vertices
        self.v1 = np.zeros(0, dtype=int)
        self.v2 = np.zeros(0, dtype=int)

        # Labels
        self.label = np.zeros(0, dtype=int)
        self.region = np.zeros(0, dtype=int)
        
        # Magnetic field
        self.aligned = np.zeros(0, dtype=int)

        # ID
        self.ID =  np.zeros(0, dtype=int)

    # Initializer
    def Initialize(self, nf):
        # Number
        self.ntot = nf

        # Vertices
        self.v1 = np.zeros(nf, dtype=int)
        self.v2 = np.zeros(nf, dtype=int)

        # Labels
        self.label = np.zeros(nf, dtype=int)
        self.region = np.zeros(nf, dtype=int)
        
        # Magnetic field
        self.aligned = np.zeros(nf, dtype=int)

        # ID
        self.ID =  np.zeros(nf, dtype=int)

# Grid cells
class Cell:
    # Definition
    def __init__(self):
        # Number
        self.ntot = 0
        self.nvert = 0
        self.nface = 0

        # Vertex pointer
        self.vp1 = np.zeros(0, dtype=int)
        self.vp2 = np.zeros(0, dtype=int)

        # Face pointer
        self.fp1 = np.zeros(0, dtype=int)
        self.fp2 = np.zeros(0, dtype=int)
        
        # Vertices
        self.vert = np.zeros(0, dtype=int)

        # Faces
        self.face = np.zeros(0, dtype=int)

        # Coordinates
        self.x  = np.zeros(0, dtype=float)
        self.y  = np.zeros(0, dtype=float)

        # Magnetic Field
        self.bt = np.zeros(0, dtype=float)
        self.bp = np.zeros(0, dtype=float)
        self.psi   = np.zeros(0, dtype=float)

        # ID
        self.ID =  np.zeros(0, dtype=int)

        # Region etc
        self.cflags = np.zeros(0, dtype=int)
        self.region = np.zeros(0, dtype=int)
        self.ft     = np.zeros(0, dtype=int)

    # Initializer
    def Initialize(self, nc, ncv, ncf):
        # Number
        self.ntot = nc
        self.nvert = ncv
        self.nface = ncf

        # Vertex pointer
        self.vp1 = np.zeros(nc, dtype=int)
        self.vp2 = np.zeros(nc, dtype=int)

        # Face pointer
        self.fp1 = np.zeros(nc, dtype=int)
        self.fp2 = np.zeros(nc, dtype=int)
        
        # Vertices
        self.vert = np.zeros(ncv, dtype=int)

        # Faces
        self.face = np.zeros(ncf, dtype=int)

        # Coordinates
        self.x  = np.zeros(nc, dtype=float)
        self.y  = np.zeros(nc, dtype=float)

        # Magnetic Field
        self.bt = np.zeros(nc, dtype=float)
        self.bp = np.zeros(nc, dtype=float)
        self.psi   = np.zeros(nc, dtype=float)

        # ID
        self.ID =  np.zeros(nc, dtype=int)

        # Region etc
        self.cflags = np.zeros(nc, dtype=int)
        self.region = np.zeros(nc, dtype=int)
        self.ft     = np.zeros(nc, dtype=int)

    # Vertex getter
    def GetVert(self, i):
        return self.vert[self.vp1[i]:self.vp1[i]+self.vp2[i]]
    
    # Face getter
    def GetFace(self, i):
        return self.face[self.fp1[i]:self.fp1[i]+self.fp2[i]]

# Grid flux surfaces
class FluxSurf:
    # Definition
    def __init__(self):
        # Number
        self.ntot = 0
        self.nface = 0

        # Faces
        self.fp1 = np.zeros(0, dtype=int)
        self.fp2 = np.zeros(0, dtype=int)
        self.psi = np.zeros(0, dtype=float)
        self.face = np.zeros(0, dtype=int)

        # ID
        self.ID = np.zeros(0, dtype=int)

    # Initializer
    def Initialize(self, nfs):
        # Number
        self.ntot = nfs 

        # Faces
        self.fp1 = np.zeros(nfs, dtype=int)
        self.fp2 = np.zeros(nfs, dtype=int)

        # ID
        self.ID = np.zeros(nfs, dtype=int)

        # Magnetic field
        self.psi = np.zeros(nfs, dtype=float)

    def InitializeFaceData(self, nfsf):
        self.nface = nfsf 
        self.face = np.zeros(nfsf, dtype=int)

    # Face getter
    def GetFace(self, i):
        return self.face[self.fp1[i]:self.fp1[i]+self.fp2[i]]
    
# Grid flux tubes
class FluxTube:
    # Definition
    def __init__(self):
        # Number
        self.ntot = 0
        self.nface = 0
        self.ncell = 0

        # Faces
        self.fp1 = np.zeros(0, dtype=int)
        self.fp2 = np.zeros(0, dtype=int)
        self.face = np.zeros(0, dtype=int)

        # Cells
        self.cp1 = np.zeros(0, dtype=int)
        self.cp2 = np.zeros(0, dtype=int)
        self.cell = np.zeros(0, dtype=int)

        # Region
        self.region = np.zeros(0, dtype=int)

        # ID
        self.ID = np.zeros(0, dtype=int)

    # Initializer
    def Initialize(self, nft):
        # Number
        self.ntot = nft

        # Faces
        self.fp1 = np.zeros(nft, dtype=int)
        self.fp2 = np.zeros(nft, dtype=int)

        # Cells
        self.cp1 = np.zeros(nft, dtype=int)
        self.cp2 = np.zeros(nft, dtype=int)

        # Region
        self.region = np.zeros(nft, dtype=int)

        # ID
        self.ID = np.zeros(nft, dtype=int)

    def InitializeFaceData(self, nftf):
        self.nface = nftf
        self.face = np.zeros(nftf, dtype=int)

    def InitializeCellData(self, nftc):
        self.ncell = nftc
        self.cell = np.zeros(nftc, dtype=int)

    # Face getter
    def GetFace(self, i):
        return self.face[self.fp1[i]:self.fp1[i]+self.fp2[i]]
    
    # Cell getter
    def GetFace(self, i):
        return self.cell[self.cp1[i]:self.cp1[i]+self.cp2[i]]

# Grid
class Grid:
    # Init
    def __init__(self):
        # Fields
        self.vert = Vert()
        self.face = Face() 
        self.cell = Cell()
        self.fs = FluxSurf()
        self.ft = FluxTube()

    # Interconnections
    def ComputeInterconnections(self):
        # Description
        #------------
        # Compute grid interconnections derived from basic quantities 
        # that are read in (cell faces etc)
