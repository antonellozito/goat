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

        # Vertices
        self.vert = np.zeros(0, dtype=int)
        self.vertP = np.zeros((0, 2), dtype=int)

        # Faces
        self.face = np.zeros(0, dtype=int)
        self.faceP = np.zeros((0, 2), dtype=int)

    # Initialization
    def Initialize(self, ntot, nvert, nface):
        # Total number of cells
        self.ntot = ntot

        # Vertices
        self.vert = np.zeros(nvert, dtype=int)
        self.vertP = np.zeros((ntot, 2), dtype=int)

        # Faces
        self.face = np.zeros(nface, dtype=int)
        self.faceP = np.zeros((ntot, 2), dtype=int)

        
# Topological mesh
class Topomesh:
    def __init__(self):
        # Fields
        self.vert = TopomeshVert()
        self.face = TopomeshFace() 
        self.cell = TopomeshCell()

        

