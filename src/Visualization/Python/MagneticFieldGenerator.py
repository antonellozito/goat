# Description
#------------
# Simple script to generate our own magnetic fields based on analytic
# descriptions. Writes out an rzpsi file using the datahandler 
# module
import Datahandler as dh 
import Plotter as pl
import numpy as np
import goat_types as gt

# Define magnetic field
def MagneticField(x, y):
    #psi = np.sin(2*np.pi*x) + np.sin(2*np.pi*y)
    #return psi
    Btor = 0.08666
    # psi = np.sqrt((x - 0.78)**2 + y**2) 
    psi = 0.001*Btor*x + 0.0000001*Btor*y
    return psi

    

writedir = './goatf/src/Visualization'


# Construct R, Z coordinates
resx = 400
resy = 400
Lx = 4.0
Ly = 4.0
Lxoffset = -0.0
Lyoffset = -2.0
R = (np.arange(0, resx+1, 1))/float(resx)*Lx + Lxoffset 
Z = (np.arange(0, resy+1, 1))/float(resy)*Ly + Lyoffset 
nR = resx + 1
nZ = resy + 1
Psi = np.zeros([nR, nZ])

Rmaj = 0.78 
Rmin = 0.26
theta = np.arange(0, 101, 1)/100.0*2*np.pi
xv = Rmaj + Rmin*np.cos(theta)
yv = Rmin*np.sin(theta)

# Construct structures
structures = [gt.Structure() for i in np.arange(0, 1, 1)]
structures[0].Initialize(len(xv), xv, yv)

# Write output
dh.WriteStructureFile(writedir, structures)

# Compute
for j in np.arange(0, nZ, 1):
    for i in np.arange(0, nR, 1):
        Psi[i, j] = MagneticField(R[i], Z[j])

# Visualize
pl.PlotStructured2DContourf(R, Z, Psi, 1)
pl.ShowFigures()

# Write
dh.WriteRZPsiFile(writedir, R, Z, Psi)





