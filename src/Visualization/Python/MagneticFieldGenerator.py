# Description
#------------
# Simple script to generate our own magnetic fields based on analytic
# descriptions. Writes out an rzpsi file using the datahandler 
# module
import Datahandler as dh 
import Plotter as pl
import numpy as np

# Define magnetic field
def MagneticField(x, y):
    psi = np.sin(2*np.pi*x) + np.sin(2*np.pi*y)
    return psi

writedir = './goatf/src/Visualization'


# Construct R, Z coordinates
resx = 100
resy = 200
Lx = 3.4
Ly = 4.4
Lxoffset = -0.2
Lyoffset = -0.2
R = (np.arange(0, resx+1, 1))/float(resx)*Lx + Lxoffset 
Z = (np.arange(0, resy+1, 1))/float(resy)*Ly + Lyoffset 
nR = resx + 1
nZ = resy + 1
Psi = np.zeros([nR, nZ])

# Compute
for j in np.arange(0, nZ, 1):
    for i in np.arange(0, nR, 1):
        Psi[i, j] = MagneticField(R[i], Z[j])

# Visualize
pl.PlotStructured2DContourf(R, Z, Psi, 1)
pl.ShowFigures()

# Write
dh.WriteRZPsiFile(writedir, R, Z, Psi)





