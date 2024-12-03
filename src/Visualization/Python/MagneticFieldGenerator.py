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
    #psi = np.sin(2*np.pi*x) + np.sin(2*np.pi*y)
    #return psi
    psi = np.sqrt((x - 0.78)**2 + y**2) 
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

# Compute
for j in np.arange(0, nZ, 1):
    for i in np.arange(0, nR, 1):
        Psi[i, j] = MagneticField(R[i], Z[j])

# Visualize
pl.PlotStructured2DContourf(R, Z, Psi, 1)
pl.ShowFigures()

# Write
dh.WriteRZPsiFile(writedir, R, Z, Psi)





