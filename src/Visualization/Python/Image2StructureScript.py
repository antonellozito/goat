# Description
#------------
# This script converts an image (path to be defined by the user) to 
# a set of contours that are written into a structure.dat file. Each 
# separate contour forms a separate structure in that file.
# Note that it is possible that a contour is generated at the edge of
# the figure. This may be manually removed. In order to rescale the 
# contours, the variables Lx, Ly, Loffsetx, Loffsety should be defined

# Load modules
import numpy as  np 
import cv2
import Datahandler as dh
import goat_types as gt

# Define inputs
impath = '/mnt/c/Users/u0110555/Desktop/code_werk/goatf/goatf/src/Visualization/PET.JPG' # image path
structuredir = '/mnt/c/Users/u0110555/Desktop/code_werk/goatf/goatf/src/Visualization' 
contourthreshold = 150
Lx = 3.0
Loffsetx = 0.0 
Ly = 4.0
Loffsety = 0.0

# Load image
img  = cv2.imread(impath)
assert img is not None, "file could not be read, check with os.path.exists()"

# Convert to gray scale
imgray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
ret, thresh = cv2.threshold(imgray, contourthreshold, 255, 0)
#cv2.imshow('binary image', thresh)
#cv2.waitKey(0)

# Compute contours - note, these should be closed by default! 
contours, hierarchy = cv2.findContours(thresh, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)
imgc = img
#cv2.drawContours(imgc, contours, -1, (0, 255, 0), 2)
#cv2.imshow('binary image', imgc)
#cv2.waitKey(0)

# Construct structures
structures = [gt.Structure() for i in np.arange(0, len(contours), 1)]
for j in np.arange(0, len(contours), 1):
    x = np.zeros(contours[j].shape[0]+1, dtype=float)
    y = np.zeros(contours[j].shape[0]+1, dtype=float)
    for k in np.arange(0, contours[j].shape[0], 1):
        x[k] = contours[j][k, 0][0]
        y[k] = contours[j][k, 0][1]
    # Need to add endpoint to close
    x[len(x)-1] = x[0]
    y[len(y)-1] = y[0]

    # mirror over x-axis
    ynew = -y

    # Shift to original position
    ynew = ynew + img.shape[0]

    y = ynew

    # Rescale
    x = x/float(img.shape[1]-1)*Lx + Loffsetx 
    y = y/float(img.shape[0]-1)*Ly + Loffsety
    structures[j].Initialize(contours[j].shape[0]+1, x, y)

# Write output
dh.WriteStructureFile(structuredir, structures)
