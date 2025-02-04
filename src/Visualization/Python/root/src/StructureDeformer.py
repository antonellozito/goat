# Description
#------------
# Simple module to deform a structure and to move points etc

import numpy as np
# Rotation around point
def RotatePoint(px, py, ox, oy, theta):
    # Description
    #------------
    # Rotate the vector OP for theta amount of degrees around point O. 
    # Returns the new coordinates of the rotated point P. Theta should 
    # be in degrees

    # Initialize
    dx = px - ox 
    dy = py - oy 
    theta = theta/180*np.pi # convert to radians 

    # Rotate
    newpx = ox + dx*np.cos(theta) - dy*np.sin(theta)
    newpy = oy + dx*np.sin(theta) + dy*np.cos(theta)

    return newpx, newpy

# Translation along a vector
def TranslatePoint(px, py, dx, dy, length):
    # Description
    #------------
    # Translate a point by a length 'length' along the vector dx, dy

    # Initialize
    dn = np.sqrt(dx**2 + dy**2)
    dx = dx/dn 
    dy = dy/dn 

    # Translate
    newpx = px + dx*length 
    newpy = py + dy*length 

    return newpx, newpy