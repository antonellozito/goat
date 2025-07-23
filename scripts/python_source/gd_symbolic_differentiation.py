# Description
#------------
# This is a simple script that can be used to analytically differentiate
# cost function or constraint contributions for the grid deformation
# algorithm (or other algorithms, if required). There are some 
# predefined cases for which first and second order derivatives are 
# computed, which one can use as a basis to construct one's own 
# derivatives. This script builds on the SymPy symbolic toolbox for 
# Python (one could do similar things in Matlab). 

# The aim is that this script provides copy-pasteable expressions for 
# first and second (and if really necessary even higher-order) derivatives
# that can be used directly in fortran. Variable naming and substition 
# is therefore important!

# It is assumed that each case defines the following:
# - diffsym: a list of symbolic variables with respect which to differentiate
# - Ji: the cost function contribution to be differentiated
# - subfrom: an expression that has to be substituted in the differentiated output
# - subto: the expression with which to substitute
# note that the order in subfrom/subto matters, since these substitutions
# are performed in sequence!

# Import libraries
import sympy 

# Initialize printer
sympy.init_printing()

# Define case
casename = 'length_distribution'

# Define output filepath
writefilepath = 'derivatives_' + casename + '.dat'

# Setup
#------
match casename:

    case 'cell_angles':

        # Cell angle cost function, where each contribution is defined
        # as the angle between two faces (they don't have to have 
        # a common vertex) and therefore four vertices with coordinates
        # x1, y1, x2, y2, x3, y3, x4, y4. The cost function contribution
        # of two faces is then:
        #
        #   Ji = 0.5*wti*(theta - theta0)**2
        #   theta = atan(cp/dp)
        #   dp = dx1*dx2 + dy1*dy2
        #   cp = dx1*dy2 - dy1*dx2
        #   dx1 = x2 - x1, dy1 = y2 - y1
        #   dx2 = x4 - x3, dy2 = y4 - y3

        # Define symbolic variables
        x1, x2, x3, x4, y1, y2, y3, y4, theta0, wti = sympy.symbols('x1, x2, x3, x4, y1, y2, y3, y4, theta0, wti')
        dx1, dx2, dy1, dy2, cp, dp, theta, Ji = sympy.symbols('dx1, dx2, dy1, dy2, cp, dp, theta, Ji')

        # Define symbols to differentiate to
        diffsym = [x1, x2, x3, x4, y1, y2, y3, y4]

        # Define substition lists (order matters!)
        subsfrom = [x2 - x1, x4 - x3, y2 - y1, y4 - y3, 
            dx1*dx2 + dy1*dy2, dx1*dy2 - dy1*dx2, sympy.atan2(cp, dp)]
        substo = [dx1, dx2, dy1, dy2, dp, cp, theta]

        # Define derived quantities
        dx1 = x2 - x1 
        dx2 = x4 - x3 
        dy1 = y2 - y1
        dy2 = y4 - y3 
        dp = dx1*dx2 + dy1*dy2
        cp = dx1*dy2 - dy1*dx2
        theta = sympy.atan2(cp, dp)
        Ji = 0.5*wti*(theta - theta0)**2

    case 'length_distribution':

        # Length (normalized and weighted) cost function 
        # contribution for each face (best to be used only for aligned
        # faces). Cost function contribution is defined as:
        #
        # Ji = 0.5*wt(xf, yf)*( (li - l0(xf, yf))/l0(xf, yf) )**2
        # li = sqrt(dx**2 + dy**2)
        # dx = x2 - x1 
        # dy = y2 - y1
        # xf = 0.5*(x1 + x2)
        # yf = 0.5*(y1 + y2)
        # l0(xf, yf): some unknown function for the length distribution,
        # but for which the first and second order derivatives exist 
        # and are finite. 

        # Define symbolic variables
        x1, x2, y1, y2, dx, dy, li, Ji, l0i, wti, rat= sympy.symbols('x1, x2, y1, y2, dx, dy, li, Ji, l0i, wti, rat')
        dl0dx, dl0dy, d2l0dx2, d2l0dxdy, d2l0dy2 = sympy.symbols('dl0dx, dl0dy, d2l0dx2, d2l0dxdy, d2l0dy2')
        dwtdx, dwtdy, d2wtdx2, d2wtdxdy, d2wtdy2 = sympy.symbols('dwtdx, dwtdy, d2wtdx2, d2wtdxdy, d2wtdy2')

        # Define undefined functions
        l0 = sympy.Function('l0')(x1, y1, x2, y2)
        wt = sympy.Function('wt')(x1, y1, x2, y2)

        # Define symbols to differentiate to
        diffsym = [x1, x2, y1, y2]

        # Define substition lists (order matters!)
        subsfrom = [x2 - x1, y2 - y1, sympy.sqrt(dx**2 + dy**2), 
            sympy.diff(l0, x1, x1), sympy.diff(l0, x2, x1), sympy.diff(l0, y1, x1), sympy.diff(l0, y2, x1), 
            sympy.diff(l0, x1, x2), sympy.diff(l0, x2, x2), sympy.diff(l0, y1, x2), sympy.diff(l0, y2, x2), 
            sympy.diff(l0, x1, y1), sympy.diff(l0, x2, y1), sympy.diff(l0, y1, y1), sympy.diff(l0, y2, y1), 
            sympy.diff(l0, x1, y2), sympy.diff(l0, x2, y2), sympy.diff(l0, y1, y2), sympy.diff(l0, y2, y2), 
            sympy.diff(l0, x1), sympy.diff(l0, x2), sympy.diff(l0, y1), sympy.diff(l0, y2), 
            sympy.diff(wt, x1, x1), sympy.diff(wt, x2, x1), sympy.diff(wt, y1, x1), sympy.diff(wt, y2, x1), 
            sympy.diff(wt, x1, x2), sympy.diff(wt, x2, x2), sympy.diff(wt, y1, x2), sympy.diff(wt, y2, x2), 
            sympy.diff(wt, x1, y1), sympy.diff(wt, x2, y1), sympy.diff(wt, y1, y1), sympy.diff(wt, y2, y1), 
            sympy.diff(wt, x1, y2), sympy.diff(wt, x2, y2), sympy.diff(wt, y1, y2), sympy.diff(wt, y2, y2), 
            sympy.diff(wt, x1), sympy.diff(wt, x2), sympy.diff(wt, y1), sympy.diff(wt, y2), 
            l0, wt, (li - l0i)/l0i]
        substo = [dx, dy, li, 
            0.25*d2l0dx2, 0.25*d2l0dx2, 0.25*d2l0dxdy, 0.25*d2l0dxdy, 
            0.25*d2l0dx2, 0.25*d2l0dx2, 0.25*d2l0dxdy, 0.25*d2l0dxdy, 
            0.25*d2l0dxdy, 0.25*d2l0dxdy, 0.25*d2l0dy2, 0.25*d2l0dy2,
            0.25*d2l0dxdy, 0.25*d2l0dxdy, 0.25*d2l0dy2, 0.25*d2l0dy2,
            0.5*dl0dx, 0.5*dl0dx, 0.5*dl0dy, 0.5*dl0dy, 
            0.25*d2wtdx2, 0.25*d2wtdx2, 0.25*d2wtdxdy, 0.25*d2wtdxdy, 
            0.25*d2wtdx2, 0.25*d2wtdx2, 0.25*d2wtdxdy, 0.25*d2wtdxdy, 
            0.25*d2wtdxdy, 0.25*d2wtdxdy, 0.25*d2wtdy2, 0.25*d2wtdy2, 
            0.25*d2wtdxdy, 0.25*d2wtdxdy, 0.25*d2wtdy2, 0.25*d2wtdy2, 
            0.5*dwtdx, 0.5*dwtdx, 0.5*dwtdy, 0.5*dwtdy, 
            l0i, wti, rat]

        # Define derived quantities
        dx = x2 - x1 
        dy = y2 - y1 
        li = sympy.sqrt(dx**2 + dy**2)
        Ji = 0.5*wt*((li - l0)/l0)**2

    case 'exp_distance':

        # Define symbolic variables
        x, y, xa, ya, c, d, d0 = sympy.symbols('x, y, xa, ya, c, d, d0')
        dx1, dx2, dy1, dy2, cp, dp, theta, Ji = sympy.symbols('dx1, dx2, dy1, dy2, cp, dp, theta, Ji')

        # Define symbols to differentiate to
        diffsym = [x, y]

        # Define substition lists (order matters!)
        subsfrom = [sympy.sqrt((x - xa)**2 + (y - ya)**2)]
        substo = [d]

        # Define derived quantities
        d = sympy.sqrt((x - xa)**2 + (y - ya)**2)
        Ji = c*sympy.exp(-d/d0) 
        
    case _: 

        raise ValueError('Case not implemented')
    

# Differentiate
#--------------
# First order
dJidx = [] 
for i in diffsym:
    # Compute
    tempdiff = sympy.diff(Ji, i)

    # Substitute
    for k in range(0, len(subsfrom)):
        tempdiff = tempdiff.subs(subsfrom[k], substo[k])
    
    # Simplify
    tempdiff.simplify()

    # Append
    dJidx.append(tempdiff)

    # Print
    print("dJid" + str(i) + " = " + str(tempdiff))

# Second order
d2Jdidx2 = []
for i in diffsym:
    for j in diffsym:
        # Compute
        tempdiff = sympy.diff(Ji, i, j)

        # Substitute
        for k in range(0, len(subsfrom)):
            tempdiff = tempdiff.subs(subsfrom[k], substo[k])
        
        # Simplify
        tempdiff.simplify()

        # Append
        d2Jdidx2.append(tempdiff)

        # Print
        print("d2Jid" + str(i) + "d" + str(j) + " = " + str(tempdiff))

# Write
#------
# Open file
with open(writefilepath, 'w') as f:
    # First order derivatives
    for i in range(0, len(diffsym)):
        # Print
        print("dJid" + str(diffsym[i]) + " = " + str(dJidx[i]) + " !" + str(diffsym[i]), file=f)

    # Second order derivatives
    k = 0
    for i in range(0, len(diffsym)):
        for j in range(0, len(diffsym)):
            # Print
            print("d2Jid" + str(diffsym[i]) + "d" + str(diffsym[j]) + " = " + str(d2Jdidx2[k])
                + " !" + str(diffsym[i]) + str(diffsym[j]), file=f)
            k = k + 1

