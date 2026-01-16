# Description
#------------
# This script converts the data available in an eqdsk file to a rzpsi.dat
# file. It is basically a script wrapper for the Datahandler routine 
# ReadRZPsiFromEqdskFile. It can take two input arguments: the first one
# should be the path to the eqdsk file. The second, optional, argument
# is then the path to where the output needs to be written (including
# the filename). If empty, the output is written to the current file in rzpsi.dat.

import GOATpy as gp 
import os
import sys 

# Check if names were given as input (first argument is python script name, 
# second is assumed to be gridname, third is topomesh name)
narg = len(sys.argv)
print("total number of arguments passed: ", narg)
if narg < 2:
    raise ValueError('ConverteqdsktoRZPsi: at least one input argument with the file location should be given')

# Set defaults
rzpsifile = './rzpsi.dat'

for i in range(1, narg):
    # Check
    if (i == 1):
        eqdskfile = sys.argv[i]
    elif (i == 2): 
        rzpsifile = sys.argv[i]

# Print
print('ConverteqdsktoRZPsi: reading eqdsk from file: ' + eqdskfile)
print('ConverteqdsktoRZPsi: writing rzpsi to file: ' + rzpsifile)

# Read
R, Z, Psi = gp.dh.ReadRZPsiFromEqdskFile(eqdskfile)

# Write 
gp.dh.WriteRZPsiFile(rzpsifile, R, Z, Psi)