import sys
import os
 
script_dir = os.path.dirname(os.path.abspath(sys.argv[0]))  # Get the script's directory
sys.path.append(script_dir)  # Add Visualization/ to sys.path

import GOATpy as gp

thisdir = os.getcwd()

Rpath = thisdir + '/R.csv'
Zpath = thisdir + '/Z.csv'
Psipath = thisdir + '/psi.csv'
separator = ','
[R, Z, PsiWb, nZ, nR] = gp.dh.ReadRZPsiFromCSV(Rpath, Zpath, Psipath, separator)
Psi = PsiWb # / 1000.0 # mWb to Wb
gp.dh.WriteRZPsiFile(thisdir, R, Z, Psi)
