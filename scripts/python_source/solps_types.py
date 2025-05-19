# This module defines some convenient types for post-processing that
# mimick the types of solps (to some extent). All in unstructured 
# format. 
import numpy as np 
from . import Datahandler as dh 
from scipy.interpolate import griddata 
from scipy import interpolate
import time 


#----------------------------------------------------------------------#
#                           CONSTANTS                                  #
#----------------------------------------------------------------------#

class PhysicalConstants:
    def __init__(self):
        self.eVperJoule = 6.241509e+18

#----------------------------------------------------------------------#
#                           GEOMETRY                                   #
#----------------------------------------------------------------------#

#----------------------------------------------------------------------#
#                             STATE                                    #
#----------------------------------------------------------------------#

class PlasmaState:
    def __init__(self):
        # Description
        #------------
        # Default initializer

        # Header
        self.header = " "
        self.label = " "

        # Dimensions
        self.ns     = 0 # number of species
        self.nCv    = 0 # number of cells 
        self.nFc    = 0 # number of faces 

        # Charges
        self.zamin  = np.zeros(self.ns, dtype=float)
        self.zamax  = np.zeros(self.ns, dtype=float)
        self.zn     = np.zeros(self.ns, dtype=float)
        self.am     = np.zeros(self.ns, dtype=float)

        # State variables
        self.na     = np.zeros((self.nCv, self.ns), dtype=float)
        self.ne     = np.zeros(self.nCv, dtype=float)
        self.ua     = np.zeros((self.nCv, self.ns), dtype=float)
        self.uadia  = np.zeros(self.nCv, dtype=float)
        self.te     = np.zeros(self.nCv, dtype=float)
        self.ti     = np.zeros(self.nCv, dtype=float)
        self.tn     = np.zeros(self.nCv, dtype=float)
        self.po     = np.zeros(self.nCv, dtype=float)
        self.kt     = np.zeros(self.nCv, dtype=float)
        self.zt     = np.zeros(self.nCv, dtype=float)
        self.kinrgy = np.zeros((self.nCv, self.ns), dtype=float)

        # Fluxes 
        self.fna    = np.zeros((self.nFc, 2, self.ns), dtype=float)
        self.fne    = np.zeros((self.nFc, 2), dtype=float)
        self.fni    = np.zeros((self.nFc, 2), dtype=float)
        self.fmo    = np.zeros((self.nFc, 2, self.ns), dtype=float)
        self.fhe    = np.zeros((self.nFc, 2), dtype=float)
        self.fhi    = np.zeros((self.nFc, 2), dtype=float)
        self.fhn    = np.zeros((self.nFc, 2), dtype=float)
        self.fch    = np.zeros((self.nFc, 2), dtype=float)
        self.fch_32 = np.zeros((self.nFc, 2), dtype=float)
        self.fch_52 = np.zeros((self.nFc, 2), dtype=float)
        self.fhm    = np.zeros((self.nFc, 2, self.ns), dtype=float)
        self.fkt    = np.zeros((self.nFc, 2), dtype=float)
        self.fzt    = np.zeros((self.nFc, 2), dtype=float)
        self.fch_p  = np.zeros((self.nFc, 2), dtype=float)

        # Additional fluxes
        self.fna_mdf     = np.zeros((self.nFc, 2, self.ns), dtype=float)
        self.fhe_mdf     = np.zeros((self.nFc, 2), dtype=float)
        self.fhi_mdf     = np.zeros((self.nFc, 2), dtype=float)
        self.fna_fcor    = np.zeros((self.nFc, 2, self.ns), dtype=float)
        self.fna_nodrift = np.zeros((self.nFc, 2, self.ns), dtype=float)
        self.fna_he      = np.zeros((self.nFc, 2, self.ns), dtype=float)
        self.fnaPSch     = np.zeros((self.nFc, 2, self.ns), dtype=float)
        self.fhePSch     = np.zeros((self.nFc, 2), dtype=float)
        self.fhiPSch     = np.zeros((self.nFc, 2), dtype=float)
        self.fna_eir     = np.zeros((self.nFc, 2, self.ns), dtype=float)
        self.fne_eir     = np.zeros((self.nFc, 2), dtype=float)
        self.fhe_eir     = np.zeros((self.nFc, 2), dtype=float)
        self.fhi_eir     = np.zeros((self.nFc, 2), dtype=float)
        self.fna_32      = np.zeros((self.nFc, 2, self.ns), dtype=float)
        self.fna_52      = np.zeros((self.nFc, 2, self.ns), dtype=float)
        self.fni_32      = np.zeros((self.nFc, 2), dtype=float)
        self.fni_52      = np.zeros((self.nFc, 2), dtype=float)
        self.fne_32      = np.zeros((self.nFc, 2), dtype=float)
        self.fne_52      = np.zeros((self.nFc, 2), dtype=float)
        self.fchdia      = np.zeros((self.nFc, 2), dtype=float)
        self.fchin       = np.zeros((self.nFc, 2), dtype=float)
        self.fchvispar   = np.zeros((self.nFc, 2), dtype=float)
        self.fchvisper   = np.zeros((self.nFc, 2), dtype=float)
        self.fchvisq     = np.zeros((self.nFc, 2), dtype=float)
        self.fchviskt   = np.zeros((self.nFc, 2), dtype=float)
        self.fchinert    = np.zeros((self.nFc, 2), dtype=float)
        self.fchanml     = np.zeros((self.nFc, 2), dtype=float)
        self.fht         = np.zeros((self.nFc, 2), dtype=float)
        self.fhj         = np.zeros((self.nFc, 2), dtype=float)
        self.fhp         = np.zeros((self.nFc, 2, self.ns), dtype=float)
        
        self.vaecrb = np.zeros((self.nFc, 2, self.ns), dtype=float)
        self.vadia  = np.zeros((self.nFc, 2, self.ns), dtype=float)
        self.wadia  = np.zeros((self.nFc, 2, self.ns), dtype=float)
        self.veecrb = np.zeros((self.nFc, 2), dtype=float)
        self.vedia  = np.zeros((self.nFc, 2), dtype=float)
        
        self.floe_noc  = np.zeros((self.nFc, 2), dtype=float)
        self.floi_noc  = np.zeros((self.nFc, 2), dtype=float)

        # Time
        self.time = 0

        # Residuals
        self.resco = np.zeros((self.nCv, self.ns), dtype=float)
        self.reshe = np.zeros((self.nCv), dtype=float)
        self.reshi = np.zeros((self.nCv), dtype=float)
        self.reshn = np.zeros((self.nCv), dtype=float)
        self.resmo = np.zeros((self.nCv, self.ns), dtype=float)
        self.resmt = np.zeros((self.nCv, self.ns), dtype=float)
        self.respo = np.zeros((self.nCv, self.ns), dtype=float)
        self.reskt = np.zeros((self.nCv, self.ns), dtype=float)
        self.reszt = np.zeros((self.nCv, self.ns), dtype=float)

        # Sources
        self.sna = np.zeros((self.nCv, 2, self.ns), dtype=float)
        self.smo = np.zeros((self.nCv, 4, self.ns), dtype=float)
        self.smq = np.zeros((self.nCv, 4, self.ns), dtype=float)
        self.shi = np.zeros((self.nCv, 4), dtype=float)
        self.she = np.zeros((self.nCv, 4), dtype=float)
        self.skt = np.zeros((self.nCv, 4), dtype=float)
        self.skt_prod = np.zeros((self.nCv), dtype=float)
        self.skt_diss = np.zeros((self.nCv), dtype=float)
        
        # Coefficients
        self.calf = np.zeros((self.nFc, 2), dtype=float)
        self.cdna = np.zeros((self.nFc, 2, self.ns), dtype=float)
        self.cdpa = np.zeros((self.nFc, 2, self.ns), dtype=float)
        self.ceqp = np.zeros((self.nCv), dtype=float)
        self.chce = np.zeros((self.nFc, 2), dtype=float)
        self.chci = np.zeros((self.nFc, 2), dtype=float)
        self.chve = np.zeros((self.nFc, 2), dtype=float)
        self.chvemx = np.zeros((self.nFc), dtype=float)
        self.chvi = np.zeros((self.nFc, 2), dtype=float)
        self.chvimx = np.zeros((self.nFc), dtype=float)
        self.csig = np.zeros((self.nFc, 2), dtype=float)
        self.cvla = np.zeros((self.nFc, 2, self.ns), dtype=float)
        self.cvsa = np.zeros((self.nFc, 2, self.ns), dtype=float)
        self.cthe = np.zeros((self.nCv, self.ns), dtype=float)
        self.cthi = np.zeros((self.nCv, self.ns), dtype=float)
        self.csigin = np.zeros((self.nFc, 2, self.ns, self.ns), dtype=float)
        self.cvsa_cl = np.zeros((self.nFc, 2, self.ns), dtype=float)
        self.fllime = np.zeros((self.nFc), dtype=float)
        self.fllimi = np.zeros((self.nFc), dtype=float)
        self.fllim0fna = np.zeros((self.nFc, 2, self.ns), dtype=float)
        self.fllim0fhi = np.zeros((self.nFc, 2, self.ns), dtype=float)
        self.fllimvisc = np.zeros((self.nFc, self.ns), dtype=float)
        self.f_luc_ke = np.zeros((self.nFc), dtype=float)
        self.f_luc_ki = np.zeros((self.nFc), dtype=float)
        self.f_luc_et = np.zeros((self.nFc), dtype=float)
        self.f_luc_sg = np.zeros((self.nFc), dtype=float)
        self.f_luc_al = np.zeros((self.nFc), dtype=float)
        self.fllim_ke = np.zeros((self.nFc), dtype=float)
        self.fllim_ki = np.zeros((self.nFc), dtype=float)
        self.fllim_et = np.zeros((self.nFc), dtype=float)
        self.fllim_sg = np.zeros((self.nFc), dtype=float)
        self.fllim_al = np.zeros((self.nFc), dtype=float)
        self.fllim_al_C = np.zeros((self.nCv), dtype=float)

        self.sig0 = np.zeros((self.nCv), dtype=float)
        self.hce0 = np.zeros((self.nCv), dtype=float)
        self.alf0 = np.zeros((self.nCv), dtype=float)
        self.hci0 = np.zeros((self.nCv), dtype=float)
        self.dpa0 = np.zeros((self.nCv, self.ns), dtype=float)
        self.dna0 = np.zeros((self.nCv, self.ns), dtype=float)
        self.vsa0 = np.zeros((self.nCv, self.ns), dtype=float)
        self.vla0 = np.zeros((self.nCv, 2, self.ns), dtype=float)
        self.dkt0 = np.zeros((self.nCv), dtype=float)
        self.dna_ExB = np.zeros((self.nCv), dtype=float)
        self.hce_ExB = np.zeros((self.nCv), dtype=float)
        self.hci_ExB = np.zeros((self.nCv), dtype=float)


    def Initialize(self, nCv, nFc, ns):
        # Description
        #------------
        # Actual initializer


        # Header (default, unless read in from an existing file)
        self.header = "VERSION03.002.000"
        self.label = "b2mn"

        # Dimensions
        self.ns     = ns # number of species
        self.nCv    = nCv # number of cells 
        self.nFc    = nFc # number of faces 

        # Charges
        self.zamin  = np.zeros(self.ns, dtype=float)
        self.zamax  = np.zeros(self.ns, dtype=float)
        self.zn     = np.zeros(self.ns, dtype=float)
        self.am     = np.zeros(self.ns, dtype=float)

        # State variables
        self.na     = np.zeros((self.nCv, self.ns), dtype=float)
        self.ne     = np.zeros(self.nCv, dtype=float)
        self.ua     = np.zeros((self.nCv, self.ns), dtype=float)
        self.uadia  = np.zeros((self.nFc, 2, self.ns), dtype=float)
        self.te     = np.zeros(self.nCv, dtype=float)
        self.ti     = np.zeros(self.nCv, dtype=float)
        self.tn     = np.zeros(self.nCv, dtype=float)
        self.po     = np.zeros(self.nCv, dtype=float)
        self.kt     = np.zeros(self.nCv, dtype=float)
        self.zt     = np.zeros(self.nCv, dtype=float)
        self.kinrgy = np.zeros((self.nCv, self.ns), dtype=float)

        # Fluxes 
        self.fna    = np.zeros((self.nFc, 2, self.ns), dtype=float)
        self.fne    = np.zeros((self.nFc, 2), dtype=float)
        self.fni    = np.zeros((self.nFc, 2), dtype=float)
        self.fmo    = np.zeros((self.nFc, 2, self.ns), dtype=float)
        self.fhe    = np.zeros((self.nFc, 2), dtype=float)
        self.fhi    = np.zeros((self.nFc, 2), dtype=float)
        self.fhn    = np.zeros((self.nFc, 2), dtype=float)
        self.fch    = np.zeros((self.nFc, 2), dtype=float)
        self.fch_32 = np.zeros((self.nFc, 2), dtype=float)
        self.fch_52 = np.zeros((self.nFc, 2), dtype=float)
        self.fhm    = np.zeros((self.nFc, 2, self.ns), dtype=float)
        self.fkt    = np.zeros((self.nFc, 2), dtype=float)
        self.fzt    = np.zeros((self.nFc, 2), dtype=float)
        self.fch_p  = np.zeros((self.nFc, 2), dtype=float)

        # Additional fluxes
        self.fna_mdf     = np.zeros((self.nFc, 2, self.ns), dtype=float)
        self.fhe_mdf     = np.zeros((self.nFc, 2), dtype=float)
        self.fhi_mdf     = np.zeros((self.nFc, 2), dtype=float)
        self.fna_fcor    = np.zeros((self.nFc, 2, self.ns), dtype=float)
        self.fna_nodrift = np.zeros((self.nFc, 2, self.ns), dtype=float)
        self.fna_he      = np.zeros((self.nFc, 2, self.ns), dtype=float)
        self.fnaPSch     = np.zeros((self.nFc, 2, self.ns), dtype=float)
        self.fhePSch     = np.zeros((self.nFc, 2), dtype=float)
        self.fhiPSch     = np.zeros((self.nFc, 2), dtype=float)
        self.fna_eir     = np.zeros((self.nFc, 2, self.ns), dtype=float)
        self.fne_eir     = np.zeros((self.nFc, 2), dtype=float)
        self.fhe_eir     = np.zeros((self.nFc, 2), dtype=float)
        self.fhi_eir     = np.zeros((self.nFc, 2), dtype=float)
        self.fna_32      = np.zeros((self.nFc, 2, self.ns), dtype=float)
        self.fna_52      = np.zeros((self.nFc, 2, self.ns), dtype=float)
        self.fni_32      = np.zeros((self.nFc, 2), dtype=float)
        self.fni_52      = np.zeros((self.nFc, 2), dtype=float)
        self.fne_32      = np.zeros((self.nFc, 2), dtype=float)
        self.fne_52      = np.zeros((self.nFc, 2), dtype=float)
        self.fchdia      = np.zeros((self.nFc, 2), dtype=float)
        self.fchin       = np.zeros((self.nFc, 2), dtype=float)
        self.fchvispar   = np.zeros((self.nFc, 2), dtype=float)
        self.fchvisper   = np.zeros((self.nFc, 2), dtype=float)
        self.fchvisq     = np.zeros((self.nFc, 2), dtype=float)
        self.fchviskt   = np.zeros((self.nFc, 2), dtype=float)
        self.fchinert    = np.zeros((self.nFc, 2), dtype=float)
        self.fchanml     = np.zeros((self.nFc, 2), dtype=float)
        self.fht         = np.zeros((self.nFc, 2), dtype=float)
        self.fhj         = np.zeros((self.nFc, 2), dtype=float)
        self.fhp         = np.zeros((self.nFc, 2, self.ns), dtype=float)

        self.vaecrb = np.zeros((self.nFc, 2, self.ns), dtype=float)
        self.vadia  = np.zeros((self.nFc, 2, self.ns), dtype=float)
        self.wadia  = np.zeros((self.nFc, 2, self.ns), dtype=float)
        self.veecrb = np.zeros((self.nFc, 2), dtype=float)
        self.vedia  = np.zeros((self.nFc, 2), dtype=float)
        
        self.floe_noc  = np.zeros((self.nFc, 2), dtype=float)
        self.floi_noc  = np.zeros((self.nFc, 2), dtype=float)

        # Residuals
        self.resco = np.zeros((self.nCv, self.ns), dtype=float)
        self.reshe = np.zeros((self.nCv), dtype=float)
        self.reshi = np.zeros((self.nCv), dtype=float)
        self.reshn = np.zeros((self.nCv), dtype=float)
        self.resmo = np.zeros((self.nCv, self.ns), dtype=float)
        self.resmt = np.zeros((self.nCv, self.ns), dtype=float)
        self.respo = np.zeros((self.nCv, self.ns), dtype=float)
        self.reskt = np.zeros((self.nCv, self.ns), dtype=float)
        self.reszt = np.zeros((self.nCv, self.ns), dtype=float)

        # Sources
        self.sna = np.zeros((self.nCv, 2, self.ns), dtype=float)
        self.smo = np.zeros((self.nCv, 4, self.ns), dtype=float)
        self.smq = np.zeros((self.nCv, 4, self.ns), dtype=float)
        self.shi = np.zeros((self.nCv, 4), dtype=float)
        self.she = np.zeros((self.nCv, 4), dtype=float)
        self.skt = np.zeros((self.nCv, 4), dtype=float)
        self.skt_prod = np.zeros((self.nCv), dtype=float)
        self.skt_diss = np.zeros((self.nCv), dtype=float)
        
        # Coefficients
        self.calf = np.zeros((self.nFc, 2), dtype=float)
        self.cdna = np.zeros((self.nFc, 2, self.ns), dtype=float)
        self.cdpa = np.zeros((self.nFc, 2, self.ns), dtype=float)
        self.ceqp = np.zeros((self.nCv), dtype=float)
        self.chce = np.zeros((self.nFc, 2), dtype=float)
        self.chci = np.zeros((self.nFc, 2), dtype=float)
        self.chve = np.zeros((self.nFc, 2), dtype=float)
        self.chvemx = np.zeros((self.nFc), dtype=float)
        self.chvi = np.zeros((self.nFc, 2), dtype=float)
        self.chvimx = np.zeros((self.nFc), dtype=float)
        self.csig = np.zeros((self.nFc, 2), dtype=float)
        self.cvla = np.zeros((self.nFc, 2, self.ns), dtype=float)
        self.cvsa = np.zeros((self.nFc, 2, self.ns), dtype=float)
        self.cthe = np.zeros((self.nCv, self.ns), dtype=float)
        self.cthi = np.zeros((self.nCv, self.ns), dtype=float)
        self.csigin = np.zeros((self.nFc, 2, self.ns, self.ns), dtype=float)
        self.cvsa_cl = np.zeros((self.nFc, 2, self.ns), dtype=float)
        self.fllime = np.zeros((self.nFc), dtype=float)
        self.fllimi = np.zeros((self.nFc), dtype=float)
        self.fllim0fna = np.zeros((self.nFc, 2, self.ns), dtype=float)
        self.fllim0fhi = np.zeros((self.nFc, 2, self.ns), dtype=float)
        self.fllimvisc = np.zeros((self.nFc, self.ns), dtype=float)
        self.f_luc_ke = np.zeros((self.nFc), dtype=float)
        self.f_luc_ki = np.zeros((self.nFc), dtype=float)
        self.f_luc_et = np.zeros((self.nFc), dtype=float)
        self.f_luc_sg = np.zeros((self.nFc), dtype=float)
        self.f_luc_al = np.zeros((self.nFc), dtype=float)
        self.fllim_ke = np.zeros((self.nFc), dtype=float)
        self.fllim_ki = np.zeros((self.nFc), dtype=float)
        self.fllim_et = np.zeros((self.nFc), dtype=float)
        self.fllim_sg = np.zeros((self.nFc), dtype=float)
        self.fllim_al = np.zeros((self.nFc), dtype=float)
        self.fllim_al_C = np.zeros((self.nCv), dtype=float)

        self.sig0 = np.zeros((self.nCv), dtype=float)
        self.hce0 = np.zeros((self.nCv), dtype=float)
        self.alf0 = np.zeros((self.nCv), dtype=float)
        self.hci0 = np.zeros((self.nCv), dtype=float)
        self.dpa0 = np.zeros((self.nCv, self.ns), dtype=float)
        self.dna0 = np.zeros((self.nCv, self.ns), dtype=float)
        self.vsa0 = np.zeros((self.nCv, self.ns), dtype=float)
        self.vla0 = np.zeros((self.nCv, 2, self.ns), dtype=float)
        self.dkt0 = np.zeros((self.nCv), dtype=float)
        self.dna_ExB = np.zeros((self.nCv), dtype=float)
        self.hce_ExB = np.zeros((self.nCv), dtype=float)
        self.hci_ExB = np.zeros((self.nCv), dtype=float)

    def ReadB2fstatefile(self, statefilepath):
        # Description
        #------------
        # Read in the plasma state from a b2fstate file 

        # Open the file
        thisfile = open(statefilepath)

        # Get the lines
        lines = thisfile.readlines()
        si = 1

        # Set the header
        self.header = lines[0]

        # Dimensions
        dim, si = dh.ReadSOLPSField('nCv,nFc,ns', lines, si, 3, 'int', 3)
        nCv     = dim[0]
        nFc     = dim[1]
        ns      = dim[2]

        # Initialize
        self.Initialize(nCv, nFc, ns)

        self.ns     = ns # number of species
        self.nCv    = nCv # number of cells 
        self.nFc    = nFc # number of faces 
        statedim    = (nCv)
        nstatedim   = nCv 
        statedims   = (nCv, ns)
        nstatedims  = nCv*ns 
        fluxdim     = (nFc, 2)
        nfluxdim    = nFc*2
        fluxdims    = (nFc, 2, ns)
        nfluxdims   = nFc*2*ns

        # label
        si = si + 1
        self.label = lines[si]

        # Charges
        self.zamin, si  = dh.ReadSOLPSField('zamin', lines, si, ns, 'float', ns)
        self.zamax, si  = dh.ReadSOLPSField('zamax', lines, si, ns, 'float', ns)
        self.zn, si     = dh.ReadSOLPSField('zn', lines, si, ns, 'float', ns)
        self.am, si     = dh.ReadSOLPSField('am', lines, si, ns, 'float', ns)

        # State variables
        self.na, si   = dh.ReadSOLPSField('na', lines, si, 
            nstatedims, 'float', statedims)
        self.ne, si     = dh.ReadSOLPSField('ne', lines, si, 
            nstatedim, 'float', statedim)
        self.ua, si     = dh.ReadSOLPSField('ua', lines, si, 
            nstatedims, 'float', statedims)
        self.uadia, si  = dh.ReadSOLPSField('uadia', lines, si, 
            nfluxdims, 'float', fluxdims)
        self.te, si     = dh.ReadSOLPSField('te', lines, si, 
            nstatedim, 'float', statedim)
        self.ti, si     = dh.ReadSOLPSField('ti', lines, si, 
            nstatedim, 'float', statedim)
        self.tn, si     = dh.ReadSOLPSField('tn', lines, si, 
            nstatedim, 'float', statedim)
        self.po, si     = dh.ReadSOLPSField('po', lines, si, 
            nstatedim, 'float', statedim)
        self.kt, si     = dh.ReadSOLPSField('kt', lines, si, 
            nstatedim, 'float', statedim)
        self.zt, si     = dh.ReadSOLPSField('zt', lines, si, 
            nstatedim, 'float', statedim)
        

        # Fluxes 
        self.fna, si    = dh.ReadSOLPSField('fna', lines, si, 
            nfluxdims, 'float', fluxdims)
        self.fhe, si    = dh.ReadSOLPSField('fhe', lines, si, 
            nfluxdim, 'float', fluxdim)
        self.fhi, si    = dh.ReadSOLPSField('fhi', lines, si, 
            nfluxdim, 'float', fluxdim)
        self.fhn, si    = dh.ReadSOLPSField('fhn', lines, si, 
            nfluxdim, 'float', fluxdim)
        self.fch, si    = dh.ReadSOLPSField('fch', lines, si, 
            nfluxdim, 'float', fluxdim)
        self.fch_32, si = dh.ReadSOLPSField('fch_32', lines, si, 
            nfluxdim, 'float', fluxdim)
        self.fch_52, si = dh.ReadSOLPSField('fch_52', lines, si, 
            nfluxdim, 'float', fluxdim)
        self.kinrgy, si = dh.ReadSOLPSField('kinrgy', lines, si, 
            nstatedims, 'float', statedims)
        self.fhm, si    = dh.ReadSOLPSField('fhm', lines, si, 
            nfluxdims, 'float', fluxdims)
        self.fkt, si    = dh.ReadSOLPSField('fkt', lines, si, 
            nfluxdim, 'float', fluxdim)
        self.fzt, si    = dh.ReadSOLPSField('fzt', lines, si, 
            nfluxdim, 'float', fluxdim)
        self.time, si   = dh.ReadSOLPSField('time', lines, si, 
            1, 'float', 1)
        self.fch_p, si  = dh.ReadSOLPSField('fch_p', lines, si, 
            nfluxdim, 'float', fluxdim)

        # Additional fluxes
        self.fna_mdf, si     = dh.ReadSOLPSField('fna_mdf', lines, si, 
            nfluxdims, 'float', fluxdims)
        self.fhe_mdf, si     = dh.ReadSOLPSField('fhe_mdf', lines, si, 
            nfluxdim, 'float', fluxdim)
        self.fhi_mdf, si     = dh.ReadSOLPSField('fhi_mdf', lines, si, 
            nfluxdim, 'float', fluxdim)
        self.fna_fcor, si    = dh.ReadSOLPSField('fna_fcor', lines, si, 
            nfluxdims, 'float', fluxdims)
        self.fna_nodrift, si = dh.ReadSOLPSField('fna_nodrift', lines, si, 
            nfluxdims, 'float', fluxdims)
        self.fna_he, si      = dh.ReadSOLPSField('fna_he', lines, si, 
            nfluxdims, 'float', fluxdims)
        self.fnaPSch, si     = dh.ReadSOLPSField('fnaPSch', lines, si, 
            nfluxdims, 'float', fluxdims)
        self.fhePSch, si     = dh.ReadSOLPSField('fhePSch', lines, si, 
            nfluxdim, 'float', fluxdim)
        self.fhiPSch, si     = dh.ReadSOLPSField('fhiPSch', lines, si, 
            nfluxdim, 'float', fluxdim)
        self.fna_eir, si     = dh.ReadSOLPSField('fna_eir', lines, si, 
            nfluxdims, 'float', fluxdims)
        self.fne_eir, si     = dh.ReadSOLPSField('fne_eir', lines, si, 
            nfluxdim, 'float', fluxdim)
        self.fhe_eir, si     = dh.ReadSOLPSField('fhe_eir', lines, si, 
            nfluxdim, 'float', fluxdim)
        self.fhi_eir, si     = dh.ReadSOLPSField('fhi_eir', lines, si, 
            nfluxdim, 'float', fluxdim)
        self.fna_32, si      = dh.ReadSOLPSField('fna_32', lines, si, 
            nfluxdims, 'float', fluxdims)
        self.fna_52, si      = dh.ReadSOLPSField('fna_52', lines, si, 
            nfluxdims, 'float', fluxdims)
        self.fni_32, si      = dh.ReadSOLPSField('fni_32', lines, si, 
            nfluxdim, 'float', fluxdim)
        self.fni_52, si      = dh.ReadSOLPSField('fni_52', lines, si, 
            nfluxdim, 'float', fluxdim)
        self.fne_32, si      = dh.ReadSOLPSField('fne_32', lines, si, 
            nfluxdim, 'float', fluxdim)
        self.fne_52, si      = dh.ReadSOLPSField('fne_52', lines, si, 
            nfluxdim, 'float', fluxdim)
        self.fchdia, si      = dh.ReadSOLPSField('fchdia', lines, si, 
            nfluxdim, 'float', fluxdim)
        self.fchin, si       = dh.ReadSOLPSField('fchin', lines, si, 
            nfluxdim, 'float', fluxdim)
        self.fchvispar, si   = dh.ReadSOLPSField('fchvispar', lines, si, 
            nfluxdim, 'float', fluxdim)
        self.fchvisper, si   = dh.ReadSOLPSField('fchvisper', lines, si, 
            nfluxdim, 'float', fluxdim)
        self.fchvisq, si     = dh.ReadSOLPSField('fchvisq', lines, si, 
            nfluxdim, 'float', fluxdim)
        self.fchviskt, si   = dh.ReadSOLPSField('fchviskt', lines, si, 
            nfluxdim, 'float', fluxdim)
        self.fchinert, si    = dh.ReadSOLPSField('fchinert', lines, si, 
            nfluxdim, 'float', fluxdim)
        
        self.vaecrb, si = dh.ReadSOLPSField('vaecrb', lines, si, 
            nfluxdims, 'float', fluxdims)
        self.vadia, si  = dh.ReadSOLPSField('vadia', lines, si, 
            nfluxdims, 'float', fluxdims)
        self.wadia, si  = dh.ReadSOLPSField('wadia', lines, si, 
            nfluxdims, 'float', fluxdims)
        self.veecrb, si = dh.ReadSOLPSField('veecrb', lines, si, 
            nfluxdim, 'float', fluxdim)
        self.vedia, si  = dh.ReadSOLPSField('vedia', lines, si, 
            nfluxdim, 'float', fluxdim)
        
        self.floe_noc, si  = dh.ReadSOLPSField('floe_noc', lines, si, 
            nfluxdim, 'float', fluxdim)
        self.floi_noc, si  = dh.ReadSOLPSField('floi_noc', lines, si, 
            nfluxdim, 'float', fluxdim)
    
    def WriteB2fstatefile(self, statefilepath):
        # Description
        #------------
        # Write the state in a b2fstate/b2fstati fashion to a file. 
        # The header in this case will 

        # Open file
        thisfile = open(statefilepath, "w")

        # Header
        if "\n" in self.header:
            thisfile.write(self.header)
        else:
            thisfile.write(self.header + "\n")
        


        # Dimensions
        dim = np.array([self.nCv, self.nFc, self.ns])
        dh.WriteSOLPSField(thisfile, dim, 'nCv,nFc,ns', 'int')

        # Label
        if "\n" in self.label:
            thisfile.write(self.label)
        else:
            thisfile.write("*cf:    char             120    label " + "\n")
            thisfile.write(self.label.ljust(120, " ") + "\n")

        # Charges
        dh.WriteSOLPSField(thisfile, self.zamin, 'zamin', 'real')
        dh.WriteSOLPSField(thisfile, self.zamax, 'zamax', 'real')
        dh.WriteSOLPSField(thisfile, self.zn, 'zn', 'real')
        dh.WriteSOLPSField(thisfile, self.am, 'am', 'real')

        # State variables
        dh.WriteSOLPSField(thisfile, self.na, 'na', 'real')
        dh.WriteSOLPSField(thisfile, self.ne, 'ne', 'real')
        dh.WriteSOLPSField(thisfile, self.ua, 'ua', 'real')
        dh.WriteSOLPSField(thisfile, self.uadia, 'uadia', 'real')
        dh.WriteSOLPSField(thisfile, self.te, 'te', 'real')
        dh.WriteSOLPSField(thisfile, self.ti, 'ti', 'real')
        dh.WriteSOLPSField(thisfile, self.tn, 'tn', 'real')
        dh.WriteSOLPSField(thisfile, self.po, 'po', 'real')
        dh.WriteSOLPSField(thisfile, self.kt, 'kt', 'real')
        dh.WriteSOLPSField(thisfile, self.zt, 'zt', 'real')
        

        # Fluxes 
        dh.WriteSOLPSField(thisfile, self.fna, 'fna', 'real')
        dh.WriteSOLPSField(thisfile, self.fhe, 'fhe', 'real')
        dh.WriteSOLPSField(thisfile, self.fhi, 'fhi', 'real')
        dh.WriteSOLPSField(thisfile, self.fhn, 'fhn', 'real')
        dh.WriteSOLPSField(thisfile, self.fch, 'fch', 'real')
        dh.WriteSOLPSField(thisfile, self.fch_32, 'fch_32', 'real')
        dh.WriteSOLPSField(thisfile, self.fch_52, 'fch_52', 'real')
        dh.WriteSOLPSField(thisfile, self.kinrgy, 'kinrgy', 'real')
        dh.WriteSOLPSField(thisfile, self.fhm, 'fhm', 'real')
        dh.WriteSOLPSField(thisfile, self.fkt, 'fkt', 'real')
        dh.WriteSOLPSField(thisfile, self.fzt, 'fzt', 'real')
        dh.WriteSOLPSField(thisfile, self.time, 'time', 'real')
        dh.WriteSOLPSField(thisfile, self.fch_p, 'fch_p', 'real')

        # Additional fluxes
        dh.WriteSOLPSField(thisfile, self.fna_mdf, 'fna_mdf', 'real')
        dh.WriteSOLPSField(thisfile, self.fhe_mdf, 'fhe_mdf', 'real')
        dh.WriteSOLPSField(thisfile, self.fhi_mdf, 'fhi_mdf', 'real')
        dh.WriteSOLPSField(thisfile, self.fna_fcor, 'fna_fcor', 'real')
        dh.WriteSOLPSField(thisfile, self.fna_nodrift, 'fna_nodrift', 'real')
        dh.WriteSOLPSField(thisfile, self.fna_he, 'fna_he', 'real')
        dh.WriteSOLPSField(thisfile, self.fnaPSch, 'fnaPSch', 'real')
        dh.WriteSOLPSField(thisfile, self.fhePSch, 'fhePSch', 'real')
        dh.WriteSOLPSField(thisfile, self.fhiPSch, 'fhiPSch', 'real')
        dh.WriteSOLPSField(thisfile, self.fna_eir, 'fna_eir', 'real')
        dh.WriteSOLPSField(thisfile, self.fne_eir, 'fne_eir', 'real')
        dh.WriteSOLPSField(thisfile, self.fhe_eir, 'fhe_eir', 'real')
        dh.WriteSOLPSField(thisfile, self.fhi_eir, 'fhi_eir', 'real')
        dh.WriteSOLPSField(thisfile, self.fna_32, 'fna_32', 'real')
        dh.WriteSOLPSField(thisfile, self.fna_52, 'fna_52', 'real')
        dh.WriteSOLPSField(thisfile, self.fni_32, 'fni_32', 'real')
        dh.WriteSOLPSField(thisfile, self.fni_52, 'fni_52', 'real')
        dh.WriteSOLPSField(thisfile, self.fne_32, 'fne_32', 'real')
        dh.WriteSOLPSField(thisfile, self.fne_52, 'fne_52', 'real')
        dh.WriteSOLPSField(thisfile, self.fchdia, 'fchdia', 'real')
        dh.WriteSOLPSField(thisfile, self.fchin, 'fchin', 'real')
        dh.WriteSOLPSField(thisfile, self.fchvispar, 'fchvispar', 'real')
        dh.WriteSOLPSField(thisfile, self.fchvisper, 'fchvisper', 'real')
        dh.WriteSOLPSField(thisfile, self.fchvisq, 'fchvisq', 'real')
        dh.WriteSOLPSField(thisfile, self.fchviskt, 'fchviskt', 'real')
        dh.WriteSOLPSField(thisfile, self.fchinert, 'fchinert', 'real')

        dh.WriteSOLPSField(thisfile, self.vaecrb, 'vaecrb', 'real')
        dh.WriteSOLPSField(thisfile, self.vadia, 'vadia', 'real')
        dh.WriteSOLPSField(thisfile, self.wadia, 'wadia', 'real')
        dh.WriteSOLPSField(thisfile, self.veecrb, 'veecrb', 'real')
        dh.WriteSOLPSField(thisfile, self.vedia, 'vedia', 'real')

        dh.WriteSOLPSField(thisfile, self.floe_noc, 'floe_noc', 'real')
        dh.WriteSOLPSField(thisfile, self.floi_noc, 'floi_noc', 'real')

    def ReadB2fplasmfFile(self, statefilepath):
        # Description
        #------------
        # Read in a b2fplasmf file. This contains additional data such 
        # as residuals, sources, etc, but also b2fstati data. 
        # Open the file
        thisfile = open(statefilepath)

        # Get the lines
        lines = thisfile.readlines()
        si = 1

        # Set the header
        self.header = lines[0]

        # Dimensions
        dim, si = dh.ReadSOLPSField('nCv,nFc,ns', lines, si, 3, 'int', 3)
        nCv     = dim[0]
        nFc     = dim[1]
        ns      = dim[2]

        # Initialize
        self.Initialize(nCv, nFc, ns)

        self.ns     = ns # number of species
        self.nCv    = nCv # number of cells 
        self.nFc    = nFc # number of faces 
        statedim    = (nCv)
        nstatedim   = nCv 
        statedims   = (nCv, ns)
        nstatedims  = nCv*ns 
        fluxdim     = (nFc, 2)
        nfluxdim    = nFc*2
        fluxdims    = (nFc, 2, ns)
        nfluxdims   = nFc*2*ns

        # No label and charges present

        # State variables
        self.na, si   = dh.ReadSOLPSField('na', lines, si, 
            nstatedims, 'float', statedims)
        self.ne, si     = dh.ReadSOLPSField('ne', lines, si, 
            nstatedim, 'float', statedim)
        self.ua, si     = dh.ReadSOLPSField('ua', lines, si, 
            nstatedims, 'float', statedims)
        self.uadia, si  = dh.ReadSOLPSField('uadia', lines, si, 
            nfluxdims, 'float', fluxdims)
        self.te, si     = dh.ReadSOLPSField('te', lines, si, 
            nstatedim, 'float', statedim)
        self.ti, si     = dh.ReadSOLPSField('ti', lines, si, 
            nstatedim, 'float', statedim)
        self.tn, si     = dh.ReadSOLPSField('tn', lines, si, 
            nstatedim, 'float', statedim)
        self.po, si     = dh.ReadSOLPSField('po', lines, si, 
            nstatedim, 'float', statedim)
        self.kt, si     = dh.ReadSOLPSField('kt', lines, si, 
            nstatedim, 'float', statedim)
        self.zt, si     = dh.ReadSOLPSField('zt', lines, si, 
            nstatedim, 'float', statedim)
        

        # Fluxes 
        self.fna, si    = dh.ReadSOLPSField('fna', lines, si, 
            nfluxdims, 'float', fluxdims)
        self.fne, si    = dh.ReadSOLPSField('fne', lines, si, 
            nfluxdim, 'float', fluxdim)
        self.fni, si    = dh.ReadSOLPSField('fni', lines, si, 
            nfluxdim, 'float', fluxdim)
        self.fmo, si    = dh.ReadSOLPSField('fmo', lines, si, 
            nfluxdims, 'float', fluxdims)
        self.fhe, si    = dh.ReadSOLPSField('fhe', lines, si, 
            nfluxdim, 'float', fluxdim)
        self.fhi, si    = dh.ReadSOLPSField('fhi', lines, si, 
            nfluxdim, 'float', fluxdim)
        self.fhn, si    = dh.ReadSOLPSField('fhn', lines, si, 
            nfluxdim, 'float', fluxdim)
        self.fch, si    = dh.ReadSOLPSField('fch', lines, si, 
            nfluxdim, 'float', fluxdim)
        self.fch_32, si = dh.ReadSOLPSField('fch_32', lines, si, 
            nfluxdim, 'float', fluxdim)
        self.fch_52, si = dh.ReadSOLPSField('fch_52', lines, si, 
            nfluxdim, 'float', fluxdim)
        self.kinrgy, si = dh.ReadSOLPSField('kinrgy', lines, si, 
            nstatedims, 'float', statedims)
        self.fkt, si    = dh.ReadSOLPSField('fkt', lines, si, 
            nfluxdim, 'float', fluxdim)
        self.fzt, si    = dh.ReadSOLPSField('fzt', lines, si, 
            nfluxdim, 'float', fluxdim)
        self.fch_p, si  = dh.ReadSOLPSField('fch_p', lines, si, 
            nfluxdim, 'float', fluxdim)
        
        # Additional fluxes
        self.fna_mdf, si     = dh.ReadSOLPSField('fna_mdf', lines, si, 
            nfluxdims, 'float', fluxdims)
        self.fhe_mdf, si     = dh.ReadSOLPSField('fhe_mdf', lines, si, 
            nfluxdim, 'float', fluxdim)
        self.fhi_mdf, si     = dh.ReadSOLPSField('fhi_mdf', lines, si, 
            nfluxdim, 'float', fluxdim)
        self.fna_fcor, si    = dh.ReadSOLPSField('fna_fcor', lines, si, 
            nfluxdims, 'float', fluxdims)
        self.fna_nodrift, si = dh.ReadSOLPSField('fna_nodrift', lines, si, 
            nfluxdims, 'float', fluxdims)
        self.fna_he, si      = dh.ReadSOLPSField('fna_he', lines, si, 
            nfluxdims, 'float', fluxdims)
        self.fnaPSch, si     = dh.ReadSOLPSField('fnaPSch', lines, si, 
            nfluxdims, 'float', fluxdims)
        self.fhePSch, si     = dh.ReadSOLPSField('fhePSch', lines, si, 
            nfluxdim, 'float', fluxdim)
        self.fhiPSch, si     = dh.ReadSOLPSField('fhiPSch', lines, si, 
            nfluxdim, 'float', fluxdim)
        self.fna_eir, si     = dh.ReadSOLPSField('fna_eir', lines, si, 
            nfluxdims, 'float', fluxdims)
        self.fne_eir, si     = dh.ReadSOLPSField('fne_eir', lines, si, 
            nfluxdim, 'float', fluxdim)
        self.fhe_eir, si     = dh.ReadSOLPSField('fhe_eir', lines, si, 
            nfluxdim, 'float', fluxdim)
        self.fhi_eir, si     = dh.ReadSOLPSField('fhi_eir', lines, si, 
            nfluxdim, 'float', fluxdim)
        self.fna_32, si      = dh.ReadSOLPSField('fna_32', lines, si, 
            nfluxdims, 'float', fluxdims)
        self.fna_52, si      = dh.ReadSOLPSField('fna_52', lines, si, 
            nfluxdims, 'float', fluxdims)
        self.fni_32, si      = dh.ReadSOLPSField('fni_32', lines, si, 
            nfluxdim, 'float', fluxdim)
        self.fni_52, si      = dh.ReadSOLPSField('fni_52', lines, si, 
            nfluxdim, 'float', fluxdim)
        self.fne_32, si      = dh.ReadSOLPSField('fne_32', lines, si, 
            nfluxdim, 'float', fluxdim)
        self.fne_52, si      = dh.ReadSOLPSField('fne_52', lines, si, 
            nfluxdim, 'float', fluxdim)
        self.fchdia, si      = dh.ReadSOLPSField('fchdia', lines, si, 
            nfluxdim, 'float', fluxdim)
        self.fchin, si       = dh.ReadSOLPSField('fchin', lines, si, 
            nfluxdim, 'float', fluxdim)
        self.fchvispar, si   = dh.ReadSOLPSField('fchvispar', lines, si, 
            nfluxdim, 'float', fluxdim)
        self.fchvisper, si   = dh.ReadSOLPSField('fchvisper', lines, si, 
            nfluxdim, 'float', fluxdim)
        self.fchvisq, si     = dh.ReadSOLPSField('fchvisq', lines, si, 
            nfluxdim, 'float', fluxdim)
        self.fchviskt, si   = dh.ReadSOLPSField('fchviskt', lines, si, 
            nfluxdim, 'float', fluxdim)
        self.fchinert, si    = dh.ReadSOLPSField('fchinert', lines, si, 
            nfluxdim, 'float', fluxdim)
        self.fchanml, si    = dh.ReadSOLPSField('fchanml', lines, si, 
            nfluxdim, 'float', fluxdim)
        self.fht, si        = dh.ReadSOLPSField('fht', lines, si, 
            nfluxdim, 'float', fluxdim)
        self.fhj, si        = dh.ReadSOLPSField('fhj', lines, si, 
            nfluxdim, 'float', fluxdim)
        self.fhm, si        = dh.ReadSOLPSField('fhm', lines, si, 
            nfluxdims, 'float', fluxdims)
        self.fhp, si        = dh.ReadSOLPSField('fhp', lines, si, 
            nfluxdims, 'float', fluxdims)
        
        self.vaecrb, si = dh.ReadSOLPSField('vaecrb', lines, si, 
            nfluxdims, 'float', fluxdims)
        self.vadia, si  = dh.ReadSOLPSField('vadia', lines, si, 
            nfluxdims, 'float', fluxdims)
        self.wadia, si  = dh.ReadSOLPSField('wadia', lines, si, 
            nfluxdims, 'float', fluxdims)
        self.veecrb, si = dh.ReadSOLPSField('veecrb', lines, si, 
            nfluxdim, 'float', fluxdim)
        self.vedia, si  = dh.ReadSOLPSField('vedia', lines, si, 
            nfluxdim, 'float', fluxdim)
        
        self.floe_noc, si  = dh.ReadSOLPSField('floe_noc', lines, si, 
            nfluxdim, 'float', fluxdim)
        self.floi_noc, si  = dh.ReadSOLPSField('floi_noc', lines, si, 
            nfluxdim, 'float', fluxdim)
        
        # Residuals
        self.resco, si  = dh.ReadSOLPSField('resco', lines, si, 
            nstatedims, 'float', statedims)
        self.reshe, si  = dh.ReadSOLPSField('reshe', lines, si, 
            nstatedim, 'float', statedim)
        self.reshi, si  = dh.ReadSOLPSField('reshi', lines, si, 
            nstatedim, 'float', statedim)
        self.reshn, si  = dh.ReadSOLPSField('reshn', lines, si, 
            nstatedim, 'float', statedim)
        self.resmo, si  = dh.ReadSOLPSField('resmo', lines, si, 
            nstatedims, 'float', statedims)
        self.respo, si  = dh.ReadSOLPSField('respo', lines, si, 
            nstatedim, 'float', statedim)
        self.reskt, si  = dh.ReadSOLPSField('reskt', lines, si, 
            nstatedim, 'float', statedim)
        self.reszt, si  = dh.ReadSOLPSField('reszt', lines, si, 
            nstatedim, 'float', statedim)

        # Sources
        self.sna, si  = dh.ReadSOLPSField('sna', lines, si, 
            nstatedims*2, 'float', (self.nCv, 2, self.ns))
        self.smo, si  = dh.ReadSOLPSField('smo', lines, si, 
            nstatedims*4, 'float', (self.nCv, 4, self.ns))
        self.smq, si  = dh.ReadSOLPSField('smq', lines, si, 
            nstatedims*4, 'float', (self.nCv, 4, self.ns))
        self.shi, si  = dh.ReadSOLPSField('shi', lines, si, 
            nstatedim*4, 'float', (self.nCv, 4))
        self.she, si  = dh.ReadSOLPSField('she', lines, si, 
            nstatedim*4, 'float', (self.nCv, 4))
        self.skt, si  = dh.ReadSOLPSField('skt', lines, si, 
            nstatedim*4, 'float', (self.nCv, 4))
        self.skt_prod, si  = dh.ReadSOLPSField('skt_prod', lines, si, 
            nstatedim, 'float', statedim)
        self.skt_diss, si  = dh.ReadSOLPSField('skt_diss', lines, si, 
            nstatedim, 'float', statedim)
        
        # Coefficients
        self.calf, si  = dh.ReadSOLPSField('calf', lines, si, 
            nfluxdim, 'float', fluxdim)
        self.cdna, si  = dh.ReadSOLPSField('cdna', lines, si, 
            nfluxdims, 'float', fluxdims)
        self.cdpa, si  = dh.ReadSOLPSField('cdpa', lines, si, 
            nfluxdims, 'float', fluxdims)
        self.ceqp, si  = dh.ReadSOLPSField('ceqp', lines, si, 
            nstatedim, 'float', statedim)
        self.chce, si  = dh.ReadSOLPSField('chce', lines, si, 
            nfluxdim, 'float', fluxdim)
        self.chci, si  = dh.ReadSOLPSField('chci', lines, si, 
            nfluxdim, 'float', fluxdim)
        self.chve, si  = dh.ReadSOLPSField('chve', lines, si, 
            nfluxdim, 'float', fluxdim)
        self.chvemx, si  = dh.ReadSOLPSField('chvemx', lines, si, 
            nFc, 'float', (nFc)) 
        self.chvi, si  = dh.ReadSOLPSField('chvi', lines, si, 
            nfluxdim, 'float', fluxdim)
        self.chvimx, si  = dh.ReadSOLPSField('chvimx', lines, si, 
            nFc, 'float', (nFc)) 
        self.csig, si  = dh.ReadSOLPSField('csig', lines, si, 
            nfluxdim, 'float', fluxdim) 
        self.cvla, si  = dh.ReadSOLPSField('cvla', lines, si, 
            nfluxdims, 'float', fluxdims) 
        self.cvsa, si  = dh.ReadSOLPSField('cvsa', lines, si, 
            nfluxdims, 'float', fluxdims) 
        self.cthe, si  = dh.ReadSOLPSField('cthe', lines, si, 
            nstatedims, 'float', statedims) 
        self.cthi, si  = dh.ReadSOLPSField('cthi', lines, si, 
            nstatedims, 'float', statedims) 
        self.csigin, si  = dh.ReadSOLPSField('csigin', lines, si, 
            nfluxdims*ns, 'float', (nFc, 2, ns, ns)) 
        self.cvsa_cl, si  = dh.ReadSOLPSField('cvsa_cl', lines, si, 
            nfluxdims, 'float', fluxdims) 
        self.fllime, si  = dh.ReadSOLPSField('fllime', lines, si, 
            nFc, 'float', (nFc)) 
        self.fllimi, si  = dh.ReadSOLPSField('fllimi', lines, si, 
            nFc, 'float', (nFc)) 
        self.fllim0fna, si  = dh.ReadSOLPSField('fllim0fna', lines, si, 
            nfluxdims, 'float', fluxdims) 
        self.fllim0fhi, si  = dh.ReadSOLPSField('fllim0fhi', lines, si, 
            nfluxdims, 'float', fluxdims) 
        self.fllimvisc, si  = dh.ReadSOLPSField('fllimvisc', lines, si, 
            nFc*ns, 'float', [nFc, ns]) 
        self.f_luc_ke, si  = dh.ReadSOLPSField('f_luc_ke', lines, si, 
            nFc, 'float', (nFc)) 
        self.f_luc_ki, si  = dh.ReadSOLPSField('f_luc_ki', lines, si, 
            nFc, 'float', (nFc)) 
        self.f_luc_et, si  = dh.ReadSOLPSField('f_luc_et', lines, si, 
            nFc, 'float', (nFc)) 
        self.f_luc_sg, si  = dh.ReadSOLPSField('f_luc_sg', lines, si, 
            nFc, 'float', (nFc)) 
        self.f_luc_al, si  = dh.ReadSOLPSField('f_luc_al', lines, si, 
            nFc, 'float', (nFc)) 
        self.fllim_ke, si  = dh.ReadSOLPSField('fllim_ke', lines, si, 
            nFc, 'float', (nFc)) 
        self.fllim_ki, si  = dh.ReadSOLPSField('fllim_ki', lines, si, 
            nFc, 'float', (nFc)) 
        self.fllim_et, si  = dh.ReadSOLPSField('fllim_et', lines, si, 
            nFc, 'float', (nFc)) 
        self.fllim_sg, si  = dh.ReadSOLPSField('fllim_sg', lines, si, 
            nFc, 'float', (nFc)) 
        self.fllim_al, si  = dh.ReadSOLPSField('fllim_al', lines, si, 
            nFc, 'float', (nFc)) 
        self.fllim_al_c, si  = dh.ReadSOLPSField('fllim_al_c', lines, si, 
            nCv, 'float', (nCv)) 
        
        self.sig0, si  = dh.ReadSOLPSField('sig0', lines, si, 
            nCv, 'float', (nCv)) 
        self.hce0, si  = dh.ReadSOLPSField('hce0', lines, si, 
            nCv, 'float', (nCv)) 
        self.alf0, si  = dh.ReadSOLPSField('alf0', lines, si, 
            nCv, 'float', (nCv)) 
        self.hci0, si  = dh.ReadSOLPSField('hci0', lines, si, 
            nCv, 'float', (nCv)) 
        self.dpa0, si  = dh.ReadSOLPSField('dpa0', lines, si, 
            nstatedims, 'float', statedims) 
        self.dna0, si  = dh.ReadSOLPSField('dna0', lines, si, 
            nstatedims, 'float', statedims) 
        self.vsa0, si  = dh.ReadSOLPSField('vsa0', lines, si, 
            nstatedims, 'float', statedims) 
        self.vla0, si  = dh.ReadSOLPSField('vla0', lines, si, 
            nCv*2*ns, 'float', (nCv, 2, ns)) 
        self.dkt0, si  = dh.ReadSOLPSField('dkt0', lines, si, 
            nstatedim, 'float', statedim) 
        self.dna_ExB, si  = dh.ReadSOLPSField('dna_ExB', lines, si, 
            nstatedim, 'float', statedim) 
        self.hce_ExB, si  = dh.ReadSOLPSField('hce_ExB', lines, si, 
            nstatedim, 'float', statedim) 
        self.hci_ExB, si  = dh.ReadSOLPSField('hci_ExB', lines, si, 
            nstatedim, 'float', statedim) 
        
    def ComparePlasmaState(self, state2):

        # Description
        #------------
        # Run a state comparison between this plasma state and another 
        # one by comparing the following quantities for each array:
        # - max. absolute error 
        # - max. relative error
        # These values are printed for each field that is compared. 

        # Print header
        print("*=======================================================*")
        print("*                state comparison                       *")
        print("*=======================================================*")
        print("| fieldname | max abserr | max relerr | mean(field) | abserr@max relerr | fieldval@max relerr|")

        # Set fields that shouldn't be compared
        skipfields = ['header', 'label', 'nCv', 'nFc', 'ns']

        # Loop over all possible fields of the class and compare
        for i in self.__dict__:
            if not (i in skipfields):
                # Get field values
                s1 = getattr(self, i)[0]
                s2 = getattr(state2, i)[0]

                # Reshape
                if isinstance(s1, np.ndarray):
                    s1 = s1.reshape(s1.size, order='F')
                    s2 = s2.reshape(s2.size, order='F')

                    # Compute errors
                    abserr = abs(s1 - s2)
                    relerr = abs(abserr/s1) 
                    try: 
                        abserratmaxrelerr = abserr[np.nanargmax(relerr)]
                        s1valatmaxrelerr = s1[np.nanargmax(relerr)]
                    except:
                        abserratmaxrelerr = abserr[np.argmax(relerr)]
                        s1valatmaxrelerr = s1[np.argmax(relerr)]

                    # Print 
                    print("|{0:16s}| {1:.2e} | {2:.2e} | {3:.2e} | {4:.2e} | {5:.2e} |".format(i, 
                        np.max(abserr), np.nanmax(relerr), np.mean(s1), abserratmaxrelerr, 
                        s1valatmaxrelerr))
                    
    def InterpolateFromPlasmaState(self, grid, statein, gridin):
        # Description
        #------------
        # Determine the plasma state of this object (self), defined on the 
        # grid 'grid', from another plasma state 'statein' defined on grid
        # 'gridin'. Interpolation is done using the built-in scipy 
        # unstructured data interpolant and is therefore based on 
        # a triangulation of the grid cell/face values (dependent on where
        # the values to be interpolated are defined). 

        # Initialize
        #-----------
        # Reinitialize self 
        self.Initialize(grid.cell.ntot, grid.face.ntot, statein.ns) # number of species shouldn't change 

        # Set fields that should be skipped
        skipfields = ['nCv', 'nFc', 'ns']

        # Set fields that should be copied
        copyfields = ['zamin', 'zamax', 'zn', 'am', 'header', 'time']

        # Set points
        cpoints = np.array((grid.cell.x, grid.cell.y)).transpose()
        fpoints = np.array((grid.face.x, grid.face.y)).transpose()
        cpointsin = np.array((gridin.cell.x, gridin.cell.y)).transpose()
        fpointsin = np.array((gridin.face.x, gridin.face.y)).transpose()

        # Initialize interpolators (nearest neighbour)
        cinterp = interpolate.NearestNDInterpolator(cpoints, np.zeros((grid.cell.ntot), dtype=float), rescale=True)
        finterp = interpolate.NearestNDInterpolator(fpoints, np.zeros((grid.face.ntot), dtype=float), rescale=True)

        # Interpolate
        #------------
        # Loop over all possible fields of the class and compare
        starttime = time.time()
        for i in self.__dict__:
            if i in skipfields: 
                pass 
            elif i in copyfields: 
                # Copy
                setattr(self, i, getattr(statein, i))
            else: 
                # Check dimensions
                s1 = getattr(statein, i)
                s2 = getattr(self, i)
                if isinstance(s1, np.ndarray):

                    # Check dimensions
                    arrayshape = s1.shape

                    if len(arrayshape) > 4:
                        raise ValueError("InterpolateFromPlasmaState: " + 
                            "only support up to four dimensional arrays, " +
                            "higher number of dimensions not supported")
                    
                    # Check first dimension
                    doskip = False 
                    if (arrayshape[0] == statein.nCv):
                        # Cell-based interpolation
                        points = cpoints
                        pointsin = cpointsin 
                        interp = cinterp 
                    elif (arrayshape[0] == statein.nFc):
                        # Face-based interpolation
                        points = fpoints 
                        pointsin = fpointsin
                        interp = finterp 
                    else: 
                        print("InterpolateFromPlasmaState: " + 
                            "first dimension of array does not correspond to " + 
                            "cell or face dimension, skipping field: " + 
                            i)
                        doskip = True 

                    if not doskip:
                        # Check second dimension
                        if len(arrayshape) == 4:
                            for j in range(0, arrayshape[1]):
                                for k in range(0, arrayshape[2]):
                                    for m in range(0, arrayshape[3]):
                                        interp.values = np.asarray(s1[:, j, k, m])
                                        s2[:, j, k, m] = interp(points)
                                        #s2[:, j, k, m] = griddata(pointsin, s1[:, j, k, m], points, method='nearest')
                        elif len(arrayshape) == 3:
                            for j in range(0, arrayshape[1]):
                                for k in range(0, arrayshape[2]):
                                    interp.values = np.asarray(s1[:, j, k])
                                    s2[:, j, k] = interp(points)
                                    #s2[:, j, k] = griddata(pointsin, s1[:, j, k], points, method='nearest')
                        elif len(arrayshape) == 2:
                            for j in range(0, arrayshape[1]):
                                interp.values = np.asarray(s1[:, j])
                                s2[:, j] = interp(points)
                                #s2[:, j] = griddata(pointsin, s1[:, j], points, method='nearest')
                        elif len(arrayshape) == 1:
                            interp.values = np.asarray(s1)
                            s2 = interp(points)
                            #s2 = griddata(pointsin, s1, points, method='nearest')
                        else:
                            raise ValueError('InterpolateFromPlasmaState: 0D array detected, cannot interpolate')
                        
                        # Set 
                        setattr(self, i, s2)
        endtime = time.time()
        print("time elapsed for interpolation: ", endtime-starttime)
                    
