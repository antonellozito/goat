# This module defines some convenient types for post-processing that
# mimick the types of solps (to some extent). All in unstructured 
# format. 
import numpy as np 
from . import Datahandler as dh 


#----------------------------------------------------------------------#
#                             STATE                                    #
#----------------------------------------------------------------------#

class PlasmaState:
    def __init__(self):
        # Description
        #------------
        # Default initializer

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
        self.fhe    = np.zeros((self.nFc, 2), dtype=float)
        self.fhi    = np.zeros((self.nFc, 2), dtype=float)
        self.fhn    = np.zeros((self.nFc, 2), dtype=float)
        self.fch    = np.zeros((self.nFc, 2), dtype=float)
        self.fch_32 = np.zeros((self.nFc, 2), dtype=float)
        self.fch_52 = np.zeros((self.nFc, 2), dtype=float)
        self.fhm    = np.zeros((self.nFc, 2), dtype=float)
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
        
        self.vaecrb = np.zeros((self.nFc, 2, self.ns), dtype=float)
        self.vadia  = np.zeros((self.nFc, 2, self.ns), dtype=float)
        self.wadia  = np.zeros((self.nFc, 2, self.ns), dtype=float)
        self.veecrb = np.zeros((self.nFc, 2), dtype=float)
        self.vedia  = np.zeros((self.nFc, 2), dtype=float)
        
        self.floe_noc  = np.zeros((self.nFc, 2), dtype=float)
        self.floi_noc  = np.zeros((self.nFc, 2), dtype=float)


    def Initialize(self, nCv, nFc, ns):
        # Description
        #------------
        # Actual initializer

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
        
        self.vaecrb = np.zeros((self.nFc, 2, self.ns), dtype=float)
        self.vadia  = np.zeros((self.nFc, 2, self.ns), dtype=float)
        self.wadia  = np.zeros((self.nFc, 2, self.ns), dtype=float)
        self.veecrb = np.zeros((self.nFc, 2), dtype=float)
        self.vedia  = np.zeros((self.nFc, 2), dtype=float)
        
        self.floe_noc  = np.zeros((self.nFc, 2), dtype=float)
        self.floi_noc  = np.zeros((self.nFc, 2), dtype=float)

    def ReadB2fstatefile(self, statefilepath):
        # Description
        #------------
        # Read in the plasma state from a b2fstate file 

        # Open the file
        thisfile = open(statefilepath)

        # Get the lines
        lines = thisfile.readlines()
        si = 1

        # Dimensions
        dim, si = dh.ReadSOLPSField('nCv,nFc,ns', lines, si, 3, 'int', 3)
        nCv     = dim[0]
        nFc     = dim[1]
        ns      = dim[2]

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

        # Charges
        self.zamin  = dh.ReadSOLPSField('zamin', lines, si, ns, 'float', ns)
        self.zamax  = dh.ReadSOLPSField('zamax', lines, si, ns, 'float', ns)
        self.zn     = dh.ReadSOLPSField('zn', lines, si, ns, 'float', ns)
        self.am     = dh.ReadSOLPSField('am', lines, si, ns, 'float', ns)

        # State variables
        self.na, si   = dh.ReadSOLPSField('na', lines, si, 
            nstatedims, 'float', statedims)
        self.ne, si     = dh.ReadSOLPSField('ne', lines, si, 
            nstatedim, 'float', statedim)
        self.ua, si     = dh.ReadSOLPSField('ua', lines, si, 
            nstatedims, 'float', statedims)
        self.uadia  = dh.ReadSOLPSField('uadia', lines, si, 
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
        self.kinrgy = dh.ReadSOLPSField('kinrgy', lines, si, 
            nstatedims, 'float', statedims)

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
        self.fhm, si    = dh.ReadSOLPSField('fhm', lines, si, 
            nfluxdims, 'float', fluxdims)
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
        skipfields = ['nCv', 'nFc', 'ns']

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