from FF_EasyPair_Cut import EasyPairCut
import json

#=========================================================================
class FF_LJ_EP_Cut(EasyPairCut):
    """
    Class for Lennard-Jones potential with cut-off.
    Inherits from EasyPairCut to handle pair interactions with a cut-off distance.
    """
    #----------------------------------------------------------------------------
    def __init__(self, **kwargs):
        """
        Initialize the FF_LJ_EP_Cut class with parameters for Lennard-Jones potential.
        """
        super().__init__(**kwargs)
        self.potential_type = 'LJ'
        self.cutoff = kwargs.get('cutoff', 10.0)  # Default cutoff distance
        self.epsilon = kwargs.get('epsilon', 1.0)  # Depth of the potential well
        self.sigma = kwargs.get('sigma', 1.0)

     #----------------------------------------------------------------------------------
    def pair_function(self, rsq, atmtype1, atmtype2):
        """Base pair function - to be overridden in child classes."""
        sigma = self.sigma[atmtype1][atmtype2]
        epsilon = self.epsilon[atmtype1][atmtype2]
        return 0.0   
    
    #----------------------------------------------------------------------------------
    def load_fpar(self, fpar):
        """
        Load force field parameters from a file or dictionary.
        """
        super().load_fpar(fpar)
        self.cutoff = fpar.get('cutoff', self.cutoff)
        self.epsilon = fpar.get('epsilon', self.epsilon)
        self.sigma = fpar.get('sigma', self.sigma)
    #----------------------------------------------------------------------------------
#==================================================================================