"""
Coordinate Types
Corresponds to CoordinateTypes.f90

Defines the displacement types used by Monte Carlo moves.
"""

import numpy as np
from .VarPrecision import dp


# =============================================================================
class Displacement:
    """
    Base displacement class for atomic moves.
    Corresponds to the Fortran Displacement type.
    """

    def __init__(self, molType: int, molIndx: int, atmIndicies: np.ndarray, newPositions: np.ndarray):
        # Molecule information
        self.molType = molType
        self.molIndx = molIndx
        self.atmIndicies = atmIndicies
        # New positions of atoms
        self.X = newPositions
# =============================================================================
        


# =============================================================================
class Addition(Displacement):
    """
    Displacement class for addition moves.
    Corresponds to the Fortran Addition type.
    """
    
    def __init__(self):
        super().__init__()
        self.addition = True
# =============================================================================


# =============================================================================
class Deletion(Displacement):
    """
    Displacement class for deletion moves.
    Corresponds to the Fortran Deletion type.
    """
    
    def __init__(self):
        super().__init__()
        self.deletion = True
# =============================================================================


# =============================================================================
class OrthoVolChange:
    """
    Displacement class for orthogonal volume changes.
    Corresponds to the Fortran OrthoVolChange type.
    """
    
    def __init__(self):
        # Volume information
        self.volOld = 0.0
        self.volNew = 0.0
        
        # Scale factors
        self.xScale = 1.0
        self.yScale = 1.0
        self.zScale = 1.0
# =============================================================================


# =============================================================================
class TriVolChange:
    """
    Displacement class for triclinic volume changes.
    Corresponds to the Fortran TriVolChange type.
    """
    
    def __init__(self):
        # Volume information
        self.volOld = 0.0
        self.volNew = 0.0
        
        # Matrix transformation
        self.matrix = np.eye(3, dtype=dp)
# =============================================================================


# =============================================================================
class AtomExchange(Displacement):
    """
    Displacement class for atom exchange moves.
    Corresponds to the Fortran AtomExchange type.
    """
    
    def __init__(self):
        super().__init__()
        
        # Exchange information
        self.oldType = 0
        self.newType = 0
        self.oldAtmIndx = 0
        self.newAtmIndx = 0 
# =============================================================================
