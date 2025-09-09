"""
Coordinate Types
Corresponds to CoordinateTypes.f90

Defines the displacement types used by Monte Carlo moves.
"""

import numpy as np
from src_python.VarPrecision import dp


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
        self.X = newPositions
# =============================================================================
class MolDisplacement():
    """
    Base displacement class for atomic moves.
    Corresponds to the Fortran Displacement type.
    """

    def __init__(self, oldpos, moldelta, molType: int, molIndx: int, atmIndicies: np.ndarray):
        self.moldelta = moldelta
        newPositions = oldpos + moldelta
        super().__init__(molType, molIndx, atmIndicies, newPositions)
        

# =============================================================================
        


# =============================================================================
class Addition:
    """
    Displacement class for addition moves.
    Corresponds to the Fortran Addition type.
    """
    #-------------------------------------------------------------
    def __init__(self, molType: int, molIndx: int, atomTypes: np.ndarray, newPositions: np.ndarray):
        self.X = newPositions
        self.molIndx = molIndx
        self.atomTypes = atomTypes
        self.molType = molType
    #-------------------------------------------------------------
# =============================================================================


# =============================================================================
class Deletion:
    """
    Displacement class for deletion moves.
    Corresponds to the Fortran Deletion type.
    """
    
    def __init__(self, molType: int, molIndx: int, atomTypes: np.ndarray, atomIndicies: np.ndarray):
        # For deletion, we don't have new positions, so pass empty array
        self.molIndx = molIndx
        self.atomTypes = atomTypes
        self.molType = molType
        self.atomIndicies = atomIndicies
# =============================================================================


# =============================================================================
class OrthoVolChange:
    """
    Displacement class for orthogonal volume changes.
    Corresponds to the Fortran OrthoVolChange type.
    """
    
    def __init__(self, scales: np.ndarray = None, volOld: float = 0.0, volNew: float = 0.0):
        # Volume information
        self.volOld = volOld
        self.volNew = volNew
        self.dim = 3  # Default to 3 dimensions
        
        
        # do no place xScale, yScale, or zScale here
        # depreciated and will not used
        
        self.scales = scales


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
