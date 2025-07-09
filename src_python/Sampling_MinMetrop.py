"""
Minimization Metropolis Sampling Rule
Corresponds to Sampling_MinMetrop.f90

This sampling style only accepts moves that lower the energy of the system.
It is essentially the "T->0" limit of the regular Metropolis algorithm.
Useful for minimization and finding low-energy configurations.
"""

import math
import sys
import numpy as np
from typing import Optional
from .VarPrecision import dp
from .Template_AcceptRule import AcceptRule

def IsNan(value):
    """Check if value is NaN"""
    return np.isnan(value)

class MinMetrop(AcceptRule):
    """
    Minimization Metropolis acceptance rule - only accepts downhill moves.
    Corresponds to the Fortran MinMetrop type.
    """
    
    def __init__(self):
        super().__init__()
    
    def make_decision(self, trial_box, e_diff: float, disp, 
                     in_prob: Optional[float] = None, 
                     log_prob: Optional[float] = None, 
                     extra_in: Optional[float] = None) -> bool:
        """
        Corresponds to MinMetrop_MakeDecision
        Accept only if energy decreases (downhill moves only)
        """
        # Check for NaN in energy difference
        if IsNan(e_diff):
            print("Error! Energy difference returned NaN", file=sys.stderr)
            print("MinMetrop Sampling can not process this move", file=sys.stderr)
            return False
        
        # Check probability if provided
        if in_prob is not None:
            if in_prob <= 0.0:
                return False
        elif log_prob is not None:
            if log_prob == float('-inf'):
                return False
        
        # Get extra terms
        extra_terms = extra_in if extra_in is not None else 0.0
        
        # Calculate total energy change
        total_e_diff = e_diff + extra_terms
        
        # Accept only if energy decreases (or stays the same)
        return total_e_diff <= 0.0
    
    def make_decision_2box(self, trial_box1, trial_box2, 
                          e_diff1: float, e_diff2: float,
                          disp1, disp2,
                          in_prob: Optional[float] = None,
                          log_prob: Optional[float] = None, 
                          extra_in: Optional[float] = None) -> bool:
        """
        Corresponds to MinMetrop_MakeDecision2Box
        Two-box minimization decision
        """
        # Check for NaN values
        if IsNan(e_diff1) or IsNan(e_diff2):
            print("Error! Energy difference returned NaN", file=sys.stderr)
            print("MinMetrop Sampling can not process this move", file=sys.stderr)
            return False
        
        # Check probability if provided
        if in_prob is not None:
            if in_prob <= 0.0:
                return False
        elif log_prob is not None:
            if log_prob == float('-inf'):
                return False
        
        # Get extra terms
        extra_terms = extra_in if extra_in is not None else 0.0
        
        # Calculate total energy change for both boxes
        total_e_diff = e_diff1 + e_diff2 + extra_terms
        
        # Accept only if total energy decreases (or stays the same)
        return total_e_diff <= 0.0
    
    def update_statistics(self, accept: bool):
        """
        Corresponds to MinMetrop_UpdateStatistics
        Update acceptance statistics (inherited from base class)
        """
        # For MinMetrop, just call the base class method
        super().update_statistics(accept)
    
    def get_extra_terms(self, disp, trial_box) -> float:
        """
        Corresponds to MinMetrop_GetExtraTerms
        Get additional energy terms for the move
        """
        # MinMetrop doesn't add any extra terms by default
        return 0.0
    
    def process_io(self, line: str) -> int:
        """
        Corresponds to MinMetrop_ProcessIO
        Process input commands (MinMetrop has no specific parameters)
        """
        # MinMetrop doesn't have specific input parameters
        # Return 0 for success, -1 for unrecognized command
        return -1
    
    def prologue(self):
        """
        Corresponds to MinMetrop_Prologue
        Initialize minimization metropolis sampling
        """
        print("Sampling Style: Minimization Metropolis")
        print("Only downhill (energy-lowering) moves will be accepted")
        super().prologue()
    
    def epilogue(self):
        """
        Corresponds to MinMetrop_Epilogue
        Final cleanup for minimization metropolis
        """
        print("Minimization Metropolis sampling completed")
        super().epilogue() 