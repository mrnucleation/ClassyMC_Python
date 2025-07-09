"""
Test Particle Sampling Rule
Corresponds to Sampling_TestParticle.f90

Sampling method for use in Widom Test Particle simulations.
"""

import math
import sys
import numpy as np
from typing import Optional
from .VarPrecision import dp
from .Template_AcceptRule import AcceptRule

def grnd():
    """Generate random number between 0 and 1"""
    return np.random.random()

class TestParticle(AcceptRule):
    """
    Test particle sampling for Widom insertion method.
    Corresponds to the Fortran TestParticle type.
    """
    
    def __init__(self):
        super().__init__()
        
        # Class attributes that correspond to Fortran variables
        self.use_diff = False
        self.e_diff = 0.0
        self.ensemb_avg = 0.0
        self.cnt = 0.0
    
    def make_decision(self, trial_box, e_diff: float, disp, 
                     in_prob: Optional[float] = None, 
                     log_prob: Optional[float] = None, 
                     extra_in: Optional[float] = None) -> bool:
        """
        Corresponds to TestParticle_MakeDecision
        Make acceptance decision for test particle insertion.
        
        Args:
            trial_box: SimBox instance
            e_diff: Energy difference 
            disp: List of Perturbation objects
            in_prob: Input probability (optional)
            log_prob: Log probability (optional)
            extra_in: Extra terms (optional)
            
        Returns:
            bool: Whether to accept the move
        """
        accept = False
        
        if in_prob is not None and in_prob <= 0.0:
            return False
        
        # Check displacement type and handle accordingly
        # Note: This is a simplified version - full implementation would need
        # proper coordinate type checking
        
        # The purpose of this section is to add any terms such as the isobaric or
        # grand canonical ensemble terms (IE the PV or chemical potential) to the
        # detailed balance condition.
        extra_terms = 0.0
        
        # For Addition moves in test particle method
        if hasattr(disp[0], 'move_type') and disp[0].move_type == 'Addition':
            self.use_diff = True
            self.e_diff = e_diff
            return False
        
        # For Deletion moves  
        elif hasattr(disp[0], 'move_type') and disp[0].move_type == 'Deletion':
            raise NotImplementedError("TestParticle Sampling Can't handle deletion just yet")
        
        # For Volume changes
        elif hasattr(disp[0], 'move_type') and disp[0].move_type == 'VolChange':
            if hasattr(disp[0], 'vol_new') and hasattr(disp[0], 'vol_old'):
                extra_terms += ((disp[0].vol_new - disp[0].vol_old) * 
                               trial_box.pressure * trial_box.beta)
        
        # Calculate acceptance probability
        if in_prob is not None:
            prob_term = math.log(in_prob)
        elif log_prob is not None:
            prob_term = log_prob
        else:
            raise ValueError("Probability has not been passed into Sampling")
        
        bias_e = -trial_box.beta * e_diff + prob_term + extra_terms
        
        if bias_e > 0.0:
            accept = True
        elif bias_e > math.log(grnd()):
            accept = True
        
        return accept
    
    def update_statistics(self, accept: bool):
        """
        Corresponds to TestParticle_UpdateStatistics
        Update test particle statistics.
        
        Args:
            accept: Whether the move was accepted
        """
        if self.use_diff:
            self.use_diff = False
            # Could add additional statistics handling here
        
        # Update ensemble averages if needed
        # This would typically involve accumulating chemical potential estimates 