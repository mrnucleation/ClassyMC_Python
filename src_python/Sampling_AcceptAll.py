"""
Accept All Sampling Rule
Corresponds to Sampling_AcceptAll.f90

A sampling rule that accepts all Monte Carlo moves regardless of energy change.
Useful for system exploration and testing purposes.
"""

import numpy as np
import sys
from typing import Optional
from .VarPrecision import dp
from .Template_AcceptRule import AcceptRule

class AcceptAll(AcceptRule):
    """
    Acceptance rule that accepts every move.
    Corresponds to the Fortran AcceptAll type.
    """
    
    def __init__(self):
        super().__init__()
    
    def make_decision(self, trial_box, e_diff: float, disp, 
                     in_prob: Optional[float] = None, 
                     log_prob: Optional[float] = None, 
                     extra_in: Optional[float] = None) -> bool:
        """
        Corresponds to AcceptAll_MakeDecision
        Always accepts the move.
        
        Args:
            trial_box: SimBox instance
            e_diff: Energy difference 
            disp: List of Perturbation objects
            in_prob: Input probability (optional)
            log_prob: Log probability (optional)
            extra_in: Extra terms (optional)
            
        Returns:
            bool: Always True
        """
        return True
    
    def make_decision_2box(self, trial_box1, trial_box2, 
                          e_diff1: float, e_diff2: float,
                          disp1, disp2,
                          in_prob: Optional[float] = None,
                          log_prob: Optional[float] = None, 
                          extra_in: Optional[float] = None) -> bool:
        """
        Corresponds to AcceptAll_MakeDecision2Box
        Always accepts the two-box move.
        
        Args:
            trial_box1: First SimBox instance
            trial_box2: Second SimBox instance
            e_diff1: Energy difference for box 1
            e_diff2: Energy difference for box 2
            disp1: List of Perturbation objects for box 1
            disp2: List of Perturbation objects for box 2
            in_prob: Input probability (optional)
            log_prob: Log probability (optional)
            extra_in: Extra terms (optional)
            
        Returns:
            bool: Always True
        """
        return True 