"""
Accept None Sampling Rule
Corresponds to Sampling_AcceptNone.f90

As its name suggests, this sampling style rejects every move. Not sure why
you will need it, but hey it's here if you do.
"""

from .Template_AcceptRule import AcceptRule
from .VarPrecision import dp
from typing import Optional

class AcceptNone(AcceptRule):
    """
    Acceptance rule that rejects every move.
    Corresponds to the Fortran AcceptNone type.
    """
    
    def __init__(self):
        super().__init__()
    
    def make_decision(self, trial_box, e_diff: float, disp, 
                     in_prob: Optional[float] = None, 
                     log_prob: Optional[float] = None, 
                     extra_in: Optional[float] = None) -> bool:
        """
        Corresponds to AcceptNone_MakeDecision
        Always rejects the move.
        
        Args:
            trial_box: SimBox instance
            e_diff: Energy difference 
            disp: List of Perturbation objects
            in_prob: Input probability (optional)
            log_prob: Log probability (optional)
            extra_in: Extra terms (optional)
            
        Returns:
            bool: Always False
        """
        return False
    
    def make_decision_2box(self, trial_box1, trial_box2, 
                          e_diff1: float, e_diff2: float,
                          disp1, disp2,
                          in_prob: Optional[float] = None,
                          log_prob: Optional[float] = None, 
                          extra_in: Optional[float] = None) -> bool:
        """
        Corresponds to AcceptNone_MakeDecision2Box
        Always rejects the two-box move.
        
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
            bool: Always False
        """
        return False 