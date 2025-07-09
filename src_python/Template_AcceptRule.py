"""
Acceptance Rule Template
Corresponds to Template_AcceptRule.f90

Base class for all Monte Carlo acceptance rules.
"""

import numpy as np
import sys
from abc import ABC, abstractmethod
from typing import Optional, List, Dict, Any
from .Template_Master import ClassyClass
from .VarPrecision import dp

class AcceptRule(ClassyClass):
    """
    Base class for all acceptance rules.
    Corresponds to the Fortran acceptrule type.
    """
    
    def __init__(self):
        super().__init__()
        
    def make_decision(self, trial_box, e_diff: float, disp, 
                     in_prob: Optional[float] = None, 
                     log_prob: Optional[float] = None, 
                     extra_in: Optional[float] = None) -> bool:
        """
        Corresponds to MakeDecision
        Make acceptance decision for a single box move.
        
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
        return True  # Default implementation accepts all
    
    def make_decision_2box(self, trial_box1, trial_box2, 
                          e_diff1: float, e_diff2: float,
                          disp1, disp2,
                          in_prob: Optional[float] = None,
                          log_prob: Optional[float] = None, 
                          extra_in: Optional[float] = None) -> bool:
        """
        Corresponds to MakeDecision2Box
        Make acceptance decision for a two-box move.
        
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
            bool: Whether to accept the move
        """
        raise NotImplementedError("This sampling procedure does not have a 2box acceptance rule defined")
    
    def update_statistics(self, accept: bool):
        """
        Corresponds to UpdateStatistics
        Update acceptance statistics.
        
        Args:
            accept: Whether the move was accepted
        """
        pass  # Default implementation does nothing
    
    def get_extra_terms(self, disp, trial_box) -> float:
        """
        Corresponds to GetExtraTerms
        Calculate extra terms for acceptance criterion (e.g., chemical potential).
        
        Args:
            disp: List of Perturbation objects
            trial_box: SimBox instance
            
        Returns:
            float: Extra terms to add to acceptance criterion
        """
        return 0.0  # Default implementation returns no extra terms
    
    def process_io(self, line: str) -> int:
        """
        Corresponds to ProcessIO
        Process input/output commands.
        
        Args:
            line: Input line to process
            
        Returns:
            int: Status code (0 = success, negative = error)
        """
        return 0  # Default implementation does nothing 