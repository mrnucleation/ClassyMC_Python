import math
import sys
import numpy as np
from src_python.VarPrecision import dp
from src_python.Template_AcceptRule import AcceptRule

# Simple implementations for missing functions
def IsNan(value):
    """Check if value is NaN"""
    return np.isnan(value)

np.random.random()

class Metropolis(AcceptRule):
    """Basic Metropolis Rule for molecular systems."""
    
    def __init__(self):
        super().__init__()
    
    def make_decision(self, trial_box, e_diff, disp, in_prob=None, log_prob=None, extra_in=None):
        """
        Make acceptance decision based on Metropolis criterion.
        
        Args:
            trial_box: SimBox instance
            e_diff: Energy difference (float)
            disp: List of Perturbation objects
            in_prob: Input probability (optional)
            log_prob: Log probability (optional)
            extra_in: Extra terms (optional)
            
        Returns:
            bool: Whether to accept the move
        """
        if IsNan(e_diff):
            print("ERROR! Invalid energy has been passed to the sampling routine!", file=sys.stderr)
            print(e_diff, file=sys.stderr)
            raise RuntimeError("Invalid energy in sampling routine")
        
        accept = False
        
        if in_prob is not None:
            if in_prob <= 0.0:
                return False
        
        if extra_in is not None:
            extra_terms = extra_in
        else:
            extra_terms = 0.0
        
        if in_prob is not None:
            prob_term = math.log(in_prob)
        elif log_prob is not None:
            prob_term = log_prob
        else:
            print("Coding Error! Probability has not been passed into Sampling", file=sys.stderr)
            raise RuntimeError("No probability provided")
        
        bias_e = -trial_box.beta * e_diff + prob_term + extra_terms
        
        if bias_e >= 0.0:
            accept = True
        elif bias_e > math.log(np.random.random()):
            accept = True
        
        return accept
    # ---------------------------------------------------------------------------
    def make_decision_2box(self, trial_box1, trial_box2, e_diff1, e_diff2, 
                          disp1, disp2, in_prob=None, log_prob=None, extra_in=None):
        """
        Make acceptance decision for two-box system.
        
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
        accept = False
        
        if in_prob is not None:
            if in_prob <= 0.0:
                return False
        
        if extra_in is not None:
            extra_terms = extra_in
        else:
            extra_terms = 0.0
        
        if in_prob is not None:
            prob_term = math.log(in_prob)
        elif log_prob is not None:
            prob_term = log_prob
        else:
            print("Coding Error! Probability has not been passed into Sampling", file=sys.stderr)
            raise RuntimeError("No probability provided")
        
        bias_e = (-trial_box1.beta * e_diff1 - 
                  trial_box2.beta * e_diff2 + 
                  prob_term + extra_terms)
        
        if bias_e > 0.0:
            accept = True
        elif bias_e > math.log(np.random.random()):
            accept = True
        
        return accept
