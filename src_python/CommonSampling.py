"""
Common Sampling Module
Corresponds to CommonSampling.f90

Manages the global sampling instance for the simulation.
"""

from .VarPrecision import dp
from .Template_AcceptRule import AcceptRule
from typing import Optional

# Global sampling instance
sampling: Optional[AcceptRule] = None

def set_sampling(sampling_instance: AcceptRule):
    """Set the global sampling instance"""
    global sampling
    sampling = sampling_instance

def get_sampling() -> Optional[AcceptRule]:
    """Get the global sampling instance"""
    return sampling

def has_sampling() -> bool:
    """Check if sampling instance is set"""
    return sampling is not None 