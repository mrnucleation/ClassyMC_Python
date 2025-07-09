"""
Sampling Input Processing
Corresponds to Script_Sampling.f90

Handles parsing input to instantiate the correct sampling type.
"""

import sys
from typing import List
from VarPrecision import dp
import CommonSampling

# Import all sampling types
from Sampling_AcceptAll import AcceptAll
from Sampling_AcceptNone import AcceptNone
from Sampling_Metropolis import Metropolis
from Sampling_TestParticle import TestParticle
from Sampling_MinMetrop import MinMetrop
from Sampling_Nested import Nested
from Sampling_Umbrella import Umbrella
from Sampling_UmbrellaWHAM import UmbrellaWHAM

def GetXCommand(line: str, position: int) -> str:
    """Simple implementation to get command at position"""
    parts = line.split()
    if position <= len(parts):
        return parts[position-1]
    return ""

def script_sampling_type(line: str) -> int:
    """
    Corresponds to Script_SamplingType
    Parse sampling type from input line and instantiate appropriate sampling method.
    
    Args:
        line: Input line containing sampling type specification
        
    Returns:
        int: Status code (0 = success, negative = error)
    """
    line_stat = 0
    
    # Get the second command which should be the sampling type
    command = GetXCommand(line, 2)
    if not command:
        return -1
    
    command = command.lower()
    
    # Instantiate the appropriate sampling type
    if command == "acceptall":
        CommonSampling.set_sampling(AcceptAll())
        
    elif command == "acceptnone":
        CommonSampling.set_sampling(AcceptNone())
        
    elif command == "metropolis":
        CommonSampling.set_sampling(Metropolis())
        
    elif command == "testparticle":
        CommonSampling.set_sampling(TestParticle())
        
    elif command == "minmetrop" or command == "min":
        CommonSampling.set_sampling(MinMetrop())
        
    elif command == "nested":
        CommonSampling.set_sampling(Nested())
        
    elif command == "umbrella":
        CommonSampling.set_sampling(Umbrella())
        
    elif command == "umbrellawham":
        CommonSampling.set_sampling(UmbrellaWHAM())
    
    else:
        print(f"Sampling Type '{command}' is not valid!", file=sys.stderr)
        return -1
    
    print(f"Sampling method set to: {command}")
    return 0

def process_sampling_input(lines: List[str]) -> int:
    """
    Process multiple lines of sampling input.
    
    Args:
        lines: List of input lines
        
    Returns:
        int: Status code (0 = success, negative = error)
    """
    for line in lines:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
            
        # Check if this is a sampling line
        first_command = GetXCommand(line, 1).lower()
        if first_command == "sampling" or first_command == "samplingtype":
            result = script_sampling_type(line)
            if result < 0:
                return result
                
        # Additional processing for sampling parameters can be added here
        elif CommonSampling.has_sampling():
            # Pass line to the sampling instance for processing
            sampling_instance = CommonSampling.get_sampling()
            if sampling_instance:
                sampling_instance.process_io(line)
    
    return 0 