"""
Force Field Type Script
Corresponds to Script_FieldType.f90

Handles selection and instantiation of force field types based on input commands.
This module acts as a dispatcher that creates the appropriate force field 
instance based on the force field name provided in the input.
"""

import sys
from typing import Optional, Dict, Type, Any

# Import all available force field implementations
from .Template_Forcefield import ForceField
from .FF_EasyPair_Cut import EasyPairCut
from .FF_LJ_Cut import LJ_Cut
from .FF_HardSphere import HardSphere

# Global force field instance
CurrentForceField: Optional[ForceField] = None

# Registry of available force field types
FORCEFIELD_REGISTRY: Dict[str, Type[ForceField]] = {
    'lj_cut': LJ_Cut,
    'ljcut': LJ_Cut,
    'lj/cut': LJ_Cut,
    'lennard-jones': LJ_Cut,
    'lennard_jones': LJ_Cut,
    'lj': LJ_Cut,
    
    'hardsphere': HardSphere,
    'hard_sphere': HardSphere,
    'hard-sphere': HardSphere,
    'hs': HardSphere,
    
    # Add more force fields as they are implemented
    # 'lj_shift': LJ_Shift,
    # 'lj_ele': LJ_Electrostatic,
    # 'tersoff': Tersoff,
    # 'einstein': Einstein,
}

def select_forcefield_type(line: str) -> int:
    """
    Corresponds to SelectFieldType_ProcessIO
    Select and instantiate a force field based on input line
    
    Args:
        line: Input line containing force field type and parameters
        
    Expected format:
        fieldtype <name> [nAtomTypes] [additional_params]
        
    Returns:
        int: Status code (0 for success, negative for error)
    """
    global CurrentForceField
    
    parts = line.strip().split()
    if len(parts) < 2:
        print("ERROR: fieldtype command requires at least a force field name")
        print("Usage: fieldtype <name> [nAtomTypes] [additional_params]")
        return -1
    
    field_name = parts[1].lower()
    
    # Parse optional parameters
    nAtomTypes = 1
    if len(parts) > 2:
        try:
            nAtomTypes = int(parts[2])
            if nAtomTypes <= 0:
                print(f"ERROR: Number of atom types must be positive, got {nAtomTypes}")
                return -1
        except ValueError:
            print(f"ERROR: Invalid number of atom types: {parts[2]}")
            return -1
    
    # Look up force field class
    if field_name not in FORCEFIELD_REGISTRY:
        print(f"ERROR: Unknown force field type: {field_name}")
        print("Available force field types:")
        for name in sorted(set(FORCEFIELD_REGISTRY.keys())):
            print(f"  {name}")
        return -1
    
    try:
        # Instantiate the force field
        forcefield_class = FORCEFIELD_REGISTRY[field_name]
        CurrentForceField = forcefield_class(nAtomTypes=nAtomTypes)
        
        # Initialize the force field
        CurrentForceField.constructor(nAtomTypes=nAtomTypes)
        
        print(f"Selected force field: {field_name}")
        print(f"Number of atom types: {nAtomTypes}")
        print(f"Force field instance: {CurrentForceField}")
        
        return 0
        
    except Exception as e:
        print(f"ERROR: Failed to initialize force field {field_name}")
        print(f"Exception: {e}")
        CurrentForceField = None
        return -1

def process_forcefield_line(line: str) -> int:
    """
    Process a line of force field input
    
    Args:
        line: Input line to process
        
    Returns:
        int: Status code (0 for success, negative for error)
    """
    global CurrentForceField
    
    if CurrentForceField is None:
        print("ERROR: No force field selected. Use 'fieldtype' command first.")
        return -1
    
    return CurrentForceField.process_io(line)

def get_current_forcefield() -> Optional[ForceField]:
    """
    Get the currently selected force field instance
    
    Returns:
        ForceField or None: Current force field instance
    """
    return CurrentForceField

def has_forcefield() -> bool:
    """
    Check if a force field is currently selected
    
    Returns:
        bool: True if force field is selected
    """
    return CurrentForceField is not None

def set_forcefield(forcefield: ForceField):
    """
    Set the current force field instance
    
    Args:
        forcefield: Force field instance to set as current
    """
    global CurrentForceField
    CurrentForceField = forcefield

def register_forcefield(name: str, forcefield_class: Type[ForceField]):
    """
    Register a new force field type
    
    Args:
        name: Name to register the force field under
        forcefield_class: Force field class to register
    """
    FORCEFIELD_REGISTRY[name.lower()] = forcefield_class
    print(f"Registered force field type: {name}")

def list_available_forcefields():
    """
    Print a list of all available force field types
    """
    print("Available force field types:")
    
    # Group by actual class to avoid duplicates
    class_to_names = {}
    for name, cls in FORCEFIELD_REGISTRY.items():
        if cls not in class_to_names:
            class_to_names[cls] = []
        class_to_names[cls].append(name)
    
    for cls, names in class_to_names.items():
        primary_name = names[0]
        aliases = names[1:] if len(names) > 1 else []
        
        print(f"  {primary_name} - {cls.__doc__.split('.')[0] if cls.__doc__ else cls.__name__}")
        if aliases:
            print(f"    Aliases: {', '.join(aliases)}")

def get_forcefield_cutoff() -> float:
    """
    Get the cutoff radius of the current force field
    
    Returns:
        float: Cutoff radius, or 0.0 if no force field selected
    """
    if CurrentForceField is not None:
        return CurrentForceField.get_cutoff()
    return 0.0

def initialize_forcefield():
    """
    Initialize the current force field for simulation
    """
    if CurrentForceField is not None:
        CurrentForceField.prologue()

def finalize_forcefield():
    """
    Finalize the current force field after simulation
    """
    if CurrentForceField is not None:
        CurrentForceField.epilogue()

def update_forcefield():
    """
    Update the current force field during simulation
    """
    if CurrentForceField is not None:
        CurrentForceField.update()

# Main processing function for script input
def process_script_line(line: str) -> int:
    """
    Process a script line that may contain force field commands
    
    Args:
        line: Input line from script
        
    Returns:
        int: Status code (0 for success, negative for error)
    """
    line = line.strip()
    if not line or line.startswith('#'):
        return 0  # Skip empty lines and comments
    
    parts = line.split()
    if not parts:
        return 0
    
    command = parts[0].lower()
    
    if command == 'fieldtype':
        return select_forcefield_type(line)
    elif command == 'listfields':
        list_available_forcefields()
        return 0
    elif has_forcefield():
        # Pass to current force field for processing
        return process_forcefield_line(line)
    else:
        # Unknown command and no force field to handle it
        print(f"Unknown command: {command}")
        return -1

# Script interface for external use
def script_fieldtype(line: str) -> int:
    """
    Corresponds to Script_FieldType_ProcessIO
    Main entry point for force field script processing
    
    Args:
        line: Input line to process
        
    Returns:
        int: Status code
    """
    return process_script_line(line)

if __name__ == "__main__":
    # Test the force field selection system
    print("Testing Force Field Selection System")
    print("=" * 40)
    
    # List available force fields
    list_available_forcefields()
    print()
    
    # Test LJ force field selection
    print("Testing LJ force field selection:")
    result = select_forcefield_type("fieldtype lj_cut 2")
    print(f"Result: {result}")
    print()
    
    # Test parameter setting
    if has_forcefield():
        print("Testing parameter setting:")
        result = process_forcefield_line("1 4.0 3.4 0.5")  # LJ parameters
        print(f"Result: {result}")
        print()
        
        # Test current force field info
        ff = get_current_forcefield()
        print(f"Current force field: {ff}")
        print(f"Cutoff: {get_forcefield_cutoff()}")
    
    print("Test completed.") 