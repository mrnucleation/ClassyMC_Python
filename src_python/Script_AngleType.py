"""
Script_AngleType - Angle Force Field Type Selection
Corresponds to Script_AngleType.f90

Handles parsing input to select and instantiate the appropriate angle force field type.
This script manages angle force field selection for harmonic, rigid, and other angle potentials
used in intramolecular energy calculations.
"""

import sys
from typing import Tuple, Dict, List, Type, Optional, Any
from .VarPrecision import dp

def GetXCommand(line: str, position: int) -> str:
    """Simple implementation of GetXCommand for parsing input"""
    parts = line.split()
    if position <= len(parts):
        return parts[position-1]
    return ""

# =============================================================================
# Placeholder Angle Force Field Classes
# =============================================================================

class PlaceholderAngleFF:
    """
    Placeholder base class for angle force fields.
    To be replaced with actual implementations when Template_AngleFF is converted.
    """
    
    def __init__(self, angle_type: str):
        self.angle_type = angle_type
        self.theta0 = 0.0  # Equilibrium angle (radians)
        self.k = 0.0       # Force constant
        
    def process_io(self, line: str) -> int:
        """
        Process input parameters for the angle force field.
        
        Args:
            line: Input line with parameters
            
        Returns:
            int: Status code (0 = success, negative = error)
        """
        try:
            parts = line.split()
            if len(parts) >= 2:
                self.theta0 = float(parts[1])
            if len(parts) >= 3:
                self.k = float(parts[2])
            return 0
        except (ValueError, IndexError):
            print(f"Error processing angle parameters: {line}", file=sys.stderr)
            return -1
    
    def compute_energy(self, theta: float) -> float:
        """
        Compute angle energy for given angle.
        
        Args:
            theta: Angle value (radians)
            
        Returns:
            float: Angle energy
        """
        return 0.0  # Placeholder implementation
    
    def __str__(self):
        return f"{self.angle_type}Angle(theta0={self.theta0}, k={self.k})"

# =============================================================================
# Specific Angle Force Field Implementations
# =============================================================================

class RigidAngle(PlaceholderAngleFF):
    """Rigid angle force field - constrains angle to fixed value"""
    
    def __init__(self):
        super().__init__("rigid")
        
    def compute_energy(self, theta: float) -> float:
        """Rigid angles have zero energy"""
        return 0.0
        
    def process_io(self, line: str) -> int:
        """Process rigid angle parameters"""
        try:
            parts = line.split()
            if len(parts) >= 2:
                self.theta0 = float(parts[1])  # Fixed angle value
            return 0
        except (ValueError, IndexError):
            print(f"Error processing rigid angle parameters: {line}", file=sys.stderr)
            return -1

class HarmonicAngle(PlaceholderAngleFF):
    """Harmonic angle force field - V = 0.5 * k * (theta - theta0)^2"""
    
    def __init__(self):
        super().__init__("harmonic")
        
    def compute_energy(self, theta: float) -> float:
        """Compute harmonic angle energy"""
        return 0.5 * self.k * (theta - self.theta0)**2
        
    def process_io(self, line: str) -> int:
        """Process harmonic angle parameters: angle_type theta0 k"""
        try:
            parts = line.split()
            if len(parts) >= 3:
                self.theta0 = float(parts[1])  # Equilibrium angle
                self.k = float(parts[2])       # Force constant
                return 0
            else:
                print(f"Harmonic angle requires theta0 and k parameters: {line}", file=sys.stderr)
                return -1
        except (ValueError, IndexError):
            print(f"Error processing harmonic angle parameters: {line}", file=sys.stderr)
            return -1

# =============================================================================
# Main Script Function
# =============================================================================

def script_angle_type(line: str, angle_num: int) -> Tuple[object, int]:
    """
    Corresponds to Script_AngleType in Fortran.
    Parse input line to select and instantiate appropriate angle force field.
    
    Args:
        line: Input line containing angle type and parameters
        angle_num: Angle type number (for reference/debugging)
        
    Returns:
        tuple: (angle_force_field_instance, status_code)
               status_code: 0 = success, negative = error
    """
    # Parse angle type from first command
    angle_type = GetXCommand(line, 1).lower().strip()
    
    if not angle_type:
        print("ERROR: No angle type specified", file=sys.stderr)
        return None, -1
    
    # Get registered angle types
    angle_registry = get_registered_angle_types()
    
    # Check if angle type is registered
    if angle_type in angle_registry:
        try:
            # Instantiate the angle force field
            angle_class = angle_registry[angle_type]
            angle_instance = angle_class()
            
            # Process input parameters
            result = angle_instance.process_io(line)
            if result < 0:
                print(f"ERROR: Failed to process angle parameters for type '{angle_type}'", file=sys.stderr)
                return None, -1
            
            print(f"Angle type {angle_num}: {angle_type} configured successfully")
            return angle_instance, 0
            
        except Exception as e:
            print(f"ERROR: Failed to instantiate angle type '{angle_type}': {e}", file=sys.stderr)
            return None, -1
    else:
        print(f"ERROR: Unknown angle type '{angle_type}'", file=sys.stderr)
        print(f"Available types: {list(angle_registry.keys())}", file=sys.stderr)
        return None, -1

# =============================================================================
# Registry and Utility Functions
# =============================================================================

# Global registry for angle force field types
_ANGLE_REGISTRY: Dict[str, Type] = {
    "rigid": RigidAngle,
    "ridgid": RigidAngle,  # Alternative spelling
    "harmonic": HarmonicAngle,
}

def register_additional_angle_types():
    """
    Register additional angle types that may be implemented later.
    Placeholders for future development.
    """
    # Placeholder registrations for future angle types
    future_types = {
        # "morse": MorseAngle,        # When implemented
        # "quartic": QuarticAngle,    # When implemented  
        # "cosine": CosineAngle,      # When implemented
        # "fourier": FourierAngle,    # When implemented
    }
    
    # Would register these when classes are implemented
    # _ANGLE_REGISTRY.update(future_types)
    pass

def list_available_angle_types() -> Dict[str, str]:
    """
    List all available angle types with descriptions.
    
    Returns:
        dict: Mapping of angle type names to descriptions
    """
    descriptions = {
        "rigid": "Fixed angle constraint - maintains constant angle",
        "ridgid": "Alternative spelling of rigid",
        "harmonic": "Harmonic angle potential - V = 0.5*k*(theta-theta0)^2",
    }
    
    # Future angle types (commented until implemented)
    future_descriptions = {
        # "morse": "Morse angle potential - anharmonic with realistic bond breaking",
        # "quartic": "Quartic angle potential - V = k*(theta-theta0)^4", 
        # "cosine": "Cosine angle potential - V = k*(1 + cos(n*theta + phi))",
        # "fourier": "Fourier series angle potential - sum of cosine terms",
    }
    
    return descriptions

def validate_angle_parameters(angle_instance) -> bool:
    """
    Validate angle force field parameters.
    
    Args:
        angle_instance: Angle force field instance
        
    Returns:
        bool: True if parameters are valid
    """
    try:
        # Basic validation - can be extended
        if hasattr(angle_instance, 'theta0'):
            if not isinstance(angle_instance.theta0, (int, float)):
                return False
                
        if hasattr(angle_instance, 'k'):
            if not isinstance(angle_instance.k, (int, float)):
                return False
            if angle_instance.k < 0:
                print("Warning: Negative force constant detected", file=sys.stderr)
                
        return True
        
    except Exception as e:
        print(f"Error validating angle parameters: {e}", file=sys.stderr)
        return False

def process_angle_input(lines: list) -> Tuple[list, int]:
    """
    Process multiple lines of angle input.
    
    Args:
        lines: List of input lines
        
    Returns:
        tuple: (list_of_angle_instances, status_code)
    """
    angle_instances = []
    
    for i, line in enumerate(lines):
        line = line.strip()
        if not line or line.startswith('#'):
            continue
            
        angle_instance, status = script_angle_type(line, i + 1)
        if status < 0:
            return [], status
            
        if angle_instance and validate_angle_parameters(angle_instance):
            angle_instances.append(angle_instance)
        else:
            print(f"ERROR: Invalid angle parameters on line {i + 1}", file=sys.stderr)
            return [], -1
    
    return angle_instances, 0

def register_angle_type(name: str, angle_class: Type):
    """
    Register a new angle force field type.
    
    Args:
        name: Name of the angle type
        angle_class: Class implementing the angle force field
    """
    global _ANGLE_REGISTRY
    _ANGLE_REGISTRY[name.lower()] = angle_class
    print(f"Registered angle type: {name}")

def get_registered_angle_types() -> Dict[str, Type]:
    """
    Get all registered angle force field types.
    
    Returns:
        dict: Mapping of type names to classes
    """
    return _ANGLE_REGISTRY.copy()

# =============================================================================
# Register additional types on import
# =============================================================================
register_additional_angle_types()

if __name__ == "__main__":
    # Test the angle type selection
    print("Available angle types:")
    for angle_type, description in list_available_angle_types().items():
        print(f"  {angle_type}: {description}")
    
    # Test some angle configurations
    test_lines = [
        "rigid 109.47",
        "harmonic 109.47 100.0",
        "ridgid 120.0",
    ]
    
    print("\nTesting angle configurations:")
    for i, line in enumerate(test_lines):
        print(f"\nTesting: {line}")
        angle_instance, status = script_angle_type(line, i + 1)
        if status == 0:
            print(f"Success: {angle_instance}")
        else:
            print(f"Failed with status: {status}") 