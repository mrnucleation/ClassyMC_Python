"""
Script_TorsionType - Torsion Force Field Type Selection
Corresponds to Script_TorsionType.f90

Handles parsing input to select and instantiate the appropriate torsion force field type.
This script manages torsion force field selection for rigid, TRAPPE, harmonic, CHARMM,
and other torsional potentials used in intramolecular energy calculations.
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
# Placeholder Torsion Force Field Classes
# =============================================================================

class PlaceholderTorsionFF:
    """
    Placeholder base class for torsion force fields.
    To be replaced with actual implementations when Template_TorsionFF is converted.
    """
    
    def __init__(self, torsion_type: str):
        self.torsion_type = torsion_type
        self.phi0 = 0.0      # Reference torsion angle (radians)
        self.k = 0.0         # Force constant
        self.n = 1           # Multiplicity
        self.parameters = []  # Additional parameters
        
    def process_io(self, line: str) -> int:
        """
        Process input parameters for the torsion force field.
        
        Args:
            line: Input line with parameters
            
        Returns:
            int: Status code (0 = success, negative = error)
        """
        try:
            parts = line.split()
            if len(parts) >= 2:
                self.phi0 = float(parts[1])
            if len(parts) >= 3:
                self.k = float(parts[2])
            return 0
        except (ValueError, IndexError):
            print(f"Error processing torsion parameters: {line}", file=sys.stderr)
            return -1
    
    def compute_energy(self, phi: float) -> float:
        """
        Compute torsion energy for given dihedral angle.
        
        Args:
            phi: Dihedral angle (radians)
            
        Returns:
            float: Torsion energy
        """
        return 0.0  # Placeholder implementation
    
    def __str__(self):
        return f"{self.torsion_type}Torsion(phi0={self.phi0}, k={self.k})"

# =============================================================================
# Specific Torsion Force Field Implementations
# =============================================================================

class RigidTorsion(PlaceholderTorsionFF):
    """Rigid torsion force field - constrains dihedral angle to fixed value"""
    
    def __init__(self):
        super().__init__("rigid")
        
    def compute_energy(self, phi: float) -> float:
        """Rigid torsions have zero energy"""
        return 0.0
        
    def process_io(self, line: str) -> int:
        """Process rigid torsion parameters: torsion_type phi0"""
        try:
            parts = line.split()
            if len(parts) >= 2:
                self.phi0 = float(parts[1])  # Fixed dihedral angle value
            return 0
        except (ValueError, IndexError):
            print(f"Error processing rigid torsion parameters: {line}", file=sys.stderr)
            return -1

class HarmonicTorsion(PlaceholderTorsionFF):
    """Harmonic torsion force field - V = 0.5 * k * (phi - phi0)^2"""
    
    def __init__(self):
        super().__init__("harmonic")
        
    def compute_energy(self, phi: float) -> float:
        """Compute harmonic torsion energy"""
        return 0.5 * self.k * (phi - self.phi0)**2
        
    def process_io(self, line: str) -> int:
        """Process harmonic torsion parameters: torsion_type phi0 k"""
        try:
            parts = line.split()
            if len(parts) >= 3:
                self.phi0 = float(parts[1])  # Equilibrium angle
                self.k = float(parts[2])     # Force constant
                return 0
            else:
                print(f"Harmonic torsion requires phi0 and k parameters: {line}", file=sys.stderr)
                return -1
        except (ValueError, IndexError):
            print(f"Error processing harmonic torsion parameters: {line}", file=sys.stderr)
            return -1

class TRAPPETorsion(PlaceholderTorsionFF):
    """
    TRAPPE torsion force field - Fourier series potential
    V = k0 + k1*(1 + cos(phi)) + k2*(1 - cos(2*phi)) + k3*(1 + cos(3*phi)) + ...
    """
    
    def __init__(self):
        super().__init__("trappe")
        self.nParameters = 0
        self.k_coeffs = []  # Fourier coefficients
        
    def compute_energy(self, phi: float) -> float:
        """Compute TRAPPE torsion energy using Fourier series"""
        import math
        
        E_tors = self.k_coeffs[0] if self.k_coeffs else 0.0  # k0 term
        
        for i in range(1, len(self.k_coeffs)):
            if i % 2 == 0:  # Even terms: 1 - cos(n*phi)
                E_tors += self.k_coeffs[i] * (1.0 - math.cos(i * phi))
            else:  # Odd terms: 1 + cos(n*phi)
                E_tors += self.k_coeffs[i] * (1.0 + math.cos(i * phi))
                
        return E_tors
        
    def process_io(self, line: str) -> int:
        """
        Process TRAPPE torsion parameters: torsion_type k0 k1 k2 k3 ...
        Number of parameters determines the order of the Fourier series.
        """
        try:
            parts = line.split()
            if len(parts) < 2:
                print(f"TRAPPE torsion requires at least k0 parameter: {line}", file=sys.stderr)
                return -1
                
            # First parameter is the torsion type, rest are coefficients
            self.nParameters = len(parts) - 1
            self.k_coeffs = []
            
            for i in range(1, len(parts)):
                self.k_coeffs.append(float(parts[i]))
                
            return 0
        except (ValueError, IndexError):
            print(f"Error processing TRAPPE torsion parameters: {line}", file=sys.stderr)
            return -1
    
    def __str__(self):
        return f"TRAPPETorsion(nParams={self.nParameters}, coeffs={self.k_coeffs})"

class CHARMMTorsion(PlaceholderTorsionFF):
    """
    CHARMM torsion force field - Multiple term potential
    V = sum_i k_i * (1 + cos(n_i * phi + delta_i))
    """
    
    def __init__(self):
        super().__init__("charmm")
        self.nParameters = 0
        self.k_terms = []      # Force constants for each term
        self.n_terms = []      # Multiplicities for each term
        self.delta_terms = []  # Phase shifts for each term
        
    def compute_energy(self, phi: float) -> float:
        """Compute CHARMM torsion energy"""
        import math
        
        E_tors = 0.0
        for i in range(len(self.k_terms)):
            k = self.k_terms[i]
            n = self.n_terms[i]
            delta = self.delta_terms[i]
            E_tors += k * (1.0 + math.cos(n * phi + delta))
            
        return E_tors
        
    def process_io(self, line: str) -> int:
        """
        Process CHARMM torsion parameters: torsion_type k1 n1 delta1 k2 n2 delta2 ...
        Parameters come in triplets: force_constant multiplicity phase_shift
        """
        try:
            parts = line.split()
            if len(parts) < 4:  # At least one triplet needed
                print(f"CHARMM torsion requires at least one k,n,delta triplet: {line}", file=sys.stderr)
                return -1
                
            # Check if we have complete triplets
            n_params = len(parts) - 1  # Exclude torsion type
            if n_params % 3 != 0:
                print(f"CHARMM torsion parameters must come in triplets (k,n,delta): {line}", file=sys.stderr)
                return -1
                
            self.nParameters = n_params // 3
            self.k_terms = []
            self.n_terms = []
            self.delta_terms = []
            
            for i in range(self.nParameters):
                base_idx = 1 + i * 3
                self.k_terms.append(float(parts[base_idx]))      # Force constant
                self.n_terms.append(int(parts[base_idx + 1]))    # Multiplicity
                self.delta_terms.append(float(parts[base_idx + 2]))  # Phase shift
                
            return 0
        except (ValueError, IndexError):
            print(f"Error processing CHARMM torsion parameters: {line}", file=sys.stderr)
            return -1
    
    def __str__(self):
        return f"CHARMMTorsion(nTerms={self.nParameters}, k={self.k_terms}, n={self.n_terms})"

# =============================================================================
# Main Script Function
# =============================================================================

def script_torsion_type(line: str, torsion_num: int) -> Tuple[object, int]:
    """
    Corresponds to Script_TorsionType in Fortran.
    Parse input line to select and instantiate appropriate torsion force field.
    
    Args:
        line: Input line containing torsion type and parameters
        torsion_num: Torsion type number (for reference/debugging)
        
    Returns:
        tuple: (torsion_force_field_instance, status_code)
               status_code: 0 = success, negative = error
    """
    # Parse torsion type from first command
    torsion_type = GetXCommand(line, 1).lower().strip()
    
    if not torsion_type:
        print("ERROR: No torsion type specified", file=sys.stderr)
        return None, -1
    
    # Get registered torsion types
    torsion_registry = get_registered_torsion_types()
    
    # Check if torsion type is registered
    if torsion_type in torsion_registry:
        try:
            # Instantiate the torsion force field
            torsion_class = torsion_registry[torsion_type]
            torsion_instance = torsion_class()
            
            # Process input parameters
            result = torsion_instance.process_io(line)
            if result < 0:
                print(f"ERROR: Failed to process torsion parameters for type '{torsion_type}'", file=sys.stderr)
                return None, -1
            
            print(f"Torsion type {torsion_num}: {torsion_type} configured successfully")
            return torsion_instance, 0
            
        except Exception as e:
            print(f"ERROR: Failed to instantiate torsion type '{torsion_type}': {e}", file=sys.stderr)
            return None, -1
    else:
        print(f"ERROR: Unknown torsion type '{torsion_type}'", file=sys.stderr)
        print(f"Available types: {list(torsion_registry.keys())}", file=sys.stderr)
        return None, -1

# =============================================================================
# Registry and Utility Functions
# =============================================================================

# Global registry for torsion force field types
_TORSION_REGISTRY: Dict[str, Type] = {
    "rigid": RigidTorsion,
    "ridgid": RigidTorsion,  # Alternative spelling
    "harmonic": HarmonicTorsion,
    "trappe": TRAPPETorsion,
    "charmm": CHARMMTorsion,
}

def register_additional_torsion_types():
    """
    Register additional torsion types that may be implemented later.
    Placeholders for future development.
    """
    # Placeholder registrations for future torsion types
    future_types = {
        # "opls": OPLSTorsion,          # When implemented
        # "amber": AmberTorsion,        # When implemented  
        # "ryckaert": RyckaertTorsion,  # When implemented
        # "fourier": FourierTorsion,    # When implemented
    }
    
    # Would register these when classes are implemented
    # _TORSION_REGISTRY.update(future_types)
    pass

def list_available_torsion_types() -> Dict[str, str]:
    """
    List all available torsion types with descriptions.
    
    Returns:
        dict: Mapping of torsion type names to descriptions
    """
    descriptions = {
        "rigid": "Fixed dihedral constraint - maintains constant dihedral angle",
        "ridgid": "Alternative spelling of rigid", 
        "harmonic": "Harmonic torsion potential - V = 0.5*k*(phi-phi0)^2",
        "trappe": "TRAPPE Fourier series potential - transferable force field",
        "charmm": "CHARMM multiple term potential - V = sum k*(1+cos(n*phi+delta))",
    }
    
    # Future torsion types (commented until implemented)
    future_descriptions = {
        # "opls": "OPLS torsion potential - optimized for liquid simulations",
        # "amber": "AMBER torsion potential - biological force field",
        # "ryckaert": "Ryckaert-Bellemans potential - cosine series expansion", 
        # "fourier": "General Fourier series torsion potential",
    }
    
    return descriptions

def validate_torsion_parameters(torsion_instance) -> bool:
    """
    Validate torsion force field parameters.
    
    Args:
        torsion_instance: Torsion force field instance
        
    Returns:
        bool: True if parameters are valid
    """
    try:
        # Basic validation - can be extended
        if hasattr(torsion_instance, 'phi0'):
            if not isinstance(torsion_instance.phi0, (int, float)):
                return False
                
        if hasattr(torsion_instance, 'k'):
            if not isinstance(torsion_instance.k, (int, float)):
                return False
                
        if hasattr(torsion_instance, 'k_coeffs'):
            if not all(isinstance(k, (int, float)) for k in torsion_instance.k_coeffs):
                return False
                
        if hasattr(torsion_instance, 'n_terms'):
            if not all(isinstance(n, int) for n in torsion_instance.n_terms):
                return False
            if any(n < 0 for n in torsion_instance.n_terms):
                print("Warning: Negative multiplicity detected", file=sys.stderr)
                
        return True
        
    except Exception as e:
        print(f"Error validating torsion parameters: {e}", file=sys.stderr)
        return False

def process_torsion_input(lines: list) -> Tuple[list, int]:
    """
    Process multiple lines of torsion input.
    
    Args:
        lines: List of input lines
        
    Returns:
        tuple: (list_of_torsion_instances, status_code)
    """
    torsion_instances = []
    
    for i, line in enumerate(lines):
        line = line.strip()
        if not line or line.startswith('#'):
            continue
            
        torsion_instance, status = script_torsion_type(line, i + 1)
        if status < 0:
            return [], status
            
        if torsion_instance and validate_torsion_parameters(torsion_instance):
            torsion_instances.append(torsion_instance)
        else:
            print(f"ERROR: Invalid torsion parameters on line {i + 1}", file=sys.stderr)
            return [], -1
    
    return torsion_instances, 0

def register_torsion_type(name: str, torsion_class: Type):
    """
    Register a new torsion force field type.
    
    Args:
        name: Name of the torsion type
        torsion_class: Class implementing the torsion force field
    """
    global _TORSION_REGISTRY
    _TORSION_REGISTRY[name.lower()] = torsion_class
    print(f"Registered torsion type: {name}")

def get_registered_torsion_types() -> Dict[str, Type]:
    """
    Get all registered torsion force field types.
    
    Returns:
        dict: Mapping of type names to classes
    """
    return _TORSION_REGISTRY.copy()

# =============================================================================
# Register additional types on import
# =============================================================================
register_additional_torsion_types()

if __name__ == "__main__":
    # Test the torsion type selection
    print("Available torsion types:")
    for torsion_type, description in list_available_torsion_types().items():
        print(f"  {torsion_type}: {description}")
    
    # Test some torsion configurations
    test_lines = [
        "rigid 180.0",
        "harmonic 180.0 50.0",
        "trappe 1.0 2.0 -0.5 0.8",
        "charmm 2.0 3 0.0 1.5 2 180.0",
        "ridgid 0.0",
    ]
    
    print("\nTesting torsion configurations:")
    for i, line in enumerate(test_lines):
        print(f"\nTesting: {line}")
        torsion_instance, status = script_torsion_type(line, i + 1)
        if status == 0:
            print(f"Success: {torsion_instance}")
        else:
            print(f"Failed with status: {status}") 