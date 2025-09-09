"""
Hard Sphere Force Field
Corresponds to FF_HardSphere.f90

Implements hard sphere interactions where particles have zero interaction energy
at distances greater than their contact distance, and infinite repulsion (overlap
rejection) at shorter distances.

This is often used for:
- Reference systems in thermodynamic calculations
- Testing Monte Carlo algorithms
- Simple excluded volume effects
- Liquid structure studies

This implementation inherits from EasyPairCut for efficient energy calculations.
"""

import numpy as np
import sys
from typing import Tuple, List, Optional
from .FF_EasyPair_Cut import EasyPairCut
from .VarPrecision import dp


class HardSphere(EasyPairCut):
    """
    Hard sphere force field implementation.
    Corresponds to the Fortran HardSphere type.

    V(r) = infinity for r < sigma
    V(r) = 0 for r >= sigma
    """

    def __init__(self, nAtomTypes=1):
        super().__init__()

        # Hard sphere specific parameters
        self.sig = None          # HS diameters per type
        self.sigTable = None     # Mixed HS diameters (stored as sigma^2)

        self.nAtomTypes = nAtomTypes
    
    def constructor(self, nAtomTypes=None):
        """
        Corresponds to Constructor_HardSphere
        Allocate and initialize HS parameter arrays
        """
        if nAtomTypes is not None:
            self.nAtomTypes = nAtomTypes
        
        # Allocate per-type arrays
        self.sig = np.full(self.nAtomTypes, 1.0, dtype=dp)

        # Allocate mixing tables (stored as sigma^2 for efficiency)
        self.sigTable = np.full((self.nAtomTypes, self.nAtomTypes), 1.0, dtype=dp)

        # Initialize rMin arrays for EasyPairCut compatibility
        self.rMin = np.full(self.nAtomTypes, 0.5, dtype=dp)  # Will be updated with actual sigma values
        self.rMinTable = np.full((self.nAtomTypes, self.nAtomTypes), 0.25, dtype=dp)  # Will be updated

        # Set default cutoff to max diameter
        self.rCut = 5.0
        self.rCutSq = 25.0
        
        print(f"HardSphere Force Field initialized with {self.nAtomTypes} atom types")

    def pair_function(self, rsq: float, atmtype1: int, atmtype2: int) -> float:
        """
        Calculate hard sphere pair energy
        Corresponds to the core logic of hard sphere interactions

        Args:
            rsq: Distance squared between atoms
            atmtype1: Type of first atom (0-indexed)
            atmtype2: Type of second atom (0-indexed)

        Returns:
            float: Hard sphere energy (0 or infinity)
        """
        return 0.0

    

    
    def single_pair(self, rsq: float, atmtype1: int, atmtype2: int) -> float:
        """
        Calculate hard sphere pair energy
        Wrapper around pair_function for compatibility

        Args:
            rsq: Distance squared between atoms
            atmtype1: Type of first atom
            atmtype2: Type of second atom

        Returns:
            float: Hard sphere energy (0 or infinity)
        """
        return self.pair_function(rsq, atmtype1, atmtype2)
    
    def process_io(self, line: str) -> int:
        """
        Process hard sphere specific input parameters

        Args:
            line: Input line to process

        Expected formats:
        - type1 sigma (single type)
        - type1 type2 sigma (pair interaction)
        - rcut value (inherited from EasyPairCut)

        Returns:
            int: Status code (0 for success, negative for error)
        """
        # Handle base EasyPair commands first (like rcut)
        result = super().process_io(line)
        if result == 0:
            return result

        parts = line.strip().split()
        if not parts:
            return 0
        
        # Ensure arrays are initialized
        if self.sig is None or self.sigTable is None:
            print("ERROR: HS arrays not initialized. Call constructor() first.")
            return -1
        
        try:
            if len(parts) == 2:
                # Single type parameter: type1 sigma
                type1 = int(parts[0]) - 1  # Convert to 0-indexed
                sigma = float(parts[1])
                
                if type1 >= self.nAtomTypes or type1 < 0:
                    print(f"ERROR: Invalid atom type {type1+1}")
                    return -1
                
                # Store single-type parameter
                self.sig[type1] = sigma
                self.rMin[type1] = sigma  # For EasyPairCut compatibility

                # Update mixing tables using arithmetic mean
                for jType in range(self.nAtomTypes):
                    sig_mix = 0.5 * (sigma + self.sig[jType])
                    self.sigTable[type1, jType] = sig_mix**2  # Store as sigma^2
                    self.sigTable[jType, type1] = sig_mix**2
                    self.rMinTable[type1, jType] = sig_mix**2  # For EasyPairCut overlap detection
                    self.rMinTable[jType, type1] = sig_mix**2
                
                return 0
                
            elif len(parts) == 3:
                # Pair interaction: type1 type2 sigma
                type1 = int(parts[0]) - 1  # Convert to 0-indexed
                type2 = int(parts[1]) - 1
                sigma = float(parts[2])
                
                if (type1 >= self.nAtomTypes or type1 < 0 or 
                    type2 >= self.nAtomTypes or type2 < 0):
                    print(f"ERROR: Invalid atom types {type1+1}, {type2+1}")
                    return -1
                
                # Set pair-specific parameter
                self.sigTable[type1, type2] = sigma**2  # Store as sigma^2
                self.sigTable[type2, type1] = sigma**2
                self.rMinTable[type1, type2] = sigma**2  # For EasyPairCut overlap detection
                self.rMinTable[type2, type1] = sigma**2
                
                return 0
                
        except (ValueError, IndexError) as e:
            print(f"ERROR: Invalid HS parameter format: {line}")
            print(f"Exception: {e}")
            return -1
        
        print(f"Unknown HS command: {line}")
        return -1
    
    def get_cutoff(self) -> float:
        """
        Get the effective cutoff radius (maximum diameter)
        Uses the base class cutoff, which is set to the maximum hard sphere diameter

        Returns:
            float: Maximum hard sphere diameter
        """
        if self.sig is not None:
            return np.max(self.sig)
        return super().get_cutoff()
    
    def prologue(self):
        """Initialize hard sphere force field before simulation"""
        print(f"Hard Sphere Force Field:")
        print(f"  Number of atom types: {self.nAtomTypes}")
        
        if self.sig is not None:
            print("  HS Parameters:")
            for i in range(self.nAtomTypes):
                print(f"    Type {i+1}: sigma={self.sig[i]:.4f}")
    
    def epilogue(self):
        """Finalize hard sphere force field after simulation"""
        print("Hard Sphere Force Field simulation completed")
    
    def __str__(self):
        return f"HardSphere(nTypes={self.nAtomTypes})" 