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
"""

import numpy as np
import sys
from typing import Tuple, List, Optional
from .Template_Forcefield import ForceField
from .VarPrecision import dp


class HardSphere(ForceField):
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
        
        # Set default cutoff to max diameter
        self.rCut = 5.0
        self.rCutSq = 25.0
        
        print(f"HardSphere Force Field initialized with {self.nAtomTypes} atom types")
    
    def detailed_calc(self, curbox) -> Tuple[float, bool]:
        """
        Corresponds to Detailed_HardSphere
        Check for hard sphere overlaps - returns infinite energy if overlapping
        
        Args:
            curbox: SimBox instance
            
        Returns:
            tuple: (total_energy, accept_flag)
        """
        atoms = curbox.get_coordinates()
        
        curbox.ETable.fill(0.0)
        
        # Check all pairs for overlaps
        for iAtom in range(curbox.nMaxAtoms - 1):
            if not curbox.is_active(iAtom):
                continue
                
            atmType1 = curbox.AtomType[iAtom]
            
            for jAtom in range(iAtom + 1, curbox.nMaxAtoms):
                if not curbox.is_active(jAtom):
                    continue
                    
                # Skip intramolecular interactions
                if curbox.MolIndx[jAtom] == curbox.MolIndx[iAtom]:
                    continue
                
                # Calculate distance
                rx = atoms[iAtom, 0] - atoms[jAtom, 0]
                ry = atoms[iAtom, 1] - atoms[jAtom, 1]
                rz = atoms[iAtom, 2] - atoms[jAtom, 2]
                
                # Apply periodic boundary conditions
                rx, ry, rz = curbox.boundary(rx, ry, rz)
                rsq = rx*rx + ry*ry + rz*rz
                
                # Check for overlap
                atmType2 = curbox.AtomType[jAtom]
                sig_sq = self.sigTable[atmType1, atmType2] if self.sigTable is not None else 1.0
                
                if rsq < sig_sq:
                    print(f"Hard sphere overlap detected!")
                    print(f"Distance: {np.sqrt(rsq):.6f}")
                    print(f"Required distance: {np.sqrt(sig_sq):.6f}")
                    print(f"Atoms: {iAtom}, {jAtom}")
                    print(f"Types: {atmType1+1}, {atmType2+1}")
                    return float('inf'), False
        
        print("Total Hard Sphere Energy: 0.0 (no overlaps)")
        return 0.0, True
    
    def diff_calc(self, curbox, disp, tempList=None, tempNNei=None) -> Tuple[float, bool]:
        """
        Corresponds to DiffECalc_HardSphere
        Check if perturbation causes hard sphere overlaps
        
        Args:
            curbox: SimBox instance
            disp: List of Perturbation objects
            tempList: Temporary neighbor list
            tempNNei: Temporary neighbor counts
            
        Returns:
            tuple: (energy_difference, accept_flag)
        """
        accept = True
        curbox.dETable.fill(0.0)
        
        # Determine perturbation type and dispatch accordingly
        if hasattr(disp[0], 'atmIndx'):
            # Displacement move
            accept = self.shift_calc_single(curbox, disp, tempList, tempNNei)
            
        elif hasattr(disp[0], 'molType') and hasattr(disp[0], 'addition'):
            # Addition move
            accept = self.new_calc(curbox, disp, tempList, tempNNei)
            
        elif hasattr(disp[0], 'molType') and hasattr(disp[0], 'deletion'):
            # Deletion move - always acceptable for hard spheres
            accept = True
            
        elif hasattr(disp[0], 'volNew'):
            # Volume change move
            accept = self.ortho_vol_calc(curbox, disp)
            
        elif hasattr(disp[0], 'newType') and hasattr(disp[0], 'oldType'):
            # Atom exchange move
            accept = self.atom_exchange(curbox, disp)
            
        else:
            print("Unknown Perturbation Type Encountered by HardSphere", file=sys.stderr)
            return 0.0, False
        
        # Hard spheres: energy is always 0 if no overlaps, infinite if overlaps
        return 0.0, accept
    
    def shift_calc_single(self, curbox, disp, tempList=None, tempNNei=None) -> bool:
        """
        Check if displacement moves cause hard sphere overlaps
        
        Args:
            curbox: SimBox instance
            disp: Displacement perturbations
            tempList: Optional neighbor list
            tempNNei: Optional neighbor counts
            
        Returns:
            bool: Accept flag (True if no overlaps)
        """
        atoms = curbox.get_coordinates()
        
        # Process each displaced atom
        for displacement in disp:
            if not hasattr(displacement, 'atmIndx'):
                continue
                
            atmIndx = displacement.atmIndx
            atmType1 = curbox.AtomType[atmIndx]
            
            # Check new position against all other atoms
            for jAtom in range(curbox.nMaxAtoms):
                if not curbox.is_active(jAtom) or jAtom == atmIndx:
                    continue
                if curbox.MolIndx[jAtom] == curbox.MolIndx[atmIndx]:
                    continue
                
                rx = displacement.x_new - atoms[jAtom, 0]
                ry = displacement.y_new - atoms[jAtom, 1]
                rz = displacement.z_new - atoms[jAtom, 2]
                rx, ry, rz = curbox.boundary(rx, ry, rz)
                rsq = rx*rx + ry*ry + rz*rz
                
                atmType2 = curbox.AtomType[jAtom]
                sig_sq = self.sigTable[atmType1, atmType2] if self.sigTable is not None else 1.0
                
                if rsq < sig_sq:
                    return False  # Overlap detected
        
        return True  # No overlaps
    
    def new_calc(self, curbox, disp, tempList=None, tempNNei=None) -> bool:
        """
        Check if addition moves cause hard sphere overlaps
        
        Args:
            curbox: SimBox instance
            disp: Addition perturbations
            tempList: Temporary neighbor list
            tempNNei: Temporary neighbor counts
            
        Returns:
            bool: Accept flag (True if no overlaps)
        """
        # For hard spheres, just check if new atoms overlap with existing ones
        atoms = curbox.get_coordinates()
        
        for addition in disp:
            if not hasattr(addition, 'molType'):
                continue
                
            # Check each atom in the new molecule
            for new_atom in addition.atoms:
                atmType1 = new_atom.atmType
                
                for jAtom in range(curbox.nMaxAtoms):
                    if not curbox.is_active(jAtom):
                        continue
                    
                    rx = new_atom.x - atoms[0, jAtom]
                    ry = new_atom.y - atoms[1, jAtom]
                    rz = new_atom.z - atoms[2, jAtom]
                    rx, ry, rz = curbox.boundary(rx, ry, rz)
                    rsq = rx*rx + ry*ry + rz*rz
                    
                    atmType2 = curbox.AtomType[jAtom]
                    sig_sq = self.sigTable[atmType1, atmType2] if self.sigTable is not None else 1.0
                    
                    if rsq < sig_sq:
                        return False  # Overlap detected
        
        return True  # No overlaps
    
    def ortho_vol_calc(self, curbox, disp) -> bool:
        """
        Check if volume change causes hard sphere overlaps
        
        Args:
            curbox: SimBox instance
            disp: Volume change perturbations
            
        Returns:
            bool: Accept flag (True if no overlaps)
        """
        # For volume changes, need to check if scaling causes overlaps
        vol_ratio = disp[0].volNew / curbox.volume
        scale_factor = vol_ratio**(1.0/3.0)
        
        atoms = curbox.get_coordinates()
        
        # Check all pairs with scaled positions
        for iAtom in range(curbox.nMaxAtoms - 1):
            if not curbox.is_active(iAtom):
                continue
                
            atmType1 = curbox.AtomType[iAtom]
            
            for jAtom in range(iAtom + 1, curbox.nMaxAtoms):
                if not curbox.is_active(jAtom):
                    continue
                    
                if curbox.MolIndx[jAtom] == curbox.MolIndx[iAtom]:
                    continue
                
                # Calculate scaled distance
                rx = (atoms[iAtom, 0] - atoms[jAtom, 0]) * scale_factor
                ry = (atoms[iAtom, 1] - atoms[jAtom, 1]) * scale_factor
                rz = (atoms[iAtom, 2] - atoms[jAtom, 2]) * scale_factor
                rx, ry, rz = curbox.boundary(rx, ry, rz)  # Apply new boundary conditions
                rsq = rx*rx + ry*ry + rz*rz
                
                atmType2 = curbox.AtomType[jAtom]
                sig_sq = self.sigTable[atmType1, atmType2] if self.sigTable is not None else 1.0
                
                if rsq < sig_sq:
                    return False  # Overlap detected
        
        return True  # No overlaps
    
    def atom_exchange(self, curbox, disp) -> bool:
        """
        Check if atom exchange causes hard sphere overlaps
        
        Args:
            curbox: SimBox instance
            disp: Exchange perturbations
            
        Returns:
            bool: Accept flag (True if no overlaps)
        """
        atoms = curbox.get_coordinates()
        
        for exchange in disp:
            if not hasattr(exchange, 'atmIndx'):
                continue
                
            atmIndx = exchange.atmIndx
            newType = exchange.newType
            
            # Check new atom type against all neighbors
            for jAtom in range(curbox.nMaxAtoms):
                if not curbox.is_active(jAtom) or jAtom == atmIndx:
                    continue
                if curbox.MolIndx[jAtom] == curbox.MolIndx[atmIndx]:
                    continue
                
                rx = atoms[atmIndx, 0] - atoms[jAtom, 0]
                ry = atoms[atmIndx, 1] - atoms[jAtom, 1]
                rz = atoms[atmIndx, 2] - atoms[jAtom, 2]
                rx, ry, rz = curbox.boundary(rx, ry, rz)
                rsq = rx*rx + ry*ry + rz*rz
                
                atmType2 = curbox.AtomType[jAtom]
                sig_sq = self.sigTable[newType, atmType2] if self.sigTable is not None else 1.0
                
                if rsq < sig_sq:
                    return False  # Overlap detected
        
        return True  # No overlaps
    
    def single_pair(self, rsq: float, atmtype1: int, atmtype2: int) -> float:
        """
        Calculate hard sphere pair energy
        
        Args:
            rsq: Distance squared between atoms
            atmtype1: Type of first atom
            atmtype2: Type of second atom
            
        Returns:
            float: Hard sphere energy (0 or infinity)
        """
        if self.sigTable is not None:
            sig_sq = self.sigTable[atmtype1, atmtype2]
        else:
            sig_sq = 1.0
            
        return float('inf') if rsq < sig_sq else 0.0
    
    def process_io(self, line: str) -> int:
        """
        Process hard sphere specific input parameters
        
        Args:
            line: Input line to process
            
        Expected formats:
        - type1 sigma (single type)
        - type1 type2 sigma (pair interaction)
        
        Returns:
            int: Status code (0 for success, negative for error)
        """
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
                
                # Update mixing tables using arithmetic mean
                for jType in range(self.nAtomTypes):
                    sig_mix = 0.5 * (sigma + self.sig[jType])
                    self.sigTable[type1, jType] = sig_mix**2  # Store as sigma^2
                    self.sigTable[jType, type1] = sig_mix**2
                
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
        
        Returns:
            float: Maximum hard sphere diameter
        """
        if self.sig is not None:
            return np.max(self.sig)
        return self.rCut
    
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