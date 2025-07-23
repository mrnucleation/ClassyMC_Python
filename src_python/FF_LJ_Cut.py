"""
Lennard-Jones 12-6 Force Field with Cutoff
Corresponds to FF_EP_LJ_Cut.f90

Implements the standard 12-6 Lennard-Jones potential with cutoff and optional tail corrections.
This is one of the most commonly used force fields in molecular simulation.

The potential is given by:
V(r) = 4*epsilon * [(sigma/r)^12 - (sigma/r)^6]

This implementation inherits from EasyPairCut and provides:
- Standard Lorentz-Berthelot mixing rules
- Optional tail corrections for long-range interactions  
- Support for multiple atom types
- Efficient distance-based cutoffs

Usage:
    lj_ff = LJ_Cut(nAtomTypes=2)
    lj_ff.constructor()
    lj_ff.process_io("1 1.0 1.0 0.5")  # type1 epsilon sigma rmin
    lj_ff.process_io("rcut 5.0")       # cutoff distance
"""

import numpy as np
import sys
from typing import Optional, Tuple, List
from .FF_EasyPair_Cut import EasyPairCut
from .VarPrecision import dp


class LJ_Cut(EasyPairCut):
    """
    Lennard-Jones 12-6 potential with cutoff.
    Corresponds to the Fortran EP_LJ_Cut type.
    
    V(r) = 4*epsilon * [(sigma/r)^12 - (sigma/r)^6]
    """
    
    #-------------------------------------------------------------------------
    def __init__(self, nAtomTypes=1):
        super().__init__()
        
        # LJ-specific parameters (will be allocated in constructor)
        self.epsilon = None          # LJ epsilon parameters per type
        self.sigma = None          # LJ sigma parameters per type
        
        # Mixing rule tables
        self.epsTable = None     # Mixed epsilon values
        self.sigTable = None     # Mixed sigma values (stored as sigma^2)
        
        # Type counting for tail corrections
        self.PerMolTypeCount = None
        self.TypeCount = None
        
        self.nAtomTypes = nAtomTypes
        
        # Allocate per-type arrays
        self.epsilon = np.full(self.nAtomTypes, 1.0, dtype=dp)  # Default epsilon = 1.0 (not 4.0)
        self.sigma = np.full(self.nAtomTypes, 1.0, dtype=dp)
        self.rMin = np.full(self.nAtomTypes, 0.5, dtype=dp)
        
        # Allocate mixing tables
        self.epsTable = np.full((self.nAtomTypes, self.nAtomTypes), 4.0, dtype=dp)  # 4*epsilon here
        self.sigTable = np.full((self.nAtomTypes, self.nAtomTypes), 1.0, dtype=dp)
        self.rMinTable = np.full((self.nAtomTypes, self.nAtomTypes), 0.25, dtype=dp)  # Stored as rMin^2
        
        # Set default cutoff
        self.rCut = 5.0
        self.rCutSq = 25.0
        
        print(f"LJ/Cut Force Field initialized with {self.nAtomTypes} atom types")
    
    #-------------------------------------------------------------------------
    def pair_function(self, rsq: float, atmtype1: int, atmtype2: np.ndarray) -> float:
        """
        Corresponds to PairFunction_EP_LJ_Cut
        Calculate Lennard-Jones 12-6 energy
        
        Args:
            rsq: Distance squared between atoms
            atmtype1: Type of first atom (0-indexed)
            atmtype2: Type of second atom (0-indexed)
            
        Returns:
            float: LJ pair energy
        """
        if self.epsTable is None or self.sigTable is None:
            print("ERROR: LJ parameters not initialized")
            return 0.0
            
        ep = self.epsTable[atmtype1, atmtype2]
        sig_sq = self.sigTable[atmtype1, atmtype2]

        print(f"sig_sq: {sig_sq}, rsq: {rsq}")

        # Calculate (sigma/r)^2
        inv_rsq = sig_sq / rsq
        
        # Calculate (sigma/r)^6
        inv_r6 = inv_rsq * inv_rsq * inv_rsq
        
        # Calculate (sigma/r)^12
        inv_r12 = inv_r6 * inv_r6
        
        # LJ potential: ep * [(sigma/r)^12 - (sigma/r)^6]
        # Note: ep already contains the factor of 4*epsilon from mixing rules
        E_Pair = ep * (inv_r12 - inv_r6)
        
        return E_Pair
    
    #-------------------------------------------------------------------------
    def tail_correction(self, curbox, disp=None) -> float:
        """
        Corresponds to TailCorrection_EP_LJ_Cut
        Calculate tail corrections for long-range LJ interactions
        
        Args:
            curbox: SimBox instance
            disp: Optional perturbations for delta corrections
            
        Returns:
            float: Tail correction energy
        """
        # Ensure type counting arrays are initialized
        if self.PerMolTypeCount is None:
            self.compute_type_count(curbox)
        
        E_Corr = 0.0
        
        if disp is None:
            # Total tail correction calculation
            volume = curbox.volume
            NMol = curbox.NMol[:curbox.nMolTypes] if hasattr(curbox, 'nMolTypes') else [curbox.nMolTotal]
            self.count_type(NMol)
            E_Corr = self.sum_tail()
            E_Corr = E_Corr / volume
            
        else:
            # Delta tail correction for moves
            volume = curbox.volume
            
            # Determine move type and calculate delta
            if hasattr(disp[0], 'atmIndx'):
                # Displacement moves don't change density
                return 0.0
                
            elif hasattr(disp[0], 'molType') and hasattr(disp[0], 'addition'):
                # Addition move
                molType = disp[0].molType
                NMol = curbox.NMol[:curbox.nMolTypes].copy()
                
                # Calculate old tail correction
                self.count_type(NMol)
                E_old = self.sum_tail() / volume
                
                # Calculate new tail correction
                NMol[molType] += 1
                self.count_type(NMol)
                E_new = self.sum_tail() / volume
                
                E_Corr = E_new - E_old
                
            elif hasattr(disp[0], 'molType') and hasattr(disp[0], 'deletion'):
                # Deletion move
                molType = disp[0].molType
                NMol = curbox.NMol[:curbox.nMolTypes].copy()
                
                # Calculate old tail correction
                self.count_type(NMol)
                E_old = self.sum_tail() / volume
                
                # Calculate new tail correction
                NMol[molType] -= 1
                self.count_type(NMol)
                E_new = self.sum_tail() / volume
                
                E_Corr = E_new - E_old
                
            elif hasattr(disp[0], 'volNew'):
                # Volume change move
                vol_new = disp[0].volNew
                NMol = curbox.NMol[:curbox.nMolTypes] if hasattr(curbox, 'nMolTypes') else [curbox.nMolTotal]
                self.count_type(NMol)
                tail_energy = self.sum_tail()
                E_Corr = tail_energy * (1.0/vol_new - 1.0/volume)
                
            else:
                print("LJ Tail Corrections not supported for this move type", file=sys.stderr)
                return 0.0
        
        return E_Corr
    
    #-------------------------------------------------------------------------
    def sum_tail(self) -> float:
        """
        Corresponds to SumTail_EP_LJ_Cut
        Sum tail corrections over all atom type pairs
        
        Returns:
            float: Total tail correction energy
        """
        pi = np.pi
        E_Corr = 0.0
        
        if self.epsTable is None or self.sigTable is None or self.TypeCount is None:
            return 0.0
            
        for iType in range(self.nAtomTypes):
            E_Type_Total = 0.0
            
            for jType in range(self.nAtomTypes):
                ep = self.epsTable[jType, iType] * 0.25  # Factor of 4 already included in epsTable
                sig = np.sqrt(self.sigTable[jType, iType])
                
                # Calculate (sigma/rcut)^3
                LJCube = (sig / self.rCut)**3
                
                # Tail correction formula:
                # u_tail = 8/3 * pi * density * eps * sig^3 * ((1/3) * (sig/rcut)^9 - (sig/rcut)^3)
                E_Type = (8.0/3.0) * pi * self.TypeCount[jType] * ep * sig**3
                E_Type = E_Type * ((1.0/3.0) * LJCube**3 - LJCube)
                
                E_Type_Total += E_Type
            
            # Multiply by number of atoms of type i
            E_Corr += self.TypeCount[iType] * E_Type_Total
        
        return E_Corr
    
    #-------------------------------------------------------------------------
    def count_type(self, NMol: List[int]):
        """
        Corresponds to CountType_EP_LJ_Cut
        Count atoms of each type based on molecule counts
        
        Args:
            NMol: Number of molecules of each type
        """
        nMolTypes = len(NMol)
        
        if self.TypeCount is None:
            self.TypeCount = np.zeros(self.nAtomTypes, dtype=int)
        
        self.TypeCount.fill(0)
        
        for iType in range(self.nAtomTypes):
            for iMolType in range(nMolTypes):
                if self.PerMolTypeCount is not None:
                    self.TypeCount[iType] += NMol[iMolType] * self.PerMolTypeCount[iType, iMolType]
    
    #-------------------------------------------------------------------------
    def compute_type_count(self, curbox):
        """
        Corresponds to ComputeTypeCount_EP_LJ_Cut
        Calculate how many atoms of each type are in each molecule type
        
        Args:
            curbox: SimBox instance with molecular data
        """
        if hasattr(curbox, 'MolData') and curbox.MolData is not None:
            nMolTypes = len(curbox.MolData)
            
            self.PerMolTypeCount = np.zeros((self.nAtomTypes, nMolTypes), dtype=int)
            self.TypeCount = np.zeros(self.nAtomTypes, dtype=int)
            
            for iMol in range(nMolTypes):
                mol_data = curbox.MolData[iMol]
                if 'atomType' in mol_data:
                    for atomType in mol_data['atomType']:
                        if atomType < self.nAtomTypes:
                            self.PerMolTypeCount[atomType, iMol] += 1
        else:
            # Fallback for simple systems
            self.PerMolTypeCount = np.ones((self.nAtomTypes, 1), dtype=int)
            self.TypeCount = np.ones(self.nAtomTypes, dtype=int)
    #-------------------------------------------------------------------------
    def process_io(self, line: str) -> int:
        """
        Corresponds to ProcessIO_EP_LJ_Cut
        Process LJ-specific input parameters
        
        Args:
            line: Input line to process
            
        Expected formats:
        - type1 epsilon sigma rmin (single type)
        - type1 type2 epsilon sigma rmin (pair interaction)
        - tailcorrection true/false
        - rcut value
        
        Returns:
            int: Status code (0 for success, negative for error)
        """
        parts = line.strip().split()
        if not parts:
            return 0
        
        command = parts[0].lower()
        
        # Handle base EasyPair commands first
        result = super().process_io(line)
        if result == 0:
            return result
        
        # Handle LJ-specific parameters
        try:
            if len(parts) == 4:
                # Single type parameters: type1 epsilon sigma rmin
                type1 = int(parts[0]) - 1  # Convert to 0-indexed
                epsilon = float(parts[1])
                sigma = float(parts[2])
                rmin = float(parts[3])
                
                if type1 >= self.nAtomTypes or type1 < 0:
                    print(f"ERROR: Invalid atom type {type1+1}")
                    return -1
                
                # Store single-type parameters
                self.epsilon[type1] = epsilon
                self.sigma[type1] = sigma
                self.rMin[type1] = rmin
                
                # Update mixing tables using Lorentz-Berthelot rules
                for jType in range(self.nAtomTypes):
                    # Geometric mean for epsilon (with factor of 4)
                    eps_mix = 4.0 * np.sqrt(epsilon * self.epsilon[jType])
                    self.epsTable[type1, jType] = eps_mix
                    self.epsTable[jType, type1] = eps_mix
                    
                    # Arithmetic mean for sigma (stored as sigma^2)
                    sig_mix = (0.5 * (sigma + self.sigma[jType]))**2
                    self.sigTable[type1, jType] = sig_mix
                    self.sigTable[jType, type1] = sig_mix
                    
                    # Maximum for rmin (stored as rmin^2)
                    rmin_mix = max(rmin, self.rMin[jType])**2
                    self.rMinTable[type1, jType] = rmin_mix
                    self.rMinTable[jType, type1] = rmin_mix
                
                return 0
                
            elif len(parts) == 5:
                # Pair interaction: type1 type2 epsilon sigma rmin
                type1 = int(parts[0]) - 1  # Convert to 0-indexed
                type2 = int(parts[1]) - 1
                epsilon = float(parts[2])
                sigma = float(parts[3])
                rmin = float(parts[4])
                
                if (type1 >= self.nAtomTypes or type1 < 0 or 
                    type2 >= self.nAtomTypes or type2 < 0):
                    print(f"ERROR: Invalid atom types {type1+1}, {type2+1}")
                    return -1
                
                # Set pair-specific parameters (already include factor of 4 for epsilon)
                self.epsTable[type1, type2] = 4.0 * epsilon
                self.epsTable[type2, type1] = 4.0 * epsilon
                
                self.sigTable[type1, type2] = sigma**2
                self.sigTable[type2, type1] = sigma**2
                
                self.rMinTable[type1, type2] = rmin**2
                self.rMinTable[type2, type1] = rmin**2
                
                return 0
                
        except (ValueError, IndexError) as e:
            print(f"ERROR: Invalid LJ parameter format: {line}")
            print(f"Exception: {e}")
            return -1
        
        print(f"Unknown LJ command: {line}")
        return -1
    
    #-------------------------------------------------------------------------
    def prologue(self):
        """Initialize LJ force field before simulation and print parameters"""
        print(f"LJ/Cut Force Field:")
        print(f"  Cutoff radius: {self.rCut}")
        print(f"  Tail corrections: {self.usetailcorrection}")
        print(f"  Number of atom types: {self.nAtomTypes}")
        
        if self.epsilon is not None and self.sigma is not None and self.rMin is not None:
            print("  LJ Parameters:")
            for i in range(self.nAtomTypes):
                print(f"    Type {i+1}: eps={self.epsilon[i]:.4f}, sig={self.sigma[i]:.4f}, rmin={self.rMin[i]:.4f}")
    
    #-------------------------------------------------------------------------
    def epilogue(self):
        """Finalize LJ force field after simulation"""
        print("LJ/Cut Force Field simulation completed")
    
    #-------------------------------------------------------------------------
    def __str__(self):
        return f"LJ_Cut(rCut={self.rCut}, nTypes={self.nAtomTypes}, tailCorr={self.usetailcorrection})" 
    #-------------------------------------------------------------------------
#================================================================================