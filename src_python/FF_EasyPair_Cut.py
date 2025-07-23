"""
Easy Pair Cut Force Field Base Class
Corresponds to FF_EasyPair_Cut.f90

This is the base class for pairwise force fields that use cutoff distances.
It provides the fundamental infrastructure for neighbor list management,
energy calculations, and move implementations that most pair force fields need.
"""

import numpy as np
import sys
from abc import ABC, abstractmethod
from typing import Tuple, List, Dict, Optional, Any
from .Template_Forcefield import ForceField
from .VarPrecision import dp
from .CoordinateTypes import Displacement

#================================================================================
class EasyPairCut(ForceField):
    """
    Template class for pairwise force fields with cutoff.
    Corresponds to the Fortran EasyPair_Cut type.
    """
    
    #-------------------------------------------------------------------------
    def __init__(self):
        super().__init__()
        
        # EasyPair specific attributes
        self.usetailcorrection = False
        self.rMin = None         
        self.rMinTable = None    
        
        self.rCut = 5.0
        self.rCutSq = 25.0
        
    
    #-------------------------------------------------------------------------
    def pair_function(self, rsq: float, atmtype1: int, atmtype2: int) -> float:
        """
        Corresponds to PairFunction_EasyPair_Cut
        Template pair function - override in subclasses
        
        Args:
            rsq: Distance squared between atoms
            atmtype1: Type of first atom
            atmtype2: Type of second atom
            
        Returns:
            float: Pairwise energy
        """
        return 0.0  # Default implementation - override in subclasses
    
    #-------------------------------------------------------------------------
    def tail_correction(self, curbox, disp=None) -> float:
        """
        Corresponds to TailCorrection_EasyPair_Cut
        Template tail correction - override in subclasses
        
        Args:
            curbox: SimBox instance
            disp: Optional perturbations for delta corrections
            
        Returns:
            float: Tail correction energy
        """
        return 0.0  # Default implementation - override in subclasses
    
    #-------------------------------------------------------------------------
    def detailed_calc_fortran(self, curbox) -> Tuple[float, bool]:
        """
        Corresponds to Detailed_EasyPair_Cut
        Compute total pairwise energy for the system
        
        Args:
            curbox: SimBox instance
            
        Returns:
            tuple: (total_energy, accept_flag)
        """
        atoms = curbox.get_coordinates()
        
        E_Total = 0.0
        #curbox.ETable.fill(0.0)
        accept = True
        
        # Double loop over all atoms
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
                rx = atoms[iAtom, :] - atoms[jAtom, :]
                # Apply periodic boundary conditions
                rx = curbox.boundary(rx)
                rsq = np.sum(rx*rx)
                
                # Check if within cutoff
                if rsq < self.rCutSq:
                    atmType2 = curbox.AtomType[jAtom]
                    
                    # Check minimum distance
                    if self.rMinTable is not None:
                        rmin_ij = self.rMinTable[atmType1, atmType2]
                        if rsq < rmin_ij:
                            print(f"ERROR! Overlapping atoms found!")
                            print(f"Distance: {np.sqrt(rsq)}")
                            print(f"Atoms: {iAtom}, {jAtom}")
                            print(f"Positions: {atoms[iAtom, :]}, {atoms[jAtom, :]}")
                            return 0.0, False
                    
                    # Calculate pair energy
                    E_Pair = self.pair_function(rsq, atmType1, atmType2)
                    E_Total += E_Pair
                    
                    # Update energy table
                    #curbox.ETable[iAtom] += E_Pair
                    #curbox.ETable[jAtom] += E_Pair
        
        print(f"Total Pair Energy: {E_Total}")
        
        # Add tail corrections if enabled
        if self.usetailcorrection:
            E_Corr = self.tail_correction(curbox)
            E_Total += E_Corr
            print(f"Total Tail Corrections: {E_Corr}")
        
        return E_Total, accept
    #-------------------------------------------------------------------------
    def detailed_calc(self, curbox) -> Tuple[float, bool]:
        '''
          Reimplimentation using Python style numpy arrays instead of Fortran style
        '''
        print("Running detailed_calc with numpy arrays...")
        atoms = curbox.get_coordinates()
        molIndx = curbox.MolIndx
        
        print(atoms.shape)
        E_Total = 0.0
        accept = True
        
        for iAtom, x_atom in enumerate(atoms[:-1]):
            cut_list = atoms[iAtom+1:,:] # Exclude intramolecular interactions
            jAtomTypes = curbox.AtomType[iAtom+1:]
            include_mask = np.where(molIndx[iAtom+1:] != molIndx[iAtom])
            cut_list = cut_list[include_mask]
            jAtomTypes = jAtomTypes[include_mask]
            

            rx = np.subtract(cut_list, x_atom.reshape(1, -1))  # Calculate distance vectors
            rx = curbox.boundary(rx) # Apply periodic boundary conditions if needed
            rsq = np.sum(rx**2, axis=1).reshape(-1)  # Calculate squared distances
            # Check if any pairs are within the auto-rection rMin
            if self.rMinTable is not None:
                rmin_ij = self.rMinTable[curbox.AtomType[iAtom], jAtomTypes]
                rmin_ij = rmin_ij.reshape(-1)
                accept = np.all(rsq >= rmin_ij)
                if not accept:
                    print(f"ERROR! Overlapping atoms found for atom {iAtom} with cutoff {self.rCut}")
                    print(f"Distance: {np.sqrt(rsq)}")
                    raise ValueError(f"Overlapping atoms found for atom {iAtom} with cutoff {self.rCut}")
            within_cutoff = rsq < self.rCutSq
            if not np.any(within_cutoff):
                continue
            rsq = rsq[within_cutoff]
            jAtomTypes = jAtomTypes[within_cutoff]
            E_pair = self.pair_function(rsq, curbox.AtomType[iAtom], jAtomTypes)
            E_Total += np.sum(E_pair)
        return E_Total, accept
    #-------------------------------------------------------------------------
    def diff_calc(self, curbox, disp, tempList=None, tempNNei=None) -> Tuple[float, bool]:
        """
        Corresponds to DiffECalc_EasyPair_Cut
        Compute energy difference due to perturbations
        
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
        E_Diff = 0.0
        
        # Determine perturbation type and dispatch accordingly
        if hasattr(disp[0], 'atmIndx'):
            # Displacement move
            if hasattr(disp[0], 'newlist') and disp[0].newlist:
                E_Diff, accept = self.shift_calc_single(curbox, disp, tempList, tempNNei)
            else:
                E_Diff, accept = self.shift_calc_single(curbox, disp)
                
        elif hasattr(disp[0], 'molType') and hasattr(disp[0], 'addition'):
            # Addition move
            E_Diff, accept = self.new_calc(curbox, disp, tempList, tempNNei)
            
        elif hasattr(disp[0], 'molType') and hasattr(disp[0], 'deletion'):
            # Deletion move
            E_Diff = self.old_calc(curbox, disp)
            
        elif hasattr(disp[0], 'volNew'):
            # Volume change move
            E_Diff, accept = self.ortho_vol_calc(curbox, disp)
            
        elif hasattr(disp[0], 'newType') and hasattr(disp[0], 'oldType'):
            # Atom exchange move
            E_Diff, accept = self.atom_exchange(curbox, disp)
            
        else:
            print("Unknown Perturbation Type Encountered by EasyPair_Cut", file=sys.stderr)
            return 0.0, False
        
        # Add tail correction differences if enabled
        if self.usetailcorrection:
            E_Corr = self.tail_correction(curbox, disp)
            E_Diff += E_Corr
        
        return E_Diff, accept
     #-------------------------------------------------------------------------
    def shift_calc_single_numpy(self, curbox, disp, tempList=None, tempNNei=None) -> Tuple[float, bool]:
        """
        Corresponds to Shift_EasyPair_Cut_Single
        Compute energy change for displacement moves
        
        Args:
            curbox: SimBox instance
            disp: Displacement perturbations
            tempList: Optional neighbor list
            tempNNei: Optional neighbor counts
            
        Returns:
            tuple: (energy_change, accept_flag)
        """
        E_Diff = 0.0
        accept = True
        
        # Get atom positions
        atoms = curbox.get_coordinates()
        molIndx = curbox.MolIndx
        
        
        # Process each displaced atom
        for displacement in disp:
            if not hasattr(displacement, 'atmIndx'):
                continue
                
            atmIndx = displacement.atmIndx
            atmType1 = curbox.AtomType[atmIndx]
            
            # Calculate old energy
            E_Old = 0.0
            for jAtom in range(curbox.nMaxAtoms):
                if not curbox.is_active(jAtom) or jAtom == atmIndx:
                    continue
                if curbox.MolIndx[jAtom] == curbox.MolIndx[atmIndx]:
                    continue
                
                rx = atoms[atmIndx, :] - atoms[jAtom, :]
                rx = curbox.boundary(rx)
                rsq = np.sum(rx**2)
                
                if rsq < self.rCutSq:
                    atmType2 = curbox.AtomType[jAtom]
                    E_Old += self.pair_function(rsq, atmType1, atmType2)
            
            # Calculate new energy
            E_New = 0.0
            for jAtom in range(curbox.nMaxAtoms):
                if not curbox.is_active(jAtom) or jAtom == atmIndx:
                    continue
                if curbox.MolIndx[jAtom] == curbox.MolIndx[atmIndx]:
                    continue
                
                rx = displacement.x_new[:] - atoms[jAtom, :]
                rx = curbox.boundary(rx)
                rsq = np.sum(rx**2)
                
                if rsq < self.rCutSq:
                    atmType2 = curbox.AtomType[jAtom]
                    
                    # Check minimum distance
                    if self.rMinTable is not None:
                        rmin_ij = self.rMinTable[atmType1, atmType2]
                        if rsq < rmin_ij:
                            return 0.0, False
                    
                    E_New += self.pair_function(rsq, atmType1, atmType2)
            
            E_Diff += E_New - E_Old
        
        return E_Diff, accept
    
   
    #-------------------------------------------------------------------------
    def shift_calc_single(self, curbox, disp, tempList=None, tempNNei=None) -> Tuple[float, bool]:
        """
        Corresponds to Shift_EasyPair_Cut_Single
        Compute energy change for displacement moves
        
        Args:
            curbox: SimBox instance
            disp: Displacement perturbations
            tempList: Optional neighbor list
            tempNNei: Optional neighbor counts
            
        Returns:
            tuple: (energy_change, accept_flag)
        """
        E_Diff = 0.0
        accept = True
        
        # Get atom positions
        atoms = curbox.get_coordinates()
        
        # Process each displaced atom
        for displacement in disp:
            if not hasattr(displacement, 'atmIndx'):
                continue
                
            atmIndx = displacement.atmIndx
            atmType1 = curbox.AtomType[atmIndx]
            
            # Calculate old energy
            E_Old = 0.0
            for jAtom in range(curbox.nMaxAtoms):
                if not curbox.is_active(jAtom) or jAtom == atmIndx:
                    continue
                if curbox.MolIndx[jAtom] == curbox.MolIndx[atmIndx]:
                    continue
                
                rx = atoms[atmIndx, :] - atoms[jAtom, :]
                rx = curbox.boundary(rx)
                rsq = np.sum(rx**2)
                
                if rsq < self.rCutSq:
                    atmType2 = curbox.AtomType[jAtom]
                    E_Old += self.pair_function(rsq, atmType1, atmType2)
            
            # Calculate new energy
            E_New = 0.0
            for jAtom in range(curbox.nMaxAtoms):
                if not curbox.is_active(jAtom) or jAtom == atmIndx:
                    continue
                if curbox.MolIndx[jAtom] == curbox.MolIndx[atmIndx]:
                    continue
                
                rx = displacement.x_new[:] - atoms[jAtom, :]
                rx = curbox.boundary(rx)
                rsq = np.sum(rx**2)
                
                if rsq < self.rCutSq:
                    atmType2 = curbox.AtomType[jAtom]
                    
                    # Check minimum distance
                    if self.rMinTable is not None:
                        rmin_ij = self.rMinTable[atmType1, atmType2]
                        if rsq < rmin_ij:
                            return 0.0, False
                    
                    E_New += self.pair_function(rsq, atmType1, atmType2)
            
            E_Diff += E_New - E_Old
        
        return E_Diff, accept
    
    #-------------------------------------------------------------------------
    def new_calc(self, curbox, disp, tempList=None, tempNNei=None) -> Tuple[float, bool]:
        """
        Compute energy for addition moves
        
        Args:
            curbox: SimBox instance
            disp: Addition perturbations
            tempList: Temporary neighbor list
            tempNNei: Temporary neighbor counts
            
        Returns:
            tuple: (energy_change, accept_flag)
        """

        atoms = curbox.get_coordinates()
        
        E_Diff = 0.0
        accept = True
        
        # Loop through all displacements
        for iDisp in range(len(disp)):
            iAtom = disp[iDisp].atmIndx
            atmType1 = curbox.AtomType[iAtom]
            
            # If neighbor lists are provided, use them; otherwise loop through all atoms
            if tempList is not None and tempNNei is not None:
                listIndx = disp[iDisp].listIndex
                maxNei = tempNNei[listIndx]
                neighbor_atoms = [tempList[jNei, listIndx] for jNei in range(maxNei)]
            else:
                # Fallback: check all atoms
                neighbor_atoms = [j for j in range(curbox.nMaxAtoms) 
                                if curbox.is_active(j) and j != iAtom 
                                and curbox.MolIndx[j] != curbox.MolIndx[iAtom]]
            
            # Loop through neighbors
            for jAtom in neighbor_atoms:
                # Calculate distance components
                rx = disp[iDisp].x_new - atoms[jAtom, 0]
                
                # Apply boundary conditions
                rx = curbox.boundary(rx)
                rsq = np.sum(rx**2) 
                
                if rsq < self.rCutSq:
                    atmType2 = curbox.AtomType[jAtom]
                    
                    # Check minimum distance if rMinTable is defined
                    if self.rMinTable is not None:
                        rmin_ij = self.rMinTable[atmType2, atmType1]
                        if rsq < rmin_ij:
                            accept = False
                            return 0.0, False
                    
                    E_Pair = self.pair_function(rsq, atmType1, atmType2)
                    E_Diff += E_Pair
                    curbox.dETable[iAtom] += E_Pair
                    curbox.dETable[jAtom] += E_Pair
        
        return E_Diff, accept
    #-------------------------------------------------------------------------
    def old_calc(self, curbox, disp) -> float:
        """
        Compute energy for deletion moves
        
        Args:
            curbox: SimBox instance  
            disp: Deletion perturbations
            
        Returns:
            float: Energy change
        """
        # Implementation would depend on the specific force field
        # This is a placeholder
        return 0.0
    #-------------------------------------------------------------------------
    def ortho_vol_calc(self, curbox, disp) -> Tuple[float, bool]:
        """
        Compute energy change for volume moves
        
        Args:
            curbox: SimBox instance
            disp: Volume change perturbations
            
        Returns:
            tuple: (energy_change, accept_flag)
        """
        # Implementation would depend on the specific force field
        # This is a placeholder
        return 0.0, True
    #-------------------------------------------------------------------------
    def atom_exchange(self, curbox, disp) -> Tuple[float, bool]:
        """
        Compute energy change for atom exchange moves
        
        Args:
            curbox: SimBox instance
            disp: Exchange perturbations
            
        Returns:
            tuple: (energy_change, accept_flag)
        """
        # Implementation would depend on the specific force field
        # This is a placeholder
        return 0.0, True
    
    #-------------------------------------------------------------------------
    def single_pair(self, rsq: float, atmtype1: int, atmtype2: int) -> float:
        """
        Corresponds to SinglePair_EasyPair_Cut
        Wrapper for pair function
        """
        return self.pair_function(rsq, atmtype1, atmtype2)
    
    #-------------------------------------------------------------------------
    def many_body(self, curbox, atmtype1: int, pos1: np.ndarray, 
                  atmtypes: np.ndarray, posN: np.ndarray) -> Tuple[float, bool]:
        """
        Corresponds to ManyBody_EasyPair_Cut
        Default implementation for pairwise potentials
        """
        return 0.0, True
    
    #-------------------------------------------------------------------------
    def process_io(self, line: str) -> int:
        """
        Corresponds to ProcessIO_EasyPair_Cut
        Base implementation - override in subclasses
        
        Args:
            line: Input line to process
            
        Returns:
            int: Status code
        """
        # Basic commands that all EasyPair force fields should support
        parts = line.strip().split()
        if not parts:
            return 0
        
        command = parts[0].lower()
        
        if command == "tailcorrection" and len(parts) > 1:
            try:
                self.usetailcorrection = parts[1].lower() in ['true', '.true.', 't', '1']
                return 0
            except (ValueError, IndexError):
                return -1
        elif command == "rcut" and len(parts) > 1:
            try:
                self.rCut = float(parts[1])
                self.rCutSq = self.rCut * self.rCut
                return 0
            except (ValueError, IndexError):
                return -1
        
        return super().process_io(line)
    
    #-------------------------------------------------------------------------
    def get_cutoff(self) -> float:
        """
        Corresponds to GetCutOff_EasyPair_Cut
        Return the cutoff radius
        """
        return self.rCut 
    #-------------------------------------------------------------------------
#================================================================================