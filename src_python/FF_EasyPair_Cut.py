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
from src_python.Template_Forcefield import ForceField
from src_python.VarPrecision import dp
from src_python.CoordinateTypes import Displacement, OrthoVolChange, Addition, Deletion

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
        
        self.checkOldrMin = True
    
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
        E_Diff = 0.0
        
        # Determine perturbation type and dispatch accordingly
        if isinstance(disp, Displacement):
            # Check if it's an addition or deletion move
            E_Diff, accept = self.shift_calc_single_numpy(curbox, disp, tempList, tempNNei)
        elif isinstance(disp, Addition):
                # Addition move
                E_Diff, accept = self.new_calc(curbox, disp, tempList, tempNNei)
        elif isinstance(disp, Deletion):
                # Deletion move
                E_Diff, accept = self.old_calc(curbox, disp)
        elif isinstance(disp, OrthoVolChange):
            # Orthorhombic volume move
            E_Diff, accept = self.ortho_vol_calc(curbox, disp)
        else:
            print("Unknown Perturbation Type Encountered by EasyPair_Cut", file=sys.stderr)
            return 0.0, False
        
        # Add tail correction differences if enabled
        if self.usetailcorrection:
            E_Corr = self.tail_correction(curbox, disp)
            E_Diff += E_Corr
        
        return E_Diff, accept
    #-------------------------------------------------------------------------
    def shift_calc_single_numpy(self, curbox, disp: Displacement, tempList=None, tempNNei=None) -> Tuple[float, bool]:
        """
        Corresponds to Shift_EasyPair_Cut_Single
        Compute energy change for displacement moves
        
        Args:
            curbox: SimBox instance
            disp: Displacement perturbations (single or list)
            tempList: Optional neighbor list
            tempNNei: Optional neighbor counts
            
        Returns:
            tuple: (energy_change, accept_flag)
        """
        assert isinstance(disp, Displacement), "shift_calc_single_numpy expects a Displacement displacement"
        
        E_Diff = 0.0
        accept = True
        
        # Get atom positions
        atoms = curbox.get_coordinates()
        molIndx = curbox.MolIndx
        #I will update this later when the neighbor list is implemented
        mask = np.where(molIndx != disp.molIndx)
        cut_list = atoms[mask]
        jAtomTypes = curbox.AtomType[mask]       
        atoms_new = disp.X
        
        #Compute the squared distances of the new positions
        rx = cut_list[None, :, :] - atoms_new[:, None, :]
        rx = rx.reshape(-1, 3)
        rx = curbox.boundary(rx)
        rsq = np.sum(rx**2, axis=1).reshape(-1)
        
        # Check if any pairs are within the cutoff
        if self.rMinTable is not None:
            rmin_ij = self.rMinTable[curbox.AtomType[disp.atmIndicies], jAtomTypes]
            rmin_ij = rmin_ij.reshape(-1)
            accept = np.all(rsq >= rmin_ij)
            if not accept:
                return 0.0, False

        
        within_cutoff = rsq < self.rCutSq
        rsq = rsq[within_cutoff]
        jAtomTypes_new = jAtomTypes[within_cutoff]

        E_pair = self.pair_function(rsq, curbox.AtomType[disp.atmIndicies], jAtomTypes_new)
        E_Diff += np.sum(E_pair)
        
        #Compute the squared distances of the old positions
        rx_old = atoms[disp.atmIndicies, None, :] - cut_list[None, :, :]
        rx_old = rx_old.reshape(-1, 3)
        rx_old = curbox.boundary(rx_old)
        rsq_old = np.sum(rx_old**2, axis=1).reshape(-1)
        
        within_cutoff_old = rsq_old < self.rCutSq
        rsq_old = rsq_old[within_cutoff_old]
        jAtomTypes_old = jAtomTypes[within_cutoff_old]
        E_pair_old = self.pair_function(rsq_old, curbox.AtomType[disp.atmIndicies], jAtomTypes_old)
        E_Diff -= np.sum(E_pair_old)
        return E_Diff, accept
   
    #-------------------------------------------------------------------------
    def new_calc(self, curbox, disp, tempList=None, tempNNei=None) -> Tuple[float, bool]:
        """
        Compute energy for addition moves using numpy vectorized operations

        Args:
            curbox: SimBox instance
            disp: Addition perturbations (list of Displacement objects)
            tempList: Temporary neighbor list
            tempNNei: Temporary neighbor counts

        Returns:
            tuple: (energy_change, accept_flag)
        """
        
        assert isinstance(disp, Addition), "new_calc expects an Addition displacement"
        
        E_Diff = 0.0
        accept = True
        
        # Get atom positions
        atoms = curbox.get_coordinates()
        molIndx = curbox.MolIndx
        atoms_new = disp.X

        iAtomTypes = disp.atomTypes
        jAtomTypes = curbox.AtomType
        
        #Compute the squared distances of the new positions
        rx = atoms[None, :, :] - atoms_new[:, None, :]
        rx = rx.reshape(-1, 3)
        rx = curbox.boundary(rx)
        rsq = np.sum(rx**2, axis=1).reshape(-1)
        print(f"rsq: {rsq}")
        
        
        within_cutoff = rsq < self.rCutSq
        rsq = rsq[within_cutoff]
        jAtomTypes_new = jAtomTypes[within_cutoff]


        # Check if any pairs are within the cutoff
        if self.rMinTable is not None:
            rmin_ij = self.rMinTable[iAtomTypes, jAtomTypes]
            rmin_ij = rmin_ij.reshape(-1)
            accept = np.all(rsq >= rmin_ij)
            if not accept:
                return 0.0, False



        E_pair = self.pair_function(rsq, iAtomTypes, jAtomTypes_new)
        E_Diff += np.sum(E_pair)
        

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
        
        assert isinstance(disp, Deletion), "old_calc expects a Deletion displacement"
        
        E_Diff = 0.0
        accept = True
        
        # Get atom positions
        atoms = curbox.get_coordinates()
        molIndx = curbox.MolIndx
        atom_indicies = disp.atomIndicies
        #I will update this later when the neighbor list is implemented
        mask = np.where(molIndx != disp.molIndx)
        cut_list = atoms[mask]
        old_atoms = atoms[atom_indicies]
        jAtomTypes = curbox.AtomType[mask]
        jAtomTypes_oldpos = curbox.AtomType[atom_indicies]
        
        #Compute the squared distances of the new positions
        rx = cut_list[None, :, :] - old_atoms[:, None, :]
        rx = rx.reshape(-1, 3)
        rx = curbox.boundary(rx)
        rsq = np.sum(rx**2, axis=1).reshape(-1)
        
        print(f"rsq: {rsq}")
        
        # This check is not needed for deletion moves since it means a failure of another move to check the rMin
        # Thus this is an optional check
        
        if self.checkOldrMin:
            if self.rMinTable is not None:
                rmin_ij = self.rMinTable[jAtomTypes_oldpos, jAtomTypes]
                rmin_ij = rmin_ij.reshape(-1)
                accept = np.all(rsq >= rmin_ij)
                if not accept:
                    return 0.0, False

        
        within_cutoff = rsq < self.rCutSq
        rsq = rsq[within_cutoff]
        jAtomTypes_new = jAtomTypes[within_cutoff]

        E_pair = self.pair_function(rsq, jAtomTypes_oldpos, jAtomTypes_new)
        E_Diff -= np.sum(E_pair)
        

        return E_Diff, accept
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
        # Accept either a single OrthoVolChange or a list with one element
        vol_disp = disp[0] if isinstance(disp, (list, tuple, np.ndarray)) else disp
        assert isinstance(vol_disp, OrthoVolChange), "ortho_vol_calc expects an OrthoVolChange displacement"

        atoms = curbox.get_coordinates()
        centermass = curbox.centermass
        mol_indx = curbox.MolIndx
        atom_types = curbox.AtomType
        
        E_total_new = 0.0
        accept = True

        # Loop over i atoms; vectorize over j>i
        num_atoms = atoms.shape[0]
        for i_atom in range(num_atoms - 1):
            # Exclude intramolecular interactions (j atoms in same molecule)
            mask_valid = (mol_indx[i_atom + 1:] != mol_indx[i_atom])
            if not np.any(mask_valid):
                continue
            molIndx1 = mol_indx[iAtom]
            neighbors = atoms[i_atom + 1:][mask_valid]
            j_types = atom_types[i_atom + 1:][mask_valid]
            jMolIndx = mol_indx[i_atom + 1:][mask_valid]
            
            # Move all the atoms by their respective molecular center mass
            dX_i = centermass[molIndx1]* (vol_disp.scales-1.0)
            dX_j = centermass[jMolIndx]* (vol_disp.scales-1.0)
            
            


            E_new = 0.0

            # Pair energy for selected pairs
            # Ensure array-like then sum
            E_total_new += float(np.sum(E_pairs))

        # Energy change is new total minus current intermolecular energy
        E_diff = E_total_new - getattr(curbox, 'E_Inter', 0.0)
        return E_diff, accept
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
