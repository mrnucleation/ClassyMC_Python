"""
Simple Simulation Box
Corresponds to Box_SimpleBox.f90

Wall-less simulation box which acts as the parent class to other simulation box types.
Unlike the base template this box can be used for either a boundless problem
or in conjunction with a hard wall constraint to make a non-periodic condition.
"""

import numpy as np
import sys
from abc import ABC, abstractmethod
from .Template_SimBox import SimBox
from .VarPrecision import dp
from typing import Optional, List, Dict, Tuple, Any, Union


class SimpleBox(SimBox):
    """
    Simple simulation box implementation that serves as the base class for most box types.
    Corresponds to the Fortran SimpleBox type.
    """
    
    def __init__(self, molData, NMolMin=None, NMolMax=None, NMol=None, nDimensions=3):
        """Initialize the SimpleBox with molecular data"""
        super().__init__()

        self.MolData = molData
        self.nMolTypes = len(molData)
        self.nDimensions = nDimensions
        assert self.nDimensions > 0, "Number of dimensions must be positive"

        # Core box properties
        self.boxStr = "NoBox"
        self.forceERecompute = False
        self.volume = 0.0
        self.pressure = 0.0
        self.HTotal = 0.0
        
        # Neighbor list tracking
        self.rebuilds = 0
        self.dangerbuilds = 0
        self.maxdr = 0.0
        self.maxdr2 = 0.0
        self.largestdr = 0.0
        self.rebuildsensitivity = 0.5
        
        # Arrays for tracking molecular structure (to be allocated)
        self.atoms = None
        self.dr = None
        self.drsq = None
        self.ETable = None
        self.dETable = None
        self.forces = None
        
        #Molecule List
        self.Molecules = []
        
        # Atom/Molecule indexing arrays
        self.AtomType = None
        self.MolType = None
        self.MolIndx = None
        self.MolSubIndx = None
        self.AtomSubIndx = None
        
        # Molecule boundaries
        self.centerMass = None
        

        
        # Chemical potentials
        self.chempot = np.zeros(self.nMolTypes, dtype=dp)
        
        # Temporary arrays
        self.newpos = None
        self.temppos = None
        
        # Constraint arrays
        self.Constrain = []
        
        # Energy calculation function
        self.EFunc = []
        
        # Neighbor lists
        self.NeighList = []
        

        if NMolMin is None or NMolMax is None or NMol is None:
            raise ValueError("Molecule bounds must be defined prior to box initialization!")
        
        self.NMolMin = np.array(NMolMin, dtype=int)
        self.NMolMax = np.array(NMolMax, dtype=int)
        self.NMol = np.array(NMol, dtype=int)
        
        
        # Allocate position and energy arrays
        self.atoms = np.zeros((self.nMaxAtoms, nDimensions), dtype=dp)
        self.dr = np.zeros((self.nMaxAtoms, nDimensions), dtype=dp)
        self.drsq = np.zeros(self.nMaxAtoms, dtype=dp)
        
        # Allocate indexing arrays
        self.AtomType = np.zeros(self.nMaxAtoms, dtype=int)
        self.MolType = np.zeros(self.nMaxAtoms, dtype=int)
        self.MolIndx = np.zeros(self.nMaxAtoms, dtype=int)
        self.MolSubIndx = np.zeros(self.nMaxAtoms, dtype=int)
        self.AtomSubIndx = np.zeros(self.nMaxAtoms, dtype=int)
        
        # Allocate molecule arrays
        self.MolStartIndx = np.zeros(self.maxMol, dtype=int)
        self.MolEndIndx = np.zeros(self.maxMol, dtype=int)
        self.centerMass = np.zeros((self.maxMol, nDimensions), dtype=dp)
        
        
        self.hasBoundary = False  # No boundary conditions by default
        
        # Build indexing arrays
        print(f"SimpleBox initialized with {self.nMaxAtoms} max atoms, {self.maxMol} max molecules")
        
    def _build_indexing_arrays(self):
        """Build indexing arrays for the box"""
        pass
    
    def load_dimension(self, line):
        """
        Corresponds to Simplebox_LoadDimension
        Default implementation for boxes without specific dimensions
        """
        return 0  # Success by default for simple box
    
    
    def boundary(self, rx):
        """
        Corresponds to SimpleBox_Boundary
        Default boundary implementation (no boundary conditions)
        """
        return rx
    
    def compute_energy(self, tablecheck=False):
        """
        Corresponds to SimpleBox_ComputeEnergy
        Compute the total energy of the system
        """
        # Compute intramolecular energy
        E_Intra, accept = self.compute_intra_energy()
        if not accept:
            return False
        
        self.E_Intra = E_Intra
        
        # Compute intermolecular energy using energy function
        if self.EFunc is not None:
            E_Inter, accept = self.EFunc.detailed_calc(self)
            if not accept:
                return False
            self.E_Inter = E_Inter
        else:
            self.E_Inter = 0.0
        
        self.ETotal = self.E_Inter + self.E_Intra
        
        if tablecheck:
            print("Energy table check would be performed here")
        
        return True
    
    def compute_intra_energy(self):
        """
        Corresponds to SimpleBox_ComputeIntraEnergy
        Compute the total intramolecular energy
        """
        E_Total = 0.0
        accept = True
        
        for iType in range(self.nMolTypes):
            for iMol in range(self.NMol[iType]):
                molIndx = self.MolGlobalIndx[iType, iMol]
                E_Mol, mol_accept = self.compute_mol_intra(iType, molIndx)
                if not mol_accept:
                    return 0.0, False
                E_Total += E_Mol
        
        return E_Total, accept
    
    def compute_mol_intra(self, molType, molIndx, newpos=None, fragment=None):
        """
        Corresponds to SimpleBox_ComputeMolIntra
        Compute intramolecular energy for a single molecule
        """
        if self.MolData[molType].get('ridgid', False):
            return 0.0, True
        
        E_Intra = 0.0
        molStart = self.MolStartIndx[molIndx]
        molEnd = self.MolEndIndx[molIndx]
        nAtoms = self.MolData[molType]['nAtoms']
        
        # Use provided positions or current positions
        if newpos is not None:
            molpos = newpos[:nAtoms, :]
        else:
            molpos = self.atoms[molStart:molEnd+1, :]
        
        # Determine which atoms to include
        if fragment is not None:
            included = np.zeros(nAtoms, dtype=bool)
            included[fragment] = True
        else:
            included = np.ones(nAtoms, dtype=bool)
        
        # Compute bond energies
        if 'bonds' in self.MolData[molType]:
            for bond in self.MolData[molType]['bonds']:
                mem1, mem2 = bond['mem1'], bond['mem2']
                if fragment is None or (included[mem1] and included[mem2]):
                    bond_coords = molpos[:, [mem1, mem2]]
                    E_Bond, accept = self._compute_bond_energy(bond, bond_coords)
                    if not accept:
                        return 0.0, False
                    E_Intra += E_Bond
        
        # Compute angle energies
        if 'angles' in self.MolData[molType]:
            for angle in self.MolData[molType]['angles']:
                mem1, mem2, mem3 = angle['mem1'], angle['mem2'], angle['mem3']
                if fragment is None or (included[mem1] and included[mem2] and included[mem3]):
                    angle_coords = molpos[:, [mem1, mem2, mem3]]
                    E_Angle, accept = self._compute_angle_energy(angle, angle_coords)
                    if not accept:
                        return 0.0, False
                    E_Intra += E_Angle
        
        # Compute torsion energies
        if 'torsions' in self.MolData[molType]:
            for torsion in self.MolData[molType]['torsions']:
                mem1, mem2, mem3, mem4 = torsion['mem1'], torsion['mem2'], torsion['mem3'], torsion['mem4']
                if fragment is None or (included[mem1] and included[mem2] and included[mem3] and included[mem4]):
                    torsion_coords = molpos[:, [mem1, mem2, mem3, mem4]]
                    E_Torsion, accept = self._compute_torsion_energy(torsion, torsion_coords)
                    if not accept:
                        return 0.0, False
                    E_Intra += E_Torsion
        
        return E_Intra, True
    
    def _compute_bond_energy(self, bond, coords):
        """Compute energy for a bond interaction"""
        # Placeholder - would call appropriate bond force field
        return 0.0, True
    
    def _compute_angle_energy(self, angle, coords):
        """Compute energy for an angle interaction"""
        # Placeholder - would call appropriate angle force field
        return 0.0, True
    
    def _compute_torsion_energy(self, torsion, coords):
        """Compute energy for a torsion interaction"""
        # Placeholder - would call appropriate torsion force field
        return 0.0, True
    
    def compute_energy_delta(self, disp, templist=None, tempnnei=None, computeintra=None):
        """
        Corresponds to SimpleBox_ComputeEnergyDelta
        Compute energy change due to a perturbation
        """
        E_Inter = 0.0
        E_Intra = 0.0
        accept = True
        
        # Compute intermolecular energy change
        if self.EFunc is not None:
            E_Inter, accept = self.EFunc.diff_calc(self, disp, templist, tempnnei)
            if not accept:
                return 0.0, 0.0, False
        
        # Compute intramolecular energy change if requested
        if computeintra is None or computeintra:
            E_Intra, accept = self.compute_intra_energy_delta(disp)
            if not accept:
                return 0.0, 0.0, False
        
        return E_Inter, E_Intra, accept
    
    def compute_intra_energy_delta(self, disp):
        """
        Corresponds to SimpleBox_ComputeIntraEnergyDelta
        Compute change in intramolecular energy due to displacement
        """
        E_Intra = 0.0
        
        if hasattr(disp[0], 'molType') and hasattr(disp[0], 'molIndx'):
            molType = disp[0].molType
            molIndx = disp[0].molIndx
            
            # Get molecule data
            molStart = self.MolStartIndx[molIndx]
            molEnd = self.MolEndIndx[molIndx]
            nAtoms = self.MolData[molType]['nAtoms']
            
            # Create new positions
            self.newpos[:nAtoms, :] = self.atoms[molStart:molEnd+1, :]
            
            # Apply displacements
            for iDisp, displacement in enumerate(disp):
                if hasattr(displacement, 'atmIndx'):
                    atomsubindx = self.AtomSubIndx[displacement.atmIndx]
                    self.newpos[atomsubindx, 0] = displacement.x_new
                    self.newpos[atomsubindx, 1] = displacement.y_new
                    self.newpos[atomsubindx, 2] = displacement.z_new
            
            # Compute energy with new positions
            E_new, accept = self.compute_mol_intra(molType, molIndx, newpos=self.newpos)
            if not accept:
                return 0.0, False
            
            # Compute energy with old positions
            E_old, accept = self.compute_mol_intra(molType, molIndx)
            if not accept:
                return 0.0, False
            
            E_Intra = E_new - E_old
        
        return E_Intra, True
    
    def get_coordinates(self, slice_range=None):
        """
        Corresponds to SimpleBox_GetCoordinates  
        Get atomic coordinates, optionally for a slice
        """
        if slice_range is None:
            return self.atoms
        else:
            start, end = slice_range
            return self.atoms[start:end+1, :]
    
    def get_mol_data(self, global_indx, **kwargs):
        """
        Corresponds to SimpleBox_GetMolData
        Get molecular data for a given global molecule index
        """
        if global_indx >= self.maxMol or global_indx < 0:
            raise IndexError(f"Invalid molecule index: {global_indx}")
        
        molStart = self.MolStartIndx[global_indx]
        molEnd = self.MolEndIndx[global_indx]
        molType = self.MolType[molStart]
        nAtoms = molEnd - molStart + 1
        
        result = {
            'molStart': molStart,
            'molEnd': molEnd,
            'molType': molType,
            'nAtoms': nAtoms
        }
        
        # Return specific values if requested
        for key, value in kwargs.items():
            if key in result:
                return result[key]
        
        return result
    
    def is_active(self, atom_indx):
        """
        Corresponds to SimpleBox_IsActive
        Check if an atom index is currently active (has a molecule)
        """
        if atom_indx >= self.nMaxAtoms or atom_indx < 0:
            return False
        
        molIndx = self.MolIndx[atom_indx]
        molType = self.MolType[atom_indx]
        molSubIndx = self.MolSubIndx[atom_indx]
        
        if molType >= self.nMolTypes or molType < 0:
            return False
        
        return molSubIndx < self.NMol[molType]
    
    def add_mol(self, molType, coords):
        """
        Corresponds to SimpleBox_AddMol
        Add a molecule to the simulation box
        """
        if self.NMol[molType] >= self.NMolMax[molType]:
            return False  # Box is full for this type
        
        # Find the molecule index
        molIndx = self.MolGlobalIndx[molType, self.NMol[molType]]
        molStart = self.MolStartIndx[molIndx]
        molEnd = self.MolEndIndx[molIndx]
        
        # Set the coordinates
        self.atoms[molStart:molEnd+1, :] = coords
        
        # Update counters
        self.NMol[molType] += 1
        self.nMolTotal += 1
        self.nAtoms += self.MolData[molType]['nAtoms']
        
        return True
    
    def delete_mol(self, molType, molSubIndx):
        """
        Corresponds to SimpleBox_DeleteMol
        Remove a molecule from the simulation box
        """
        if molSubIndx >= self.NMol[molType] or molSubIndx < 0:
            return False
        
        # Find molecule to delete
        molIndx = self.MolGlobalIndx[molType, molSubIndx]
        
        # If not the last molecule, swap with last molecule
        if molSubIndx < self.NMol[molType] - 1:
            lastMolIndx = self.MolGlobalIndx[molType, self.NMol[molType] - 1]
            
            # Swap molecule positions
            molStart = self.MolStartIndx[molIndx]
            molEnd = self.MolEndIndx[molIndx]
            lastMolStart = self.MolStartIndx[lastMolIndx]
            lastMolEnd = self.MolEndIndx[lastMolIndx]
            
            self.atoms[molStart:molEnd+1, :] = self.atoms[lastMolStart:lastMolEnd+1, :]
            
            # Update global index
            self.MolGlobalIndx[molType, molSubIndx] = molIndx
            self.MolGlobalIndx[molType, self.NMol[molType] - 1] = lastMolIndx
        
        # Update counters
        self.NMol[molType] -= 1
        self.nMolTotal -= 1
        self.nAtoms -= self.MolData[molType]['nAtoms']
        
        return True
    
    def update_position(self, disp):
        """
        Corresponds to SimpleBox_UpdatePosition
        Update atomic positions based on displacement
        """
        for displacement in disp:
            if hasattr(displacement, 'atmIndx'):
                atmIndx = displacement.atmIndx
                self.atoms[atmIndx, 0] = displacement.x_new
                self.atoms[atmIndx, 1] = displacement.y_new
                self.atoms[atmIndx, 2] = displacement.z_new
    
    def check_constraint(self, disp=None):
        """
        Corresponds to SimpleBox_CheckConstraint
        Check if constraints are satisfied
        """
        if self.Constrain is None:
            return True
        
        for constraint in self.Constrain:
            if disp is not None:
                if not constraint.diff_check(self, disp):
                    return False
            else:
                if not constraint.check_initial_constraint(self):
                    return False
        
        return True
    
    def check_post_energy(self, disp, E_Diff):
        """
        Corresponds to SimpleBox_CheckPostEnergy
        Check constraints that depend on energy
        """
        if self.Constrain is None:
            return True
        
        for constraint in self.Constrain:
            if hasattr(constraint, 'post_energy'):
                if not constraint.post_energy(self, disp, E_Diff):
                    return False
        
        return True
    
    def compute_cm(self, molIndx):
        """
        Corresponds to SimpleBox_ComputeCM
        Compute center of mass for a molecule
        """
        molStart = self.MolStartIndx[molIndx]
        molEnd = self.MolEndIndx[molIndx]
        molType = self.MolType[molStart]
        
        total_mass = 0.0
        cm = np.zeros(3)
        
        for iAtom in range(molStart, molEnd + 1):
            atomType = self.AtomType[iAtom]
            mass = self.MolData[molType]['masses'][iAtom - molStart]  # Assuming masses are stored
            total_mass += mass
            cm += mass * self.atoms[iAtom, :]
        
        if total_mass > 0:
            self.centerMass[molIndx, :] = cm / total_mass
    
    def get_thermo(self, thermo_id):
        """
        Get thermodynamic property by ID
        """
        thermo_map = {
            0: self.ETotal,
            1: self.volume,
            2: self.pressure,
            3: self.temperature,
            4: self.HTotal
        }
        return thermo_map.get(thermo_id, 0.0)
    
    def thermo_lookup(self, property_name):
        """
        Get thermodynamic property ID by name
        """
        name_map = {
            'energy': 0,
            'volume': 1,
            'pressure': 2,
            'temperature': 3,
            'enthalpy': 4
        }
        return name_map.get(property_name.lower(), -1)
    
    def count_atoms(self, molType=None):
        """
        Corresponds to SimpleBox_CountAtoms
        Count atoms of a specific type or all atoms
        """
        if molType is None:
            return self.nAtoms
        
        count = 0
        for iType in range(self.nMolTypes):
            if iType == molType:
                count += self.NMol[iType] * self.MolData[iType]['nAtoms']
        
        return count
    
    def get_max_atoms(self):
        """Get maximum number of atoms"""
        return self.nMaxAtoms
    
    def get_max_mol(self):
        """Get maximum number of molecules"""  
        return self.maxMol
    
    def get_min_mol(self, molType):
        """Get minimum molecules for a type"""
        return self.NMolMin[molType]
    
    def update_volume(self, new_volume):
        """
        Corresponds to SimpleBox_UpdateVolume
        Update the box volume
        """
        self.volume = new_volume
    
    def update_energy(self, E_Diff):
        """
        Corresponds to SimpleBox_UpdateEnergy
        Update the total energy
        """
        self.ETotal += E_Diff
    
    def update_neigh_lists(self, disp):
        """
        Corresponds to SimpleBox_UpdateNeighLists
        Update neighbor lists after a move
        """
        if self.NeighList is not None:
            for neighList in self.NeighList:
                neighList.update(disp)
    
    def prologue(self):
        """
        Corresponds to SimpleBox_Prologue
        Initialize the simulation box
        """
        print(f"Initializing SimpleBox {self.boxID}")
        
        # Check constraints
        if not self.check_constraint():
            raise RuntimeError("Initial constraints not satisfied!")
        
        # Compute initial energy
        if not self.compute_energy():
            raise RuntimeError("Initial energy computation failed!")
        
        # Build neighbor lists
        if self.NeighList is not None:
            for i, neighList in enumerate(self.NeighList):
                neighList.build_list(i)
        
        # Compute center of mass for all molecules
        for iMol in range(self.maxMol):
            molStart = self.MolStartIndx[iMol]
            if self.MolSubIndx[molStart] < self.NMol[self.MolType[molStart]]:
                self.compute_cm(iMol)
        
        print(f"Box {self.boxID} initialized with {self.nMolTotal} molecules")
        print(f"Box {self.boxID} Total Energy: {self.ETotal}")
    
    def epilogue(self):
        """
        Corresponds to SimpleBox_Epilogue
        Finalize the simulation box
        """
        print(f"Finalizing SimpleBox {self.boxID}")
        print(f"Final energy: {self.ETotal}")
        print(f"Neighbor list rebuilds: {self.rebuilds}")
    
    def maintenance(self):
        """
        Corresponds to SimpleBox_Maintenance
        Perform maintenance operations
        """
        # Check if neighbor lists need rebuilding
        if self.largestdr > self.rebuildsensitivity:
            if self.NeighList is not None:
                for neighList in self.NeighList:
                    neighList.rebuild()
            self.rebuilds += 1
            self.largestdr = 0.0
    
    def update(self):
        """
        Corresponds to SimpleBox_Update
        General update operations
        """
        # Update displacement tracking
        self.largestdr = np.max(self.drsq) if self.drsq is not None else 0.0
    
    def safety_check(self):
        """
        Corresponds to SimpleBox_EnergySafetyCheck
        Perform energy safety checks
        """
        if self.forceERecompute:
            print("Forcing energy recomputation...")
            self.compute_energy()
            self.forceERecompute = False
    
    def process_io(self, line):
        """
        Corresponds to SimpleBox_ProcessIO
        Process input/output commands
        """
        parts = line.strip().split()
        if not parts:
            return 0
        
        command = parts[0].lower()
        
        # Handle common box commands
        if command == "volume" and len(parts) > 1:
            try:
                self.volume = float(parts[1])
                return 0
            except ValueError:
                return -1
        elif command == "pressure" and len(parts) > 1:
            try:
                self.pressure = float(parts[1])
                return 0
            except ValueError:
                return -1
        elif command == "temperature" and len(parts) > 1:
            try:
                self.temperature = float(parts[1])
                return 0
            except ValueError:
                return -1
        
        # If not handled here, return to parent class
        return super().process_io(line)
    
    def dump_data(self, filename):
        """
        Corresponds to SimpleBox_DumpData
        Write simulation data to file
        """
        try:
            with open(filename, 'w') as f:
                f.write(f"boxtype simplebox\n")
                f.write(f"volume {self.volume}\n")
                f.write(f"pressure {self.pressure}\n")
                f.write(f"temperature {self.temperature}\n")
                f.write(f"energy {self.ETotal}\n")
                f.write(f"nmol {' '.join(map(str, self.NMol))}\n")
                
                # Write atomic coordinates
                for iAtom in range(self.nMaxAtoms):
                    if self.is_active(iAtom):
                        molType = self.MolType[iAtom]
                        molSubIndx = self.MolSubIndx[iAtom]
                        atomSubIndx = self.AtomSubIndx[iAtom]
                        coords = self.atoms[iAtom, :]
                        f.write(f"{molType+1} {molSubIndx+1} {atomSubIndx+1} "
                               f"{coords[0]} {coords[1]} {coords[2]}\n")
                        
        except IOError as e:
            print(f"Error writing dump file {filename}: {e}", file=sys.stderr)
    
    def get_reduced_coords(self, real_coords):
        """
        Corresponds to SimpleBox_GetReducedCoords
        Convert real coordinates to reduced coordinates (default: no transformation)
        """
        return np.array(real_coords)
    
    def get_real_coords(self, reduced_coords):
        """
        Corresponds to SimpleBox_GetRealCoords  
        Convert reduced coordinates to real coordinates (default: no transformation)
        """
        return np.array(reduced_coords)
    
    def __del__(self):
        """Destructor equivalent"""
        # Python handles memory management automatically
        pass