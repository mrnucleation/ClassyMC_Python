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


# =============================================================================
class SimpleBox(SimBox):
    """
    Simple simulation box implementation that serves as the base class for most box types.
    Corresponds to the Fortran SimpleBox type.
    """
    
    def __init__(self, molData, NMolMin=None, NMolMax=None, NMol=None, nDimensions=3):
        """Initialize the SimpleBox with molecular data"""
        super().__init__()

        self.MolData = molData
        assert isinstance(self.MolData, list), "Molecular data must be a list of Molecule_Type objects"
        
        self.atomtype_dict = {}
        for mol_type in self.MolData:
            for atom_type in mol_type.atomtypes:
                if atom_type not in self.atomtype_dict:
                    self.atomtype_dict[atom_type] = len(self.atomtype_dict)
        
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
        

        if NMolMin is None:
            NMolMin = [0] * self.nMolTypes
        else:
            self.NMolMin = np.array(NMolMin, dtype=int)
            
        if NMolMax is None:
            NMolMax = [1000] * self.nMolTypes
        else:
            self.NMolMax = np.array(NMolMax, dtype=int)
        if NMol is None:
            NMol = [0] * self.nMolTypes
        
        
        # Allocate molecule arrays
        self.hasBoundary = False  # No boundary conditions by default
        
        # Build indexing arrays
        print(f"SimpleBox initialized with {self.nMaxAtoms} max atoms, {self.maxMol} max molecules")
        
        
    # -------------------------------------------------------------------------
    def collapse_mols_to_array(self):
        """Collapse molecules to array representation"""
        pass  # Implementation would go here

    
    # -------------------------------------------------------------------------
    def load_dimension(self, boxlengths: List[float]) -> bool:
        """
        Corresponds to Simplebox_LoadDimension
        Default implementation for boxes without specific dimensions
        """
        return True
    
    # -------------------------------------------------------------------------
    def boundary(self, rx):
        """
        Corresponds to SimpleBox_Boundary
        Default boundary implementation (no boundary conditions)
        """
        return rx
    
    # -------------------------------------------------------------------------
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
            print("Computing intermolecular energy...")
            print(f"Number of energy functions: {self.EFunc[0]}")
            E_Inter, accept = self.EFunc[0].detailed_calc(self)
            print(E_Inter, accept)
            if not accept:
                return False
            self.E_Inter = E_Inter
        else:
            self.E_Inter = 0.0
        
        self.ETotal = self.E_Inter + self.E_Intra
        
        if tablecheck:
            print("Energy table check would be performed here")
        
        return True
    
    # -------------------------------------------------------------------------
    def compute_intra_energy(self):
        """
        Corresponds to SimpleBox_ComputeIntraEnergy
        Compute the total intramolecular energy
        """
        E_Total = 0.0
        accept = True
        return E_Total, accept
        
        for iType in range(self.nMolTypes):
            for iMol in range(self.NMol[iType]):
                molIndx = self.MolGlobalIndx[iType, iMol]
                E_Mol, mol_accept = self.compute_mol_intra(iType, molIndx)
                if not mol_accept:
                    return 0.0, False
                E_Total += E_Mol
        
        return E_Total, accept
    
    # -------------------------------------------------------------------------
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

    
    # -------------------------------------------------------------------------
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
    
    # -------------------------------------------------------------------------
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
    
    # -------------------------------------------------------------------------
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
    
    # -------------------------------------------------------------------------
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
    # -------------------------------------------------------------------------
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
    
    # -------------------------------------------------------------------------
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
    
    # -------------------------------------------------------------------------
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
    
    # -------------------------------------------------------------------------
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
    
    # -------------------------------------------------------------------------
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
    
    # -------------------------------------------------------------------------
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
    
    
    # -------------------------------------------------------------------------
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
    
    
    
    # -------------------------------------------------------------------------
    def update_energy(self, E_Diff):
        """
        Corresponds to SimpleBox_UpdateEnergy
        Update the total energy
        """
        self.ETotal += E_Diff
    
    # -------------------------------------------------------------------------
    def update_neigh_lists(self, disp):
        """
        Corresponds to SimpleBox_UpdateNeighLists
        Update neighbor lists after a move
        """
        if self.NeighList is not None:
            for neighList in self.NeighList:
                neighList.update(disp)
    
    # -------------------------------------------------------------------------
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
    
    # -------------------------------------------------------------------------
    def epilogue(self):
        """
        Corresponds to SimpleBox_Epilogue
        Finalize the simulation box
        """
        print(f"Finalizing SimpleBox {self.boxID}")
        print(f"Final energy: {self.ETotal}")
        print(f"Neighbor list rebuilds: {self.rebuilds}")
    
    # -------------------------------------------------------------------------
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
    
    # -------------------------------------------------------------------------
    def update(self):
        """
        Corresponds to SimpleBox_Update
        General update operations
        """
        # Update displacement tracking
        self.largestdr = np.max(self.drsq) if self.drsq is not None else 0.0
    
    # -------------------------------------------------------------------------
    def safety_check(self):
        """
        Corresponds to SimpleBox_EnergySafetyCheck
        Perform energy safety checks
        """
        if self.forceERecompute:
            print("Forcing energy recomputation...")
            self.compute_energy()
            self.forceERecompute = False
    
    # -------------------------------------------------------------------------
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
    
    # -------------------------------------------------------------------------
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
                               f"{coords[0]} {coords[1]} {coords[2]}\\n")
                        
        except IOError as e:
            print(f"Error writing dump file {filename}: {e}", file=sys.stderr)
    
    # -------------------------------------------------------------------------
    def load_coordinate(self, lines, moldata):
        """
        Load coordinates from input lines.
        Format:
        - First line: NMol (mol 1) (mol 2) .....
        - Second line: NMax (mol 1) (mol 2) .....
        - Third line: NMin (mol 1) (mol 2) .....
        - Following lines: 1 1 1 2.0 1.0 2.0 ..... (mol_type mol_index atom_index coord1 coord2 coord3 ...)
        
        Args:
            lines: List of strings containing coordinate data
            
        Returns:
            int: 0 for success, -1 for failure
        """
        if len(lines) < 3:
            print("Error: Insufficient lines for coordinate loading")
            return -1
        
        # Parse NMol, NMax, NMin lines
        nmol_line = lines[0].strip().split()
        nmax_line = lines[1].strip().split()
        nmin_line = lines[2].strip().split()
        
        # Remove the first word (NMol, NMax, NMin) and convert to integers
        if len(nmol_line) > 1:
            self.NMol = np.array([int(x) for x in nmol_line[1:]], dtype=int)
        if len(nmax_line) > 1:
            self.NMolMax = np.array([int(x) for x in nmax_line[1:]], dtype=int)
        if len(nmin_line) > 1:
            self.NMolMin = np.array([int(x) for x in nmin_line[1:]], dtype=int)
        
        # Initialize counters
        current_atom = 0
        current_mol = 0
        
        self.atoms = None
        self.MolIndx = None
        self.MolType = None
        self.MolSubIndx = None
        self.AtomSubIndx = None
        
        
        newatoms = []
        newMolType = []
        newMolIndx = []
        newMolSubIndx = []
        newAtomSubIndx = []
        newAtomType = []
        
        
        seen_mol_types = {}
        
        
        # Process coordinate lines (starting from line 3)
        for line in lines[3:]:
            line = line.strip()
            print(line)
            if not line:  # Skip empty lines
                continue
                
            parts = line.split()
            if len(parts) < 3 + self.nDimensions:  # Need at least mol_type, mol_index, atom_index + coords
                continue
                
            # Parse the line: mol_type mol_index atom_index coord1 coord2 coord3 ...
            mol_type = int(parts[0]) - 1  # Convert to 0-based indexing
            submol_index = int(parts[1]) - 1  # Convert to 0-based indexing
            atom_index = int(parts[2]) - 1  # Convert to 0-based indexing
            
            # Parse coordinates
            coords = [float(x) for x in parts[3:3+self.nDimensions]]
            coords = np.array(coords, dtype=dp).reshape(1, -1) 
            
            if coords.shape[1] != self.nDimensions:
                raise IOError(f"Error: Expected {self.nDimensions} coordinates, got {coords.shape[1]}")

            # Store coordinates
            newatoms.append(coords)
            
            if (mol_type, submol_index) not in seen_mol_types:
                seen_mol_types[(mol_type, submol_index)] = current_mol
                current_mol += 1
            
            newMolType.append(mol_type)
            newMolIndx.append(seen_mol_types[(mol_type, submol_index)])
            newMolSubIndx.append(submol_index) 
            newAtomSubIndx.append(atom_index)
            
            # Get the atom type from the molecule data
            if mol_type < len(self.MolData):
                atom_type_str = self.MolData[mol_type].atomtypes[atom_index]
                if atom_type_str not in self.atomtype_dict:
                    raise IOError(f"Error: Atom type '{atom_type_str}' not found in atomtype_dict")
                atom_type = self.atomtype_dict[atom_type_str]
            newAtomType.append(atom_type)
            current_atom += 1
        
        self.atoms = np.vstack(newatoms) 
        self.MolType = np.array(newMolType, dtype=int)
        self.MolIndx = np.array(newMolIndx, dtype=int)
        self.MolSubIndx = np.array(newMolSubIndx, dtype=int)
        self.AtomSubIndx = np.array(newAtomSubIndx, dtype=int)
        self.AtomType = np.array(newAtomType, dtype=int)
        
        self.nAtoms = current_atom
        self.nMolTotal = current_mol
        
        print(f"Loaded {self.nAtoms} atoms and {self.nMolTotal} molecules")
    # -------------------------------------------------------------------------
    def get_reduced_coords(self, real_coords):
        """
        Corresponds to SimpleBox_GetReducedCoords
        Convert real coordinates to reduced coordinates (default: no transformation)
        """
        return np.array(real_coords)
    # -------------------------------------------------------------------------
    def get_real_coords(self, reduced_coords):
        """
        Corresponds to SimpleBox_GetRealCoords  
        Convert reduced coordinates to real coordinates (default: no transformation)
        """
        return np.array(reduced_coords)
    
    # -------------------------------------------------------------------------
    def pick_random_molecule(self):
        """
        Randomly pick one molecule from the simulation box.
        
        This function implements the requested functionality to randomly select
        a single molecule from the available molecules in the box. It uses
        uniform random selection and returns comprehensive molecule information.
        
        Returns:
            dict: Dictionary containing molecule information with keys:
                - 'mol_index': Global molecule index (0-based)
                - 'mol_type': Molecule type index
                - 'mol_sub_index': Molecule sub-index within its type
                - 'atom_start': Starting atom index
                - 'atom_end': Ending atom index
                - 'n_atoms': Number of atoms in molecule
                - 'coordinates': Array of atom coordinates
                - 'center_of_mass': Center of mass coordinates
        """
        import random
        
        # Check if there are any molecules in the box
        if self.nMolTotal <= 0:
            raise ValueError("Cannot pick random molecule: box is empty")
        
        # Generate random molecule index (0 to nMolTotal-1)
        mol_global_index = random.randint(0, self.nMolTotal - 1)
        
        # Get molecule data using existing get_mol_data method
        mol_data = self.get_mol_data(mol_global_index)
        
        # Get additional information
        mol_start = mol_data['molStart']
        mol_type = mol_data['molType']
        
        # Find molecule sub-index within its type
        mol_sub_index = -1
        for i in range(self.NMol[mol_type]):
            if self.MolGlobalIndx[mol_type, i] == mol_global_index:
                mol_sub_index = i
                break
        
        # Get coordinates and center of mass
        coordinates = self.get_coordinates(slice_range=(mol_start, mol_data['molEnd']))
        center_of_mass = self.centerMass[mol_global_index, :] if self.centerMass is not None else None
        
        # Return comprehensive molecule information
        return {
            'mol_index': mol_global_index,
            'mol_type': mol_type,
            'mol_sub_index': mol_sub_index,
            'atom_start': mol_start,
            'atom_end': mol_data['molEnd'],
            'n_atoms': mol_data['nAtoms'],
            'coordinates': coordinates,
            'center_of_mass': center_of_mass
        }
# =============================================================================
