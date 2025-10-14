"""
Cubic Simulation Box
Corresponds to Box_CubicBox.f90

Implements a cubic periodic boundary condition simulation box.
"""

import numpy as np
import sys
from src_python.Box_SimpleBox import SimpleBox
from src_python.VarPrecision import dp
from src_python.CoordinateTypes import OrthoVolChange

def GetXCommand(line, position):
    """Simple implementation of GetXCommand for parsing input"""
    parts = line.split()
    if position <= len(parts):
        return parts[position-1]
    return ""

class CubeBox(SimpleBox):
    """
    Cubic simulation box with periodic boundary conditions.
    Corresponds to the Fortran CubeBox type.
    """
    #--------------------------------------------------------------------------------------------------
    def __init__(self, MolData, NMolMin=None, NMolMax=None, NMol=None, nDimensions=3):
        """Initialize cubic box, calling parent constructor"""
        super().__init__(MolData, NMolMin, NMolMax, NMol, nDimensions)
        
        # Cubic box specific attributes
        self.boxL = 0.0      # Box length
        self.boxL2 = 0.0     # Half box length (for efficiency)
        self.boxStr = "Cube"
        
    #--------------------------------------------------------------------------------------------------
    @property
    def set_boxL(self, boxL):
        """
        Corresponds to Cube_SetBoxL
        Set the box length for the cubic box
        """
        self.boxL = boxL
        self.boxL2 = boxL / 2.0
        self.volume = np.prod(self.boxL)
    #--------------------------------------------------------------------------------------------------
    @property
    def get_boxL(self):
        """
        Corresponds to Cube_GetBoxL
        Get the box length for the cubic box
        """
        return self.boxL
    #--------------------------------------------------------------------------------------------------

        
    #--------------------------------------------------------------------------------------------------   
    def load_dimension(self, boxlengths: list[float]) -> bool:
        """
        Corresponds to Cube_LoadDimension
        Parse box length from input line
        """
        try:
            self.boxL = np.full(self.nDimensions, boxlengths[0], dtype=dp)
            self.boxL2 = self.boxL / 2.0
            self.volume = np.prod(self.boxL)
            return True
        except (ValueError, IndexError):
            return False
    #--------------------------------------------------------------------------------------------------
    def get_dimensions(self):
        """
        Corresponds to Cube_GetDimensions
        Returns the box boundaries as a list of [min, max] pairs for each dimension
        """
        low = np.full(self.nDimension, -self.boxL2)
        high = np.full(self.nDimension, self.boxL2)
        return np.stack([low, high], axis=1)
    #--------------------------------------------------------------------------------------------------
    def boundary(self, rx: np.ndarray):
        """
        Corresponds to Cube_Boundary
        Apply periodic boundary conditions for cubic box

        Args:
            rx: Displacement vectors to apply PBC to

        Returns:
            Modified displacement vectors with minimum image convention
        """
        # Apply minimum image convention: map to [-L/2, L/2]
        rx = rx - self.boxL * np.round(rx / self.boxL)
        return rx
        
    #--------------------------------------------------------------------------------------------------
    
    def boundary_new(self, rx, disp):
        """
        Corresponds to Cube_BoundaryNew
        Apply boundary conditions with scaling from volume changes
        """

        #Check if all scales are the same
        if not np.all(disp.scales == disp.scales[0]):
            print("Warning! CubeBox boundary_new: scales are not all the same", file=sys.stderr)
            print("Cubic box only supports isotropic scaling", file=sys.stderr)
            raise RuntimeError("scales are not all the same")

        scale_factor = disp.scales[0]

        # Apply minimum image convention with scaled box
        scaled_boxL = self.boxL * scale_factor
        rx = rx - scaled_boxL * np.round(rx / scaled_boxL)
        return rx
    #--------------------------------------------------------------------------------------------------
    def process_io(self, line):
        """
        Corresponds to Cube_ProcessIO
        Process input commands specific to cubic box
        """
        try:
            command = GetXCommand(line, 4)  # Get 4th command
            
            if command.lower() == "boxlength":
                value_str = GetXCommand(line, 5)  # Get 5th command (value)
                self.boxL = float(value_str)
                self.boxL2 = self.boxL / 2.0
                return 0
            else:
                # Fall back to parent class processing
                return super().process_io(line)
                
        except (ValueError, IndexError):
            return super().process_io(line)
    #--------------------------------------------------------------------------------------------------
    def dump_data(self, filename):
        """
        Corresponds to Cube_DumpData
        Write box configuration to file
        """
        try:
            with open(filename, 'w') as f:
                f.write("boxtype cube\n")
                f.write(f"dimension {self.boxL}\n")
                f.write(f"molmin {' '.join(map(str, self.NMolMin))}\n")
                f.write(f"molmax {' '.join(map(str, self.NMolMax))}\n")
                f.write(f"mol {' '.join(map(str, self.NMol))}\n")
                
                # Write atom coordinates
                # Note: MolData and nMolTypes would need to be passed or accessed differently
                nMolTypes = len(self.NMol) if self.NMol is not None else 0
                
                # TODO: Need to implement proper MolData access
                # for iType in range(nMolTypes):
                #     for iMol in range(self.NMol[iType]):
                #         for iAtom in range(MolData[iType].nAtoms):
                #             # Calculate array index
                #             sub_indx = sum(self.NMolMax[:iType]) + iMol
                #             array_indx = self.MolStartIndx[sub_indx] + iAtom
                #             
                #             f.write(f"{iType+1} {iMol+1} {iAtom+1} "
                #                    f"{self.atoms[array_indx, 0]} "
                #                    f"{self.atoms[array_indx, 1]} "
                #                    f"{self.atoms[array_indx, 2]}\n")
                pass  # Placeholder until MolData is available
        except IOError as e:
            print(f"Error writing dump file {filename}: {e}", file=sys.stderr)
    #--------------------------------------------------------------------------------------------------
    def prologue(self):
        """
        Corresponds to Cube_Prologue
        Initialize and validate the cubic box setup
        """
       
        # Check initial constraints
        if self.Constrain is not None:
            for constraint in self.Constrain:
                constraint.method.prologue()
                accept = constraint.method.check_initial_constraint(self)
                if not accept:
                    print("Initial Constraints are not satisfied!", file=sys.stderr)
                    raise RuntimeError("Initial constraints failed")
        
        # Calculate total molecules
        if self.NMol is not None:
            self.nMolTotal = sum(self.NMol)
        else:
            self.nMolTotal = 0
        
        
        # Build neighbor lists
        if self.NeighList is not None:
            for i, neigh_list in enumerate(self.NeighList):
                neigh_list.build_list(i)
        
        # Compute initial energy
        self.compute_energy()
        
        # Print initialization info
        print(f"Box {self.boxID} Total Molecule Count: {self.nMolTotal}")
        print(f"Box {self.boxID} Volume: {self.volume}")
        print(f"Box {self.boxID} Number Density: {self.nMolTotal/self.volume}")
        
        # Compute center of mass for molecules
        #for iMol in range(self.maxMol):
        #    mol_data = self.get_mol_data(iMol)
        #    mol_start = mol_data['molStart']
        #    if self.MolSubIndx[mol_start] <= self.NMol[self.MolType[mol_start]]:
        #        self.compute_cm(iMol)
    
    def get_reduced_coords(self, real_coords):
        """
        Corresponds to Cube_GetReducedCoords
        Convert real coordinates to reduced (fractional) coordinates
        """
        reduced_coords = np.zeros(3)
        for i in range(self.nDimension):
            reduced_coords[i] = (real_coords[i] + self.boxL2) / self.boxL
        return reduced_coords

    def get_real_coords(self, reduced_coords):
        """
        Corresponds to Cube_GetRealCoords
        Convert reduced (fractional) coordinates to real coordinates
        """
        real_coords = np.zeros(3)
        for i in range(self.nDimension):
            real_coords[i] = self.boxL * reduced_coords[i] - self.boxL2
        return real_coords
    #--------------------------------------------------------------------------------------------------
    def update_volume(self, vol_disp):
        """
        Corresponds to Cube_UpdateVolume
        Update box dimensions after volume change
        """

        assert isinstance(vol_disp, OrthoVolChange), "CubeBox update_volume: disp must be an OrthoVolChange"

        if isinstance(vol_disp, OrthoVolChange):
            vol_ratio = vol_disp.volNew / vol_disp.volOld
            self.volume = vol_disp.volNew
            self.boxL = self.boxL * (vol_ratio ** (1.0/3.0))
            self.boxL2 = self.boxL / 2.0
    #--------------------------------------------------------------------------------------------------

    @property
    def cell(self):
        """
        Return the ASE cell matrix for the cubic box.
        For ASE, the cell matrix should be a 3x3 matrix where the diagonal elements
        are the box lengths and off-diagonal elements are the tilt factors.
        For a cubic box, we use an orthogonal cell matrix.
        """
        if self.boxL is not None:
            # Handle both scalar and array cases for boxL
            if np.isscalar(self.boxL):
                boxL_val = self.boxL
            else:
                boxL_val = self.boxL[0] if len(self.boxL) > 0 else 0.0

            if boxL_val > 0:
                if self.nDimensions == 3:
                    return np.diag([boxL_val, boxL_val, boxL_val])
                else:
                    # For other dimensions, create appropriate cell matrix
                    cell_matrix = np.zeros((3, 3))
                    for i in range(min(self.nDimensions, 3)):
                        if np.isscalar(self.boxL):
                            cell_matrix[i, i] = self.boxL
                        else:
                            cell_matrix[i, i] = self.boxL[i]
                    return cell_matrix
            else:
                return None
        else:
            return None
    #--------------------------------------------------------------------------------------------------
#==================================================================================================