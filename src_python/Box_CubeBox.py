"""
Cubic Simulation Box
Corresponds to Box_CubicBox.f90

Implements a cubic periodic boundary condition simulation box.
"""

import numpy as np
import sys
from .Box_SimpleBox import SimpleBox
from .VarPrecision import dp

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
    
    def __init__(self, MolData, NMolMin=None, NMolMax=None, NMol=None, nDimensions=3):
        """Initialize cubic box, calling parent constructor"""
        super().__init__(MolData, NMolMin, NMolMax, NMol, nDimensions)
        
        # Cubic box specific attributes
        self.boxL = 0.0      # Box length
        self.boxL2 = 0.0     # Half box length (for efficiency)
        self.boxStr = "Cube"
        
        
        
    #--------------------------------------------------------------------------   
    def load_dimension(self, boxlengths: list[float]) -> bool:
        """
        Corresponds to Cube_LoadDimension
        Parse box length from input line
        """
        try:
            self.boxL = np.full(self.nDimensions, boxlengths[0], dtype=dp)
            self.boxL2 = self.boxL / 2.0
            self.volume = self.boxL ** (self.nDimensions)
            return True
        except (ValueError, IndexError):
            return False
    
    def get_dimensions(self):
        """
        Corresponds to Cube_GetDimensions
        Returns the box boundaries as a list of [min, max] pairs for each dimension
        """
        dimensions = []
        for i in range(self.nDimension):
            dimensions.append([-self.boxL2, self.boxL2])
        return dimensions
    
    def boundary(self, rx: np.ndarray):
        """
        Corresponds to Cube_Boundary
        Apply periodic boundary conditions for cubic box
        
        Args:
            rx, ry, rz: Coordinates to apply PBC to
            
        Returns:
            Modified coordinates (rx, ry, rz) or just rx if others are None
        """
        rx = np.where(rx > self.boxL2, rx - self.boxL, rx)
        rx = np.where(rx < -self.boxL2, rx + self.boxL, rx)
        return rx
        

    
    def boundary_new(self, rx, disp):
        """
        Corresponds to Cube_BoundaryNew
        Apply boundary conditions with scaling from volume changes
        """
        # Extract scale factor from displacement
        scale_factor = 1.0
        if hasattr(disp[0], 'xScale'):
            scale_factor = disp[0].xScale
        
        rx = np.where(rx > self.boxL2*scale_factor, rx - self.boxL*scale_factor, rx)
        rx = np.where(rx < -self.boxL2*scale_factor, rx + self.boxL*scale_factor, rx)
        return rx
    
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
    
    def prologue(self):
        """
        Corresponds to Cube_Prologue
        Initialize and validate the cubic box setup
        """
        # Check if particles are within bounds
        for iAtom in range(self.nMaxAtoms):
            if not self.is_active(iAtom):
                continue
                
            for iDim in range(self.nDimension):
                if abs(self.atoms[iAtom, iDim]) > self.boxL2:
                    print(f"Warning! Particle out of bounds!", file=sys.stderr)
                    print(f"Particle Number: {iAtom}", file=sys.stderr)
                    print(f"Box Length: {self.boxL}", file=sys.stderr)
                    print(f"Position: {self.atoms[iAtom, :]}", file=sys.stderr)
                    raise RuntimeError("Particle out of bounds")
        
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
        
        # Set volume
        self.volume = self.boxL ** 3
        
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
        for iMol in range(self.maxMol):
            mol_data = self.get_mol_data(iMol)
            mol_start = mol_data['molStart']
            if self.MolSubIndx[mol_start] <= self.NMol[self.MolType[mol_start]]:
                self.compute_cm(iMol)
    
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
    
    def update_volume(self, disp):
        """
        Corresponds to Cube_UpdateVolume
        Update box dimensions after volume change
        """
        if hasattr(disp[0], 'volNew') and hasattr(disp[0], 'volOld'):
            vol_ratio = disp[0].volNew / disp[0].volOld
            self.volume = disp[0].volNew
            self.boxL = self.boxL * (vol_ratio ** (1.0/3.0))
            self.boxL2 = self.boxL / 2.0
