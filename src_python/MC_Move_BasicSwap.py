"""
Basic Swap Monte Carlo Move
Corresponds to Move_MC_BasicSwap.f90

A simple Monte Carlo move that adds or removes molecules from the simulation box.
This is the most basic form of grand canonical Monte Carlo.
"""

import numpy as np
import sys
import math
from random import random
from Template_MCMove import MCMove
from CoordinateTypes import Addition, Deletion
from VarPrecision import dp


class BasicSwap(MCMove):
    """
    Basic swap Monte Carlo move for grand canonical simulations.
    Corresponds to the Fortran BasicSwap type.
    """
    
    def __init__(self):
        super().__init__()
        
        # Basic swap parameters
        self.verbose = True
        self.proportional = True
        self.tuneMax = True
        self.limit = 5.0
        self.targAccpt = 50.0
        
        # Rejection counters
        self.ovlaprej = 0
        self.constrainrej = 0
        self.detailedrej = 0
        
        # Displacement arrays
        self.newPart = []  # Will be allocated in prologue
        self.oldPart = [Deletion()]
    
    def constructor(self):
        """
        Corresponds to BasicSwap_Constructor
        Initialize the move
        """
        # Initialize box probabilities if not set
        if len(self.boxProb) == 0:
            # Default to single box
            self.boxProb = np.array([1.0], dtype=dp)
    
    def full_move(self, trial_box):
        """
        Corresponds to BasicSwap_FullMove
        Perform a basic swap move (either add or remove molecule)
        
        Args:
            trial_box: SimBox instance
            
        Returns:
            bool: Whether the move was accepted
        """
        if random() > 0.5:
            return self.swap_in(trial_box)
        else:
            return self.swap_out(trial_box)
    
    def swap_in(self, trial_box):
        """
        Corresponds to BasicSwap_SwapIn
        Perform a swap-in move (add molecule)
        
        Args:
            trial_box: SimBox instance
            
        Returns:
            bool: Whether the move was accepted
        """
        self.atmps += 1.0
        accept = True
        
        # Select molecule type to add
        nMolTypes = len(trial_box.MolData) if hasattr(trial_box, 'MolData') else 1
        nType = int(nMolTypes * random())
        
        if trial_box.NMol[nType] + 1 > trial_box.NMolMax[nType]:
            return False
        
        # Get molecule data
        nMove = trial_box.NMol[nType] + 1
        nMove = trial_box.MolGlobalIndx[nType, nMove - 1]
        mol_data = trial_box.get_mol_data(nMove)
        molType = mol_data['molType']
        molStart = mol_data['molStart']
        molEnd = mol_data['molEnd']
        
        # Set up new molecule atoms
        nAtoms = molEnd - molStart + 1
        for iAtom in range(nAtoms):
            atomIndx = molStart + iAtom - 1
            self.newPart[iAtom].molType = nType
            self.newPart[iAtom].molIndx = nMove
            self.newPart[iAtom].atmIndx = atomIndx
        
        # Generate random position for new molecule
        # This is a simplified implementation
        for iAtom in range(nAtoms):
            self.newPart[iAtom].x_new = trial_box.boxL * (2.0 * random() - 1.0)
            self.newPart[iAtom].y_new = trial_box.boxL * (2.0 * random() - 1.0)
            self.newPart[iAtom].z_new = trial_box.boxL * (2.0 * random() - 1.0)
        
        # Check constraints
        if not trial_box.check_constraint(self.newPart[:nAtoms]):
            self.constrainrej += 1
            return False
        
        # Calculate energy
        e_inter, e_intra, e_diff, accept = trial_box.compute_energy_delta(
            self.newPart[:nAtoms], self.tempList, self.tempNnei, computeintra=True
        )
        if not accept:
            self.ovlaprej += 1
            return accept
        
        # Check post-energy constraints
        if not trial_box.check_post_energy(self.newPart[:nAtoms], e_diff):
            self.constrainrej += 1
            return False
        
        # Calculate acceptance probability
        # For basic swap, this is just the standard grand canonical acceptance
        Prob = 1.0  # Simplified - would need proper calculation
        
        # Get extra terms (chemical potential)
        from CommonSampling import sampling
        if sampling is not None:
            extra_terms = sampling.get_extra_terms(self.newPart[:nAtoms], trial_box)
        else:
            extra_terms = 0.0
        
        # Accept/reject
        if sampling is not None:
            accept = sampling.make_decision(trial_box, e_diff, self.newPart[:nAtoms], in_prob=Prob, extra_in=extra_terms)
        else:
            accept = True
        
        if accept:
            self.accpt += 1.0
            trial_box.update_energy(e_diff, e_inter, e_intra)
            trial_box.update_position(self.newPart[:nAtoms], self.tempList, self.tempNnei)
        else:
            self.detailedrej += 1
        
        return accept
    
    def swap_out(self, trial_box):
        """
        Corresponds to BasicSwap_SwapOut
        Perform a swap-out move (remove molecule)
        
        Args:
            trial_box: SimBox instance
            
        Returns:
            bool: Whether the move was accepted
        """
        self.atmps += 1.0
        accept = True
        
        # Select molecule to remove
        nTarget = self.uniform_molecule_select(trial_box)
        mol_data = trial_box.get_mol_data(nTarget)
        molType = mol_data['molType']
        
        if trial_box.NMol[molType] - 1 < trial_box.NMolMin[molType]:
            return False
        
        self.oldPart[0].molType = molType
        self.oldPart[0].molIndx = nTarget
        
        # Check constraints
        if not trial_box.check_constraint(self.oldPart):
            self.constrainrej += 1
            return False
        
        # Calculate energy
        e_inter, e_intra, e_diff, accept = trial_box.compute_energy_delta(
            self.oldPart, self.tempList, self.tempNnei, computeintra=True
        )
        if not accept:
            self.ovlaprej += 1
            return accept
        
        # Check post-energy constraints
        if not trial_box.check_post_energy(self.oldPart, e_diff):
            self.constrainrej += 1
            return False
        
        # Calculate acceptance probability
        Prob = 1.0  # Simplified - would need proper calculation
        
        # Get extra terms (chemical potential)
        from CommonSampling import sampling
        if sampling is not None:
            extra_terms = sampling.get_extra_terms(self.oldPart, trial_box)
        else:
            extra_terms = 0.0
        
        # Accept/reject
        if sampling is not None:
            accept = sampling.make_decision(trial_box, e_diff, self.oldPart, in_prob=Prob, extra_in=extra_terms)
        else:
            accept = True
        
        if accept:
            self.accpt += 1.0
            trial_box.update_energy(e_diff, e_inter, e_intra)
            trial_box.delete_mol(self.oldPart[0].molIndx)
        else:
            self.detailedrej += 1
        
        return accept
    
    def prologue(self):
        """
        Corresponds to BasicSwap_Prologue
        Initialize the move before simulation
        """
        # Allocate newPart array
        maxAtoms = 10  # Default value - would need proper calculation
        # Note: In a real implementation, this would be calculated from MolData
        
        self.newPart = [Addition() for _ in range(maxAtoms)]
        self.create_temp_array(maxAtoms)
        
        print(f"(Basic Swap) Move initialized")
    
    def epilogue(self):
        """
        Corresponds to BasicSwap_Epilogue
        Finalize the move after simulation
        """
        print(f"Basic Swap Moves Accepted: {int(self.accpt):15d}")
        print(f"Basic Swap Moves Attempted: {int(self.atmps):15d}")
        
        accpt_rate = self.get_accept_rate()
        print(f"Basic Swap Acceptance Rate: {accpt_rate:15.8f}")
        
        if self.verbose:
            print(f"Basic Swap, Rejections due to overlap: {self.ovlaprej:15d}")
            print(f"Basic Swap, Rejections due to constraint: {self.constrainrej:15d}")
            print(f"Basic Swap, Rejections due to detailed balance: {self.detailedrej:15d}")
    
    def process_io(self, line: str) -> int:
        """
        Corresponds to BasicSwap_ProcessIO
        Process input commands
        """
        parts = line.strip().split()
        if len(parts) < 5:
            return -1
        
        command = parts[3].lower()
        value_str = parts[4]
        
        try:
            if command == "verbose":
                self.verbose = value_str.lower() in ['true', '.true.', 't', '1']
            elif command == "proportional":
                self.proportional = value_str.lower() in ['true', '.true.', 't', '1']
            elif command == "updatefreq":
                self.maintFreq = int(value_str)
            else:
                return -1
        except (ValueError, IndexError):
            return -1
        
        return 0 