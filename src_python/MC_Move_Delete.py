"""
Delete Monte Carlo Move
Corresponds to Move_MC_Delete.f90

A simple Monte Carlo move that removes molecules from the simulation box.
This is useful for grand canonical simulations where you want to control
the deletion process separately from insertion.
"""

import numpy as np
import sys
import math
from random import random
from Template_MCMove import MCMove
from CoordinateTypes import Deletion
from VarPrecision import dp


class Delete(MCMove):
    """
    Delete Monte Carlo move for removing molecules.
    Corresponds to the Fortran Delete type.
    """
    
    def __init__(self):
        super().__init__()
        
        # Delete move parameters
        self.verbose = True
        self.proportional = True
        
        # Rejection counters
        self.ovlaprej = 0
        self.constrainrej = 0
        self.detailedrej = 0
        
        # Displacement array
        self.oldPart = [Deletion()]
    
    def constructor(self):
        """
        Corresponds to Delete_Constructor
        Initialize the move
        """
        # Initialize box probabilities if not set
        if len(self.boxProb) == 0:
            # Default to single box
            self.boxProb = np.array([1.0], dtype=dp)
    
    def full_move(self, trial_box):
        """
        Corresponds to Delete_FullMove
        Perform a delete move (remove molecule)
        
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
        
        # Check if deletion is allowed
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
        # For delete moves, this is typically just the standard grand canonical acceptance
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
        Corresponds to Delete_Prologue
        Initialize the move before simulation
        """
        # Create temporary arrays
        self.create_temp_array(1)
        
        print(f"(Delete) Move initialized")
    
    def epilogue(self):
        """
        Corresponds to Delete_Epilogue
        Finalize the move after simulation
        """
        print(f"Delete Moves Accepted: {int(self.accpt):15d}")
        print(f"Delete Moves Attempted: {int(self.atmps):15d}")
        
        accpt_rate = self.get_accept_rate()
        print(f"Delete Acceptance Rate: {accpt_rate:15.8f}")
        
        if self.verbose:
            print(f"Delete, Rejections due to overlap: {self.ovlaprej:15d}")
            print(f"Delete, Rejections due to constraint: {self.constrainrej:15d}")
            print(f"Delete, Rejections due to detailed balance: {self.detailedrej:15d}")
    
    def process_io(self, line: str) -> int:
        """
        Corresponds to Delete_ProcessIO
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