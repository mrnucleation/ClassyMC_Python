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
from .Template_MCMove import MCMove
from .CoordinateTypes import Addition, Deletion
from .VarPrecision import dp
from .CommonSampling import sampling

#===============================================================================================================
class BasicSwap(MCMove):
    """
    Basic swap Monte Carlo move for grand canonical simulations.
    Corresponds to the Fortran BasicSwap type.
    """
    #-----------------------------------------------------------------------------------------------------------
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

           
   
    #-----------------------------------------------------------------------------------------------------------
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
    #-----------------------------------------------------------------------------------------------------------
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
        
        # Select molecule type to add using the new random selection function
        type_data = trial_box.pick_random_molecule_type()
        if type_data is None:
            return False
        
        molType = type_data['type_number']
        molData = type_data['type_object']
        
        nAtoms = molData.nAtoms
        nMove = trial_box.nMolTotal  # Next available molecule index
        
        box_dimensions = trial_box.get_dimensions()
        
        low, high = box_dimensions[:,0], box_dimensions[:, 1]
        
        # Generate random position for the insertion point of the new molecule
        newPositions = np.random.uniform(low, high)
        
        #Create new positions for the molecule (unimplemented currently)
        # Temporaily assumes a single atom till proper constructors are integrated
        coords = newPositions.reshape(1, newPositions.shape)
        
        
        # Create coordinates for the new molecule (would need proper placement logic)
        atom_types = molData.atomtypes
        
        # Create Addition displacement
        disp = Addition(molType=molType, 
                        molIndx=nMove, 
                        atomTypes=atom_types, 
                        newPositions=coords)
        
        
        # Check constraints
        if not trial_box.check_constraint(disp):
            self.constrainrej += 1
            return False
        
        # Calculate energy
        e_inter, e_intra, e_diff = trial_box.compute_energy_delta(
            disp, self.tempList, self.tempNnei, computeintra=True
        )
        if e_diff is False:  # energy calculation failed
            self.ovlaprej += 1
            return False
        
        # Check post-energy constraints
        if not trial_box.check_post_energy(disp, e_diff):
            self.constrainrej += 1
            return False
        
        # Calculate acceptance probability
        # For basic swap, this is just the standard grand canonical acceptance
        vol = trial_box.volume
        Prob = vol/(trial_box.nMolTotal+1 )
        #Prob = GasProb*Prob/GenProb       
        # Get extra terms (chemical potential)
        if sampling is not None:
            extra_terms = sampling.get_extra_terms(disp, trial_box)
        else:
            extra_terms = 0.0
        
        # Accept/reject
        if sampling is not None:
            accept = sampling.make_decision(trial_box, e_diff, disp, in_prob=Prob, extra_in=extra_terms)
        else:
            accept = True
        
        if accept:
            self.accpt += 1.0
            trial_box.update_energy(e_diff)
            trial_box.update_position(disp)
        else:
            self.detailedrej += 1
        
        return accept
    #-----------------------------------------------------------------------------------------------
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
        molinfo = trial_box.pick_random_molecule()
        molType = molinfo['mol_type']
        
        if trial_box.NMol[molType] - 1 < trial_box.NMolMin[molType]:
            return False
        
        
        disp = Deletion(molType=molinfo['mol_type'],
                        molIndx=molinfo['mol_index'], 
                        atomTypes=molinfo['atom_types'], 
                        atomIndicies=molinfo['atomindicies']
                        )
        
        
        # Check constraints
        if not trial_box.check_constraint(disp):
            self.constrainrej += 1
            return False
        
        # Calculate energy
        e_inter, e_intra, e_diff, accept = trial_box.compute_energy_delta(
            disp, self.tempList, self.tempNnei, computeintra=True
        )
        if not accept:
            self.ovlaprej += 1
            return accept
        
        # Check post-energy constraints
        if not trial_box.check_post_energy(diso, e_diff):
            self.constrainrej += 1
            return False
        
        # Calculate acceptance probability
        vol = trial_box.volume
        Prob = trial_box.nMolTotal/vol
        
        # Get extra terms (chemical potential)
        if sampling is not None:
            extra_terms = sampling.get_extra_terms(disp, trial_box)
        else:
            extra_terms = 0.0
        
        # Accept/reject
        if sampling is not None:
            accept = sampling.make_decision(trial_box, e_diff, disp, in_prob=Prob, extra_in=extra_terms)
        else:
            accept = True
        
        if accept:
            self.accpt += 1.0
            trial_box.update_energy(e_diff, e_inter, e_intra)
            trial_box.update_position(disp)
        else:
            self.detailedrej += 1
        
        return accept
    #-----------------------------------------------------------------------------------------------------------
    def get_accept_rate(self):
        """
        Corresponds to BasicSwap_GetAcceptRate
        Get the acceptance rate for the move
        """
        return self.accpt / self.atmps
   
    #-----------------------------------------------------------------------------------------------------------
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
    
    #-----------------------------------------------------------------------------------------------------------
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
    
    #-----------------------------------------------------------------------------------------------------------
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
    #-----------------------------------------------------------------------------------------------------------
#===============================================================================================================