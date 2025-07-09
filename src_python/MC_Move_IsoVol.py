"""
Isochoric Volume Monte Carlo Move
Corresponds to Move_MC_IsoVol.f90

Monte Carlo move for the Isobaric ensemble. This move scales the dimensions
of the simulation box uniformly in all directions.
"""

import numpy as np
import sys
import math
from random import random
from Template_MCMove import MCMove
from CoordinateTypes import OrthoVolChange, TriVolChange
from VarPrecision import dp


class IsoVol(MCMove):
    """
    Isochoric volume change Monte Carlo move.
    Corresponds to the Fortran IsoVol type.
    """
    
    def __init__(self):
        super().__init__()
        
        # Move parameters
        self.verbose = True
        self.style = 1  # 1=log scale, 2=linear scale
        self.maxDv = 0.01
        self.tuneMax = True
        self.limit = 3.00
        self.targAccpt = 50.0
        
        # Rejection counters
        self.volrejection = 0
        self.ovlaprej = 0
        self.constrainrej = 0
        self.detailedrej = 0
        
        # Volume change displacement
        self.nDim = 3
        self.disp = [OrthoVolChange()]
        self.disptri = [TriVolChange()]
    
    def constructor(self):
        """
        Corresponds to IsoVol_Constructor
        Initialize the move
        """
        self.create_temp_array(1)
    
    def full_move(self, trial_box):
        """
        Corresponds to IsoVol_FullMove
        Perform an isochoric volume change move
        
        Args:
            trial_box: SimBox instance
            
        Returns:
            bool: Whether the move was accepted
        """
        self.atmps += 1.0
        self.load_box_info(trial_box, self.disp)
        
        # Propose volume change
        if self.style == 1:  # Log scale
            dV = self.maxDv * (2.0 * random() - 1.0)
            self.disp[0].volNew = trial_box.volume * math.exp(dV)
            self.disp[0].volOld = trial_box.volume
        elif self.style == 2:  # Linear scale
            dV = self.maxDv * (2.0 * random() - 1.0)
            self.disp[0].volNew = trial_box.volume + dV
            self.disp[0].volOld = trial_box.volume
        else:
            raise ValueError(f"Invalid volume change style: {self.style}")
        
        # Check for negative volume
        if self.disp[0].volNew < 0.0:
            return False
        
        # Calculate scale factors based on box type
        scale_factor = (self.disp[0].volNew / self.disp[0].volOld) ** (1.0 / 3.0)
        self.disp[0].xScale = scale_factor
        self.disp[0].yScale = scale_factor
        self.disp[0].zScale = scale_factor
        
        # Check constraints
        if not trial_box.check_constraint(self.disp):
            self.constrainrej += 1
            return False
        
        # Calculate energy change
        e_inter, e_intra, e_diff, accept = trial_box.compute_energy_delta(
            self.disp, self.tempList, self.tempNnei, computeintra=False
        )
        if not accept:
            self.ovlaprej += 1
            return False
        
        # Check post-energy constraints
        if not trial_box.check_post_energy(self.disp, e_diff):
            self.constrainrej += 1
            return False
        
        # Calculate probability term
        if self.style == 1:  # Log scale
            prob = (trial_box.nMolTotal + 1) * math.log(self.disp[0].volNew / self.disp[0].volOld)
        else:  # Linear scale
            prob = trial_box.nMolTotal * math.log(self.disp[0].volNew / self.disp[0].volOld)
        
        # Get extra terms (PV term)
        from CommonSampling import sampling
        if sampling is not None:
            extra_terms = sampling.get_extra_terms(self.disp, trial_box)
        else:
            extra_terms = 0.0
        
        # Accept/reject
        if sampling is not None:
            accept = sampling.make_decision(trial_box, e_diff, self.disp, log_prob=prob, extra_in=extra_terms)
        else:
            # Default acceptance if no sampling rule is set
            accept = True
        
        if accept:
            self.accpt += 1.0
            trial_box.update_energy(e_diff, e_inter, e_intra)
            trial_box.update_position(self.disp, self.tempList, self.tempNnei)
        else:
            self.detailedrej += 1
        
        return accept
    
    def maintenance(self):
        """
        Corresponds to IsoVol_Maintenance
        Tune the maximum volume change to achieve target acceptance rate
        """
        if not self.tuneMax:
            return
        
        if self.atmps < 0.5:
            return
        
        if self.get_accept_rate() > self.targAccpt:
            if self.maxDv * 1.01 < self.limit:
                self.maxDv = self.maxDv * 1.01
            else:
                self.maxDv = self.limit
        else:
            self.maxDv = self.maxDv * 0.99
    
    def prologue(self):
        """
        Corresponds to IsoVol_Prologue
        Initialize the move before simulation
        """
        # Initialize box probabilities if not set
        if len(self.boxProb) == 0:
            # Default to single box if BoxArray not available
            self.boxProb = np.array([1.0], dtype=dp)
        
        print(f"(Iso-Volume) Maximum Volume Change: {self.maxDv:15.8f}")
    
    def epilogue(self):
        """
        Corresponds to IsoVol_Epilogue
        Finalize the move after simulation
        """
        print(f"Iso-Volume  Moves Accepted: {int(self.accpt):15d}")
        print(f"Iso-Volume  Moves Attempted: {int(self.atmps):15d}")
        
        accpt_rate = self.get_accept_rate()
        print(f"Iso-Volume Acceptance Rate: {accpt_rate:15.8f}")
        
        if self.tuneMax:
            print(f"Final Maximum Volume Change: {self.maxDv:15.8f}")
        
        if self.verbose:
            print(f"Iso-Volume, Rejections due to overlap: {self.ovlaprej:15d}")
            print(f"Iso-Volume, Rejections due to constraint: {self.constrainrej:15d}")
            print(f"Iso-Volume, Rejections due to detailed balance: {self.detailedrej:15d}")
    
    def process_io(self, line: str) -> int:
        """
        Corresponds to IsoVol_ProcessIO
        Process input commands
        
        Args:
            line: Input line to process
            
        Returns:
            int: Status code (0 for success, negative for error)
        """
        parts = line.strip().split()
        if len(parts) < 5:
            return -1
        
        command = parts[3].lower()
        value_str = parts[4]
        
        try:
            if command == "dynamiclimit":
                self.limit = float(value_str)
            elif command == "dynamictarget":
                self.targAccpt = float(value_str)
            elif command == "maxdv":
                self.maxDv = float(value_str)
            elif command == "style":
                if value_str.lower() == "log":
                    self.style = 1
                elif value_str.lower() == "linear":
                    self.style = 2
                else:
                    return -1
            elif command == "tunemax":
                self.tuneMax = value_str.lower() in ['true', '.true.', 't', '1']
            elif command == "updatefreq":
                self.maintFreq = int(value_str)
            elif command == "verbose":
                self.verbose = value_str.lower() in ['true', '.true.', 't', '1']
            else:
                return -1
        except (ValueError, IndexError):
            return -1
        
        return 0 