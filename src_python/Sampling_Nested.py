"""
Nested Sampling Rule
Corresponds to Sampling_Nested.f90

Implements nested sampling for computing partition functions and free energies.
This is an advanced sampling method that explores phase space systematically
by maintaining a collection of live points.
"""

import math
import numpy as np
import sys
from typing import Optional, List
from .VarPrecision import dp
from .Template_AcceptRule import AcceptRule

def grnd():
    """Generate random number between 0 and 1"""
    return np.random.random()

class Nested(AcceptRule):
    """
    Nested sampling for systematically exploring low-energy configurations.
    Corresponds to the Fortran Nested type.
    """
    
    def __init__(self):
        super().__init__()
        
        # Core nested sampling attributes
        self.parallel = True
        self.firstpass = True
        self.binmiss = 0
        
        # Energy tracking
        self.E_Median = 0.0
        self.E_Min = 0.0  
        self.E_Max = 0.0
        self.dE = 0.0
        self.E_Hist = np.zeros(1001, dtype=dp)  # 0 to 1000 bins
        
        # File and ensemble settings
        self.logfile = "Nested.dat"
        self.ensemble = "canonical"
        
    def make_decision(self, trial_box, e_diff: float, disp, 
                     in_prob: Optional[float] = None, 
                     log_prob: Optional[float] = None, 
                     extra_in: Optional[float] = None) -> bool:
        """
        Corresponds to Nested_MakeDecision
        Make acceptance decision based on nested sampling criterion
        """
        # Check for ensemble restrictions
        if self.ensemble.strip().lower() == 'canonical':
            # Only allow displacement moves in canonical ensemble
            if not all(hasattr(d, 'atmIndx') for d in disp):
                print("In order to use moves that change volume or particle counts", file=sys.stderr)
                print("the flag canonical must be set to false", file=sys.stderr)
                print("Command => modify sampling canonical .false.", file=sys.stderr)
                raise RuntimeError("Ensemble restriction violated")
        
        # Get extra terms for the move
        extra_terms = self.get_extra_terms(disp, trial_box)
        extra_terms_old = self.get_extra_terms_old(disp, trial_box)
        
        # Calculate total energy after the move
        E_Total = trial_box.ETotal + e_diff + extra_terms
        
        # Calculate energy per atom based on move type
        if hasattr(disp[0], 'molType') and hasattr(disp[0], 'addition_type'):
            # Addition move
            n_atoms_added = trial_box.MolData[disp[0].molType]['nAtoms']
            E_PerAtom = E_Total / (trial_box.nAtoms + n_atoms_added)
        elif hasattr(disp[0], 'molType') and hasattr(disp[0], 'deletion_type'):
            # Deletion move
            n_atoms_removed = trial_box.MolData[disp[0].molType]['nAtoms']
            E_PerAtom = E_Total / (trial_box.nAtoms - n_atoms_removed)
        else:
            # Displacement or other move
            E_PerAtom = E_Total / trial_box.nAtoms
        
        # Current energy per atom
        E_PerAtomOld = (trial_box.ETotal + extra_terms_old) / trial_box.nAtoms
        
        # If current system is above median, accept any downhill move
        if E_PerAtomOld > self.E_Median:
            if E_PerAtom < E_PerAtomOld:
                return True
        
        # Check if new energy is below median
        if E_PerAtom <= self.E_Median:
            accept = True
            # Update histogram if within bounds
            if self.E_Min < E_PerAtom < self.E_Max:
                ebin = int(math.floor((E_PerAtom - self.E_Min) * self.dE))
                if 0 <= ebin < len(self.E_Hist):
                    self.E_Hist[ebin] += 1.0
            else:
                self.binmiss += 1
        else:
            # Reject move but still update histogram with old energy
            if self.E_Min < E_PerAtomOld < self.E_Max:
                ebin = int(math.floor((E_PerAtomOld - self.E_Min) * self.dE))
                if 0 <= ebin < len(self.E_Hist):
                    self.E_Hist[ebin] += 1.0
            else:
                self.binmiss += 1
            accept = False
        
        return accept
    
    def make_decision_2box(self, trial_box1, trial_box2, 
                          e_diff1: float, e_diff2: float,
                          disp1, disp2,
                          in_prob: Optional[float] = None,
                          log_prob: Optional[float] = None, 
                          extra_in: Optional[float] = None) -> bool:
        """
        Corresponds to Nested_MakeDecision2Box
        Two-box ensemble support (currently not implemented)
        """
        raise NotImplementedError("Two Box Ensembles are not currently implemented for NestedSampling")
    
    def get_extra_terms(self, disp, trial_box) -> float:
        """
        Corresponds to Nested_GetExtraTerms
        Calculate extra thermodynamic terms based on ensemble
        """
        if self.ensemble.strip().lower() == 'canonical':
            return 0.0
        
        extra_terms = 0.0
        pot_terms = 0.0
        
        if self.ensemble.strip().lower() == 'grand':
            # Grand canonical: add chemical potential terms
            for iType in range(trial_box.nMolTypes):
                pot_terms += trial_box.chempot[iType] * trial_box.NMol[iType]
            
            # Adjust for move type
            if hasattr(disp[0], 'molType') and hasattr(disp[0], 'addition_type'):
                # Addition move
                pot_terms += trial_box.chempot[disp[0].molType]
            elif hasattr(disp[0], 'molType') and hasattr(disp[0], 'deletion_type'):
                # Deletion move
                pot_terms -= trial_box.chempot[disp[0].molType]
            elif hasattr(disp[0], 'newType') and hasattr(disp[0], 'oldType'):
                # Atom exchange move
                mol_new = trial_box.MolIndx[disp[0].newAtmIndx]
                mol_old = trial_box.MolIndx[disp[0].oldAtmIndx]
                type_new = trial_box.MolType[mol_new]
                type_old = trial_box.MolType[mol_old]
                pot_terms += trial_box.chempot[type_new] * trial_box.NMol[type_new]
                pot_terms -= trial_box.chempot[type_old] * trial_box.NMol[type_old]
            elif hasattr(disp[0], 'volNew'):
                # Volume change not allowed in grand ensemble
                raise RuntimeError("Volume change moves are not allowed in Grand Ensemble Mode!")
                
        elif self.ensemble.strip().lower() == 'isobaric':
            # Isobaric: add PV terms
            if hasattr(disp[0], 'volNew'):
                # Volume change move
                extra_terms += disp[0].volNew * trial_box.pressure
            elif hasattr(disp[0], 'atmIndx'):
                # Displacement move
                extra_terms = 0.0
            else:
                raise RuntimeError("Moves besides volume change and translation moves are not allowed in isobaric mode!")
        else:
            raise RuntimeError("Invalid Ensemble was given to the Nested Sampling Algorithm!")
        
        return extra_terms + pot_terms
    
    def get_extra_terms_old(self, disp, trial_box) -> float:
        """
        Corresponds to Nested_GetExtraTermsOld
        Calculate extra terms for the current (old) state
        """
        if self.ensemble.strip().lower() == 'canonical':
            return 0.0
        
        extra_terms = 0.0
        
        if self.ensemble.strip().lower() == 'grand':
            # Grand canonical: current chemical potential terms
            for iType in range(trial_box.nMolTypes):
                extra_terms += trial_box.chempot[iType] * trial_box.NMol[iType]
                
        elif self.ensemble.strip().lower() == 'isobaric':
            # Isobaric: current PV term
            extra_terms += trial_box.volume * trial_box.pressure
        
        return extra_terms
    
    def process_io(self, line: str) -> int:
        """
        Corresponds to Nested_ProcessIO
        Process input commands for nested sampling
        """
        parts = line.strip().split()
        if len(parts) < 4:
            return -1
        
        command = parts[2].lower()
        
        try:
            if command == "adjustfreq":
                value = float(parts[3])
                self.maintFreq = int(math.floor(value))
                
            elif command == "ensemble":
                self.ensemble = parts[3]
                
            elif command == "emax":
                self.E_Max = float(parts[3])
                # Apply energy unit conversion if needed
                # self.E_Max *= outEngUnit
                self.dE = 1000.0 / (self.E_Max - self.E_Min)
                
            elif command == "emin":
                self.E_Min = float(parts[3])
                # Apply energy unit conversion if needed
                # self.E_Min *= outEngUnit
                self.dE = 1000.0 / (self.E_Max - self.E_Min)
                
            elif command == "parallel":
                self.parallel = parts[3].lower() in ['true', '.true.', 't']
                
            elif command == "filename":
                self.logfile = parts[3]
                
            else:
                return -1
                
        except (ValueError, IndexError):
            return -1
        
        return 0
    
    def prologue(self):
        """
        Corresponds to Nested_Prologue
        Initialize nested sampling
        """
        self.E_Median = self.E_Max
        self.E_Hist.fill(0.0)
    
    def maintenance(self):
        """
        Corresponds to Nested_Maintenance
        Update the median energy value and adjust sampling bounds
        """
        if self.firstpass:
            self.E_Hist.fill(0.0)
            self.firstpass = False
            return
        
        print("Updating Nested Sampling...")
        
        # For non-parallel or single processor, work with local histogram
        temp_hist = self.E_Hist.copy()
        
        # Find median energy
        norm = np.sum(temp_hist) * 0.5
        if norm == 0.0:
            print("Zero norm encountered in Nested Sampling")
            print("Simulation is likely Trapped")
            raise RuntimeError("Zero norm in nested sampling")
        
        sum_int = 0.0
        n_median = -1
        for i in range(len(temp_hist)):
            n_median = i
            sum_int += temp_hist[i]
            if sum_int > norm:
                break
        
        # Calculate new median value
        self.E_Median = 0.5 * ((n_median / self.dE) + ((n_median - 1) / self.dE) + 2.0 * self.E_Min)
        self.E_Max = self.E_Median
        self.dE = 1000.0 / (self.E_Max - self.E_Min)
        
        print(f"New Median Value: {self.E_Median}")
        self.E_Hist.fill(0.0) 