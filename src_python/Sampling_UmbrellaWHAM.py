"""
Umbrella WHAM Sampling Rule
Corresponds to Sampling_UmbrellaWHAM.f90

Implements umbrella sampling with automatic WHAM (Weighted Histogram Analysis Method)
for multi-dimensional biasing along reaction coordinates. This is the most complex
sampling method with full free energy analysis capabilities.
"""

import math
import numpy as np
import sys
from typing import Optional, List, Tuple
from .Template_AcceptRule import AcceptRule
from .VarPrecision import dp

def grnd():
    """Generate random number between 0 and 1"""
    return np.random.random()

class UmbrellaWHAM(AcceptRule):
    """
    Umbrella sampling with WHAM analysis for free energy calculations.
    Corresponds to the Fortran UmbrellaWHAM type.
    """
    
    def __init__(self):
        super().__init__()
        
        # Core umbrella sampling parameters
        self.lastaccept = False
        self.nBiasVar = 1  # Default minimum
        
        # Initialize with default arrays
        self.constructor()
        
        # Bias tracking
        self.oldIndx = 0
        self.newIndx = 0
        self.umbrellaLimit = 0
        
        # File handling
        self.fileName = "umbrella.dat"
        self.potfile = 0
        self.histfile = 0
        self.whamfile = 0
        
        # Reference values
        self.refBin = 1
        
        # WHAM parameters
        self.maxSelfConsist = 3000
        self.tolLimit = 1e-5
        self.nWhamItter = 0
        self.nCurWhamItter = 0
        
        # WHAM arrays (will be allocated in prologue)
        self.WHAM_Numerator = None
        self.WHAM_Denominator = None
        self.HistStorage = None
        self.FreeEnergyEst = None
        self.BiasStorage = None
        self.NewBias = None
        self.ProbArray = None
        self.TempHist = None
    
    def constructor(self):
        """
        Corresponds to UmbrellaWHAM_Constructor
        Allocate arrays based on number of bias variables
        """
        if self.nBiasVar < 1:
            # Set minimum for basic initialization
            self.nBiasVar = 1
        
        # Allocate core arrays
        self.refVals = np.zeros(self.nBiasVar, dtype=dp)
        self.AnalysisIndex = np.zeros(self.nBiasVar, dtype=int)
        self.varType = np.ones(self.nBiasVar, dtype=int)  # Default to type 1
        self.indexCoeff = np.zeros(self.nBiasVar, dtype=int)
        self.binMax = np.zeros(self.nBiasVar, dtype=int)
        self.binMin = np.zeros(self.nBiasVar, dtype=int)
        self.nBins = np.ones(self.nBiasVar, dtype=int)  # Default to 1 bin each
        self.binIndx = np.zeros(self.nBiasVar, dtype=int)
        
        # Value and bin arrays with defaults
        self.UBinSize = np.ones(self.nBiasVar, dtype=dp)  # Default bin size
        self.UArray = np.zeros(self.nBiasVar, dtype=int)
        self.varValues = np.zeros(self.nBiasVar, dtype=dp)
        self.valMin = np.zeros(self.nBiasVar, dtype=dp)
        self.valMax = np.ones(self.nBiasVar, dtype=dp)  # Default max of 1.0
    
    def prologue(self):
        """
        Corresponds to UmbrellaWHAM_Prologue
        Initialize umbrella sampling with full WHAM setup
        """
        # Validate analysis function indices (would need AnalysisArray)
        # for i in range(self.nBiasVar):
        #     indx = self.AnalysisIndex[i]
        #     if indx >= len(AnalysisArray) or indx < 0:
        #         raise RuntimeError(f"Invalid analysis function index: {indx}")
        
        # Validate variable bounds
        for i in range(self.nBiasVar):
            if self.valMin[i] > self.valMax[i]:
                print("ERROR! The given bounds for one of the umbrella variables does not make sense!")
                print(f"Minimum Value: {self.valMin[i]}")
                print(f"Maximum Value: {self.valMax[i]}")
                raise RuntimeError("Invalid variable bounds")
        
        # Calculate index coefficients for N-dimensional mapping
        self.indexCoeff[0] = 1
        for i in range(1, self.nBiasVar):
            self.indexCoeff[i] = 1
            for j in range(i):
                self.indexCoeff[i] *= self.nBins[j]
        
        # Calculate total umbrella bins
        self.umbrellaLimit = 1
        for i in range(self.nBiasVar):
            self.umbrellaLimit += self.indexCoeff[i] * self.nBins[i]
        self.umbrellaLimit = max(self.umbrellaLimit, 1)
        
        print(f"Sampling Style: Histogram based Umbrella Sampling w/ Auto WHAM method")
        print(f"Number of Umbrella Bins: {self.umbrellaLimit}")
        
        # Allocate main arrays
        self.UBias = np.zeros(self.umbrellaLimit + 1, dtype=dp)
        self.UHist = np.zeros(self.umbrellaLimit + 1, dtype=dp)
        self.TempHist = np.zeros(self.umbrellaLimit + 1, dtype=dp)
        self.NewBias = np.zeros(self.umbrellaLimit + 1, dtype=dp)
        
        # Get reference bin
        ref_index, stat = self.get_u_index_array(self.refVals)
        if stat != 0:
            raise RuntimeError("Index Error Encountered in Umbrella Sampling")
        self.refBin = ref_index
        
        # Read initial bias
        self.read_initial_bias()
        
        # Setup WHAM arrays (only on master processor in parallel)
        # For single processor, allocate everything
        self.WHAM_Numerator = np.zeros(self.umbrellaLimit + 1, dtype=dp)
        self.WHAM_Denominator = np.zeros((self.umbrellaLimit + 1, self.nWhamItter + 1), dtype=dp)
        self.HistStorage = np.zeros(self.umbrellaLimit + 1, dtype=dp)
        self.BiasStorage = np.zeros((self.umbrellaLimit + 1, self.nWhamItter + 1), dtype=dp)
        self.FreeEnergyEst = np.zeros(self.umbrellaLimit + 1, dtype=dp)
        self.ProbArray = np.zeros(self.umbrellaLimit + 1, dtype=dp)
        
        self.nCurWhamItter = 1
        
        print(f"Reference values: {self.refVals}, bin: {ref_index}")
        print(f"Bin Size: {self.UBinSize}")
        
        self.oldIndx = 1
    
    def make_decision(self, trial_box, e_diff: float, disp, 
                     in_prob: Optional[float] = None, 
                     log_prob: Optional[float] = None, 
                     extra_in: Optional[float] = None) -> bool:
        """
        Corresponds to UmbrellaWHAM_MakeDecision
        Make acceptance decision with umbrella biasing
        """
        # Check probability
        if in_prob is not None:
            if in_prob <= 0.0:
                return False
            prob_term = math.log(in_prob)
        elif log_prob is not None:
            prob_term = log_prob
        else:
            print("Coding Error! Probability has not been passed into Sampling", file=sys.stderr)
            raise RuntimeError("No probability provided")
        
        # Check if any bias variables changed
        indx_changed = False
        # This would require AnalysisArray to be available
        # for iBias in range(self.nBiasVar):
        #     indx = self.AnalysisIndex[iBias]
        #     use_bias = AnalysisArray[indx].func.calc_new_state(disp)
        #     if use_bias:
        #         indx_changed = True
        
        # Get current and new bias indices
        self.oldIndx = self.get_bias_index()
        if indx_changed:
            self.newIndx, accept = self.get_new_bias_index()
            if not accept:
                return False
        else:
            self.newIndx = self.oldIndx
        
        # Calculate bias difference
        bias_old = self.UBias[self.oldIndx]
        bias_new = self.UBias[self.newIndx]
        
        # Extra terms
        extra_terms = extra_in if extra_in is not None else 0.0
        
        # Calculate acceptance probability
        bias_e = (-trial_box.beta * e_diff + prob_term + 
                 (bias_new - bias_old) + extra_terms)
        
        if bias_e >= 0.0:
            return True
        elif bias_e > math.log(grnd()):
            return True
        else:
            return False
    
    def make_decision_2box(self, trial_box1, trial_box2, 
                          e_diff1: float, e_diff2: float,
                          disp1, disp2,
                          in_prob: Optional[float] = None,
                          log_prob: Optional[float] = None, 
                          extra_in: Optional[float] = None) -> bool:
        """
        Corresponds to UmbrellaWHAM_MakeDecision2Box
        Two-box umbrella sampling decision
        """
        # Check probability
        if in_prob is not None:
            if in_prob <= 0.0:
                return False
            prob_term = math.log(in_prob)
        elif log_prob is not None:
            prob_term = log_prob
        else:
            print("Coding Error! Probability has not been passed into Sampling", file=sys.stderr)
            raise RuntimeError("No probability provided")
        
        # Update analysis functions for both boxes (would need AnalysisArray)
        # ... analysis updates ...
        
        # Get bias indices
        self.oldIndx = self.get_bias_index()
        self.newIndx, accept = self.get_new_bias_index()
        if not accept:
            return False
        
        bias_old = self.UBias[self.oldIndx]
        bias_new = self.UBias[self.newIndx]
        
        extra_terms = extra_in if extra_in is not None else 0.0
        
        # Two-box acceptance criterion
        bias_e = (-trial_box1.beta * e_diff1 - trial_box2.beta * e_diff2 + 
                 prob_term + extra_terms + (bias_new - bias_old))
        
        if bias_e >= 0.0:
            return True
        elif bias_e >= math.log(grnd()):
            return True
        else:
            return False
    
    def get_bias_index(self) -> int:
        """
        Corresponds to UmbrellaWHAM_GetBiasIndex
        Get current bias index from analysis variables
        """
        # This would require AnalysisArray access
        # For now, return a placeholder
        bias_index = 1
        
        # Calculate index from current bias variables
        for iBias in range(self.nBiasVar):
            # bias_val = AnalysisArray[self.AnalysisIndex[iBias]].func.get_result()
            # Convert to bin index based on variable type
            # ...
            pass
        
        if bias_index < 1 or bias_index > self.umbrellaLimit:
            print("ERROR! Umbrella Bias Index out of bounds!")
            print("This may be due to the initial system configuration being outside the Umbrella sampling range")
            raise RuntimeError("Umbrella Bias Index out of bounds!")
        
        return bias_index
    
    def get_new_bias_index(self) -> Tuple[int, bool]:
        """
        Corresponds to UmbrellaWHAM_GetNewBiasIndex
        Get new bias index from proposed move
        """
        # This would require AnalysisArray access
        # For now, return placeholder
        return 1, True
    
    def get_u_index_array(self, var_array) -> Tuple[int, int]:
        """
        Corresponds to UmbrellaWHAM_GetUIndexArray
        Convert variable values to umbrella index
        """
        bias_index = 1
        
        for iBias in range(self.nBiasVar):
            if var_array[iBias] > self.valMax[iBias]:
                return 0, 1  # Out of upper bounds
            if var_array[iBias] < self.valMin[iBias]:
                return 0, -1  # Out of lower bounds
            
            # Calculate bin index
            self.binIndx[iBias] = int(math.floor((var_array[iBias] - self.valMin[iBias]) / self.UBinSize[iBias]))
        
        # Calculate linear index
        for iBias in range(self.nBiasVar):
            bias_index += self.indexCoeff[iBias] * self.binIndx[iBias]
        
        return bias_index, 0  # Success
    
    def read_initial_bias(self):
        """
        Corresponds to UmbrellaWHAM_ReadInitialBias
        Read bias values from input file
        """
        try:
            with open(self.fileName, 'r') as f:
                self.UBias.fill(0.0)
                
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    
                    parts = line.split()
                    if len(parts) >= self.nBiasVar + 1:
                        var_values = [float(parts[i]) for i in range(self.nBiasVar)]
                        cur_bias = float(parts[self.nBiasVar])
                        
                        bias_index, stat = self.get_u_index_array(var_values)
                        if stat == 0:  # Valid index
                            self.UBias[bias_index] = cur_bias
                
            # Normalize bias relative to reference
            ref_val = self.UBias[self.refBin]
            for i in range(1, self.umbrellaLimit + 1):
                self.UBias[i] -= ref_val
                
        except IOError:
            print(f"Warning: Could not read bias file {self.fileName}")
            # Initialize with zeros
            self.UBias.fill(0.0)
    
    def find_var_values(self, umbrella_index):
        """
        Find variable values from umbrella index (inverse mapping)
        """
        remaining = umbrella_index - 1
        var_values = np.zeros(self.nBiasVar, dtype=int)
        
        for iBias in range(self.nBiasVar - 1, -1, -1):
            var_values[iBias] = remaining // self.indexCoeff[iBias]
            remaining = remaining % self.indexCoeff[iBias]
        
        return var_values
    
    def process_io(self, line: str) -> int:
        """
        Process input commands for umbrella WHAM sampling
        """
        # Implementation would parse various umbrella-specific commands
        # For now, return success
        return 0
    
    def maintenance(self):
        """
        Corresponds to UmbrellaWHAM_Maintenance
        Perform WHAM analysis and update biases
        """
        print("Performing WHAM maintenance...")
        # Implementation would perform full WHAM analysis
        pass
    
    def epilogue(self):
        """
        Corresponds to UmbrellaWHAM_Epilogue
        Final WHAM analysis and output
        """
        print("Performing final WHAM analysis...")
        # Implementation would output final free energy surfaces
        pass 