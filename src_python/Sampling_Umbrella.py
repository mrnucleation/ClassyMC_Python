"""
Basic Umbrella Sampling Rule
Corresponds to Sampling_Umbrella.f90

Implements basic umbrella sampling with biasing along reaction coordinates.
This is a simplified version that extends UmbrellaWHAM but without the 
complex WHAM analysis.
"""

import numpy as np
import sys
from typing import Optional, List
from VarPrecision import dp
from Sampling_UmbrellaWHAM import UmbrellaWHAM

class Umbrella(UmbrellaWHAM):
    """
    Basic umbrella sampling for biased sampling along reaction coordinates.
    Extends UmbrellaWHAM but provides simpler histogram-based functionality.
    Corresponds to the Fortran Umbrella type.
    """
    
    def __init__(self):
        super().__init__()
        self.histunit = -1  # File unit for histogram output
        
    def prologue(self):
        """
        Corresponds to Umbrella_Prologue
        Initialize umbrella sampling with analysis functions and bias matrix
        """
        # Import here to avoid circular imports
        from AnalysisData import AnalysisArray
        from SimControl import nCycles
        
        # Validate analysis function indices
        for i in range(self.nBiasVar):
            indx = self.AnalysisIndex[i]
            if indx >= len(AnalysisArray) or indx < 0:
                print(f"ERROR! The Umbrella Sampling routine has been directed to an invalid Analysis function")
                print(f"Chosen Function: {indx}")
                print(f"Number of Analysis Functions: {len(AnalysisArray)}")
                raise RuntimeError("Error detected in Umbrella Sampling")
        
        # Validate variable bounds
        for i in range(self.nBiasVar):
            if self.valMin[i] > self.valMax[i]:
                print("ERROR! The given bounds for one of the umbrella variables does not make sense!")
                print("Smallest bin is larger than the largest bin")
                print(f"Minimum Value: {self.valMin[i]}")
                print(f"Maximum Value: {self.valMax[i]}")
                raise RuntimeError("Invalid umbrella variable bounds")
        
        # Set analysis functions to be used per move
        for i in range(self.nBiasVar):
            indx = self.AnalysisIndex[i]
            AnalysisArray[indx].func.usedInMove = True
            AnalysisArray[indx].func.perMove = True
        
        # Calculate index coefficients for N-dimensional to 1D mapping
        # U = a1*x1 + a2*x2 + a3*x3 + ...
        self.indexCoeff[0] = 1
        for i in range(1, self.nBiasVar):
            self.indexCoeff[i] = 1
            for j in range(i):
                self.indexCoeff[i] *= self.nBins[j]
        
        # Calculate total number of umbrella bins
        self.umbrellaLimit = 1
        for i in range(self.nBiasVar):
            self.umbrellaLimit += self.indexCoeff[i] * self.nBins[i]
        
        print(f"Sampling Style: Histogram based Umbrella Sampling")
        print(f"Number of Umbrella Bins: {self.umbrellaLimit}")
        
        # Allocate bias and histogram arrays
        self.UBias = np.zeros(self.umbrellaLimit + 1, dtype=dp)
        self.UHist = np.zeros(self.umbrellaLimit + 1, dtype=dp)
        
        # Get reference bin index
        stat = self.get_u_index_array(self.refVals)
        if stat[1] != 0:
            raise RuntimeError("Error getting reference bin index")
        self.refBin = stat[0]
        
        # Read initial bias values from file
        self.read_initial_bias()
        
        # Setup for WHAM iterations
        self.nWhamItter = int(np.ceil(nCycles / self.maintFreq))
        
        # Open histogram output file
        with open("Umbrella_Histogram.dat", 'w') as f:
            self.histunit = f  # In Python, we'll handle file I/O differently
        
        # Allocate temporary arrays
        self.TempHist = np.zeros(self.umbrellaLimit, dtype=dp)
        self.NewBias = np.zeros(self.umbrellaLimit, dtype=dp)
        
        print(f"Reference values: {self.refVals}, bin: {self.refBin}")
        print(f"Bin Size: {self.UBinSize}")
        
        self.oldIndx = 1
    
    def maintenance(self):
        """
        Corresponds to Umbrella_Maintenance
        Collect and output histogram data
        """
        self.collect_hist()
    
    def collect_hist(self):
        """
        Corresponds to Umbrella_CollectHist
        Collect histogram data and write to file
        """
        print("Collecting Umbrella Histogram...")
        
        # For single processor version, copy histogram
        self.TempHist = self.UHist.copy()
        
        # Write histogram to file
        output_format = "  " + "  ".join([f"{{{j}:10.4f}}" for j in range(self.nBiasVar)]) + "  {hist:22.1f}"
        
        with open("Umbrella_Histogram.dat", 'w') as f:
            for i in range(1, self.umbrellaLimit + 1):
                if self.HistStorage[i] != 0.0:
                    var_values = self.find_var_values(i)
                    # Convert bin indices to actual values
                    actual_values = [var_values[j] * self.UBinSize[j] + self.valMin[j] 
                                   for j in range(self.nBiasVar)]
                    
                    # Format output string
                    format_dict = {str(j): actual_values[j] for j in range(self.nBiasVar)}
                    format_dict['hist'] = self.TempHist[i]
                    
                    line = output_format.format(**format_dict)
                    f.write(line + "\n")
    
    def find_var_values(self, umbrella_index):
        """
        Find the variable values corresponding to a given umbrella index.
        Inverse of the N-dimensional to 1D mapping.
        """
        remaining = umbrella_index - 1
        var_values = np.zeros(self.nBiasVar, dtype=int)
        
        for i in range(self.nBiasVar - 1, -1, -1):
            var_values[i] = remaining // self.indexCoeff[i]
            remaining = remaining % self.indexCoeff[i]
        
        return var_values
    
    def epilogue(self):
        """
        Corresponds to Umbrella_Epilogue
        Final histogram collection and output
        """
        print("Writing Umbrella Sampling Histogram...")
        self.collect_hist() 