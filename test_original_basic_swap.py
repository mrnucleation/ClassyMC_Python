#!/usr/bin/env python3
"""
Test script for the original BasicSwap move class
Tests the swap_in and swap_out methods with the actual implementation
"""

import numpy as np
import sys
import os

from src_python.MC_Move_BasicSwap import BasicSwap
from src_python.CoordinateTypes import Addition, Deletion
from src_python.Box_CubeBox import CubeBox
from src_python.Molecule_Definition import Molecule_Type
from src_python.FF_LJ_Cut import LJ_Cut
from src_python.Script_LoadCoordinates import load_coords
from src_python.Sampling_AcceptAll import AcceptAll

#-----------------------------------------------------------------------------------------------------------
def test_original_with_fixes():
    """Test the original BasicSwap with fixes applied"""
    print("=== Testing Original BasicSwap with Bug Fixes ===")
    mock_lj_data = {
        "atoms": [("Ar", "LJ")],
    }
    atomtypes = ["LJ"]
    LJ_type = Molecule_Type(mock_lj_data, atomtypes=atomtypes)

    box = load_coords("SimpleLJ.clssy", [LJ_type])

    # Define LJ Forcefield and attach to box
    lj_ff = LJ_Cut(nAtomTypes=1)
    lj_ff.rMin = np.zeros(1)
    lj_ff.rMinTable = np.zeros((1,1))
    box.EFunc.append(lj_ff)

    # Create sampling object
    sampling = AcceptAll()

    # Create BasicSwap object
    basic_swap = BasicSwap()

    # Test swap_in
    result = basic_swap.swap_in(box, sampling)
    print(f"Swap in result: {result}")

    # Test swap_out
    result = basic_swap.swap_out(box, sampling)
    print(f"Swap out result: {result}")

#-----------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    test_original_with_fixes()
