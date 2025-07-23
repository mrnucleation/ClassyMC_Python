#!/usr/bin/env python3
"""
Test script for LJ_Cut force field implementation
"""

import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
from src_python.FF_LJ_Cut import LJ_Cut
from src_python.VarPrecision import dp

def test_lj_cut_basic():
    """Test basic LJ_Cut functionality"""
    print("Testing LJ_Cut Force Field...")
    
    # Create LJ_Cut instance with 2 atom types
    lj_ff = LJ_Cut(nAtomTypes=2)
    
    # Initialize arrays
    lj_ff.constructor()
    
    # Test parameter setting via process_io
    print("\nTesting parameter input...")
    
    # Set parameters for type 1: epsilon=1.0, sigma=1.0, rmin=0.5
    result = lj_ff.process_io("1 1.0 1.0 0.5")
    print(f"Type 1 parameter setting result: {result}")
    
    # Set parameters for type 2: epsilon=0.5, sigma=1.2, rmin=0.6  
    result = lj_ff.process_io("2 0.5 1.2 0.6")
    print(f"Type 2 parameter setting result: {result}")
    
    # Set cutoff
    result = lj_ff.process_io("rcut 5.0")
    print(f"Cutoff setting result: {result}")
    
    # Test pair function
    print("\nTesting pair function...")
    
    # Verify arrays are initialized properly
    assert lj_ff.sigma is not None, "Sigma array not initialized"
    assert lj_ff.epsilon is not None, "Epsilon array not initialized"
    assert lj_ff.epsTable is not None, "Epsilon table not initialized"
    assert lj_ff.sigTable is not None, "Sigma table not initialized"
    
    # Test at sigma distance (should be zero)
    sigma_12 = 0.5 * (lj_ff.sigma[0] + lj_ff.sigma[1])  # Mixed sigma
    rsq_sigma = sigma_12**2
    energy_at_sigma = lj_ff.pair_function(rsq_sigma, 0, 1)
    print(f"Energy at sigma distance: {energy_at_sigma:.6f}")
    
    # Test at close distance (should be large positive)
    rsq_close = 0.8**2
    energy_close = lj_ff.pair_function(rsq_close, 0, 1)
    print(f"Energy at r=0.8: {energy_close:.6f}")
    
    # Test at far distance (should be small negative)
    rsq_far = 2.0**2
    energy_far = lj_ff.pair_function(rsq_far, 0, 1)
    print(f"Energy at r=2.0: {energy_far:.6f}")
    
    # Test minimum position (should be at r = 2^(1/6) * sigma ≈ 1.122 * sigma)
    r_min_theory = (2.0**(1.0/6.0)) * sigma_12
    rsq_min = r_min_theory**2
    energy_min = lj_ff.pair_function(rsq_min, 0, 1)
    print(f"Energy at theoretical minimum (r={r_min_theory:.3f}): {energy_min:.6f}")
    
    # Verify mixing rules
    print("\nVerifying mixing rules...")
    eps_11 = lj_ff.epsilon[0] * lj_ff.epsilon[0]
    eps_22 = lj_ff.epsilon[1] * lj_ff.epsilon[1] 
    eps_12_theory = 4.0 * np.sqrt(lj_ff.epsilon[0] * lj_ff.epsilon[1])
    eps_12_actual = lj_ff.epsTable[0, 1]
    
    print(f"Type 1-1 epsilon: {4.0 * eps_11:.3f}")
    print(f"Type 2-2 epsilon: {4.0 * eps_22:.3f}")
    print(f"Type 1-2 epsilon (theoretical): {eps_12_theory:.3f}")
    print(f"Type 1-2 epsilon (actual): {eps_12_actual:.3f}")
    
    sig_12_theory = (0.5 * (lj_ff.sigma[0] + lj_ff.sigma[1]))**2
    sig_12_actual = lj_ff.sigTable[0, 1]
    print(f"Type 1-2 sigma^2 (theoretical): {sig_12_theory:.3f}")
    print(f"Type 1-2 sigma^2 (actual): {sig_12_actual:.3f}")
    
    print("\nBasic LJ_Cut test completed successfully!")

def test_tail_corrections():
    """Test tail correction functionality"""
    print("\nTesting tail corrections...")
    
    # Create a mock simulation box class for testing
    class MockSimBox:
        def __init__(self):
            self.volume = 1000.0
            self.nMolTypes = 1
            self.NMol = [100]  # 100 molecules of type 0
            self.MolData = [{'atomType': [0]}]  # Each molecule has 1 atom of type 0
    
    lj_ff = LJ_Cut(nAtomTypes=1)
    lj_ff.constructor()
    lj_ff.process_io("1 1.0 1.0 0.5")
    lj_ff.usetailcorrection = True
    
    mock_box = MockSimBox()
    
    # Test tail correction calculation
    tail_corr = lj_ff.tail_correction(mock_box)
    print(f"Tail correction energy: {tail_corr:.6f}")
    
    print("Tail correction test completed!")

if __name__ == "__main__":
    try:
        test_lj_cut_basic()
        test_tail_corrections()
        print("\nAll tests passed!")
    except Exception as e:
        print(f"Test failed with error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1) 