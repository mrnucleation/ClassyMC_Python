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
from src_python.Box_SimpleBox import SimpleBox
from src_python.CoordinateTypes import Displacement
from src_python.Molecule_Definition import Molecule_Type

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

def test_displacement_energy_consistency():
    """Test that delta energy calculations match full energy recalculations"""
    print("\nTesting displacement energy consistency...")
    
    # Set random seed for reproducibility
    np.random.seed(42)
    
    # Create a simple molecular system with LJ atoms
    mock_lj_data = {
        "atoms": [("Ar", "LJ")],
    }
    atomtypes = ["LJ"]
    LJ_type = Molecule_Type(mock_lj_data, atomtypes=atomtypes)
    
    # Create simulation box with 10 LJ atoms
    n_atoms = 10
    box_length = 10.0
    
    # Create box with proper molData
    box = SimpleBox(molData=[LJ_type], NMol=[n_atoms], nDimensions=3)
    box.boxID = 1
    box.boxL = np.array([box_length, box_length, box_length], dtype=dp)
    
    # Generate random positions for atoms
    box.atoms = np.random.uniform(0, box_length, size=(n_atoms, 3)).astype(dp)
    
    # Set up atom indexing
    box.MolIndx = np.arange(n_atoms, dtype=int)
    box.MolSubIndx = np.zeros(n_atoms, dtype=int)
    box.AtomType = np.zeros(n_atoms, dtype=int)
    box.AtomSubIndx = np.zeros(n_atoms, dtype=int)
    
    # Initialize force field
    lj_ff = LJ_Cut(nAtomTypes=1)
    lj_ff.constructor()
    lj_ff.process_io("1 1.0 1.0 0.0")  # epsilon=1.0, sigma=1.0, rmin=0.0
    lj_ff.process_io("rcut 5.0")
    lj_ff.rMin = np.zeros(1)
    lj_ff.rMinTable = np.zeros((1, 1))
    
    box.EFunc = [lj_ff]
    
    # Compute initial energy
    success = box.compute_energy()
    assert success, "Initial energy computation failed"
    print(f"  Initial total energy: {box.ETotal:.6f}")
    
    # Perform multiple random displacements and test energy consistency
    n_moves = 300
    max_displacement = 0.5
    energy_tolerance = 1e-10
    
    max_error = 0.0
    errors = []
    
    print(f"  Performing {n_moves} random displacement moves...")
    
    for move_idx in range(n_moves):
        # Select random atom to move
        atom_idx = np.random.randint(0, n_atoms)
        mol_idx = atom_idx  # Since each molecule has 1 atom
        mol_type = 0
        
        # Generate random displacement
        delta = np.random.uniform(-max_displacement, max_displacement, size=3).astype(dp)
        
        # Store original position
        original_pos = box.atoms[atom_idx].copy()
        
        # Create new position (with periodic boundary conditions)
        new_pos = original_pos + delta
        new_pos = new_pos % box_length  # Apply PBC
        new_pos = new_pos.reshape(1, -1)
        
        # Store energy before move
        energy_before = box.ETotal
        
        # Create displacement object
        disp = Displacement(mol_type, mol_idx, np.array([atom_idx]), new_pos)
        
        # Compute delta energy using incremental method
        E_Inter_delta, E_Intra_delta, accept = box.compute_energy_delta(disp)
        assert accept, f"Energy delta computation failed for move {move_idx}"
        
        # Apply displacement
        box.update_position(disp)
        
        # Compute full energy from scratch
        success = box.compute_energy()
        assert success, f"Full energy computation failed for move {move_idx}"
        
        energy_after_full = box.ETotal
        
        # Calculate expected energy using delta
        energy_after_predicted = energy_before + E_Inter_delta
        
        # Compare energies
        energy_error = abs(energy_after_full - energy_after_predicted)
        errors.append(energy_error)
        max_error = max(max_error, energy_error)
        
        # Check if within tolerance
        if energy_error > energy_tolerance:
            print(f"    WARNING: Move {move_idx} energy mismatch!")
            print(f"      Atom {atom_idx}: {original_pos} -> {new_pos.flatten()}")
            print(f"      Energy before: {energy_before:.10f}")
            print(f"      Delta energy: {E_Inter_delta:.10f}")
            print(f"      Predicted after: {energy_after_predicted:.10f}")
            print(f"      Actual after: {energy_after_full:.10f}")
            print(f"      Error: {energy_error:.2e}")
            
        # Revert the move for next test (move atom back)
        revert_disp = Displacement(mol_type, mol_idx, np.array([atom_idx]), 
                                   original_pos.reshape(1, -1))
        box.update_position(revert_disp)
        
        # Verify we're back to original energy
        success = box.compute_energy()
        assert success, f"Energy computation failed after reverting move {move_idx}"
        revert_error = abs(box.ETotal - energy_before)
        
        if revert_error > energy_tolerance:
            print(f"    WARNING: Move {move_idx} revert energy mismatch!")
            print(f"      Expected: {energy_before:.10f}")
            print(f"      Got: {box.ETotal:.10f}")
            print(f"      Error: {revert_error:.2e}")
    
    # Print statistics
    errors = np.array(errors)
    print(f"\n  Statistics for {n_moves} moves:")
    print(f"    Maximum error: {max_error:.2e}")
    print(f"    Mean error: {np.mean(errors):.2e}")
    print(f"    Median error: {np.median(errors):.2e}")
    print(f"    Std dev error: {np.std(errors):.2e}")
    
    # Assert that all moves had acceptable error
    assert max_error < 1e-8, f"Maximum energy error {max_error:.2e} exceeds tolerance"
    
    print("\n  Displacement energy consistency test completed successfully!")
    print(f"  All {n_moves} moves had delta energy matching full recalculation within tolerance.")

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
        test_displacement_energy_consistency()
        test_tail_corrections()
        print("\nAll tests passed!")
    except Exception as e:
        print(f"Test failed with error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1) 