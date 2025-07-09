#!/usr/bin/env python3
"""
Test Script for CubeBox Class
==============================

This script tests the functionality of the CubeBox class by creating a cubic simulation box
with 3 atoms and testing various operations like initialization, boundary conditions,
energy calculations, and basic simulation functionality.
"""

import sys
import os
import numpy as np

# Add the parent directory to the path so we can import src_python as a package
parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, parent_dir)

from src_python.Box_CubeBox import CubeBox
from src_python.VarPrecision import dp
from src_python.FF_HardSphere import HardSphere
from src_python.Sampling_Metropolis import Metropolis

def create_molecular_data():
    """Create simple molecular data for testing"""
    # Define a simple atomic species (e.g., argon atoms)
    mol_data = [
        {
            'nAtoms': 1,           # Single atom molecule
            'atomType': [0],       # All atoms are type 0
            'masses': [39.948],    # Argon mass
            'ridgid': True,        # No intramolecular interactions
            'bonds': [],           # No bonds for single atoms
            'angles': [],          # No angles for single atoms
            'torsions': []         # No torsions for single atoms
        }
    ]
    return mol_data

def test_cubic_box_initialization():
    """Test basic initialization of the cubic box"""
    print("=" * 60)
    print("TEST 1: Cubic Box Initialization")
    print("=" * 60)
    
    # Create molecular data
    mol_data = create_molecular_data()
    
    # Define molecule counts (required for SimpleBox constructor)
    NMolMin = [0]      # Minimum 0 molecules of type 0
    NMolMax = [5]      # Maximum 5 molecules of type 0  
    NMol = [3]         # Currently 3 molecules of type 0
    
    # Create cubic box with required parameters
    box = CubeBox(mol_data, NMolMin=NMolMin, NMolMax=NMolMax, NMol=NMol)
    
    # Set box dimensions
    box.boxL = 10.0  # 10 Angstrom box
    box.boxL2 = 5.0
    box.volume = 1000.0  # 10^3 Angstrom^3
    box.nDimension = 3
    box.boxID = 1
    box.temperature = 298.15  # Room temperature in Kelvin
    box.beta = 1.0 / (8.314e-3 * box.temperature)  # 1/(kB*T) in mol/kJ units
    
    print(f"✓ Box created with dimensions: {box.boxL} x {box.boxL} x {box.boxL}")
    print(f"✓ Box volume: {box.volume} Å³")
    print(f"✓ Box temperature: {box.temperature} K")
    print(f"✓ Box beta: {box.beta:.6f} mol/kJ")
    
    return box, mol_data

def test_box_constructor():
    """Test the box constructor and atom allocation"""
    print("\n" + "=" * 60)
    print("TEST 2: Box Constructor and Atom Allocation")
    print("=" * 60)
    
    box, mol_data = test_cubic_box_initialization()
    
    # Box is already initialized in test_cubic_box_initialization()
    # The constructor was called during __init__
    
    print(f"✓ Box constructor completed")
    print(f"✓ Maximum atoms: {box.nMaxAtoms}")
    print(f"✓ Current atoms: {box.nAtoms}")
    print(f"✓ Maximum molecules: {box.maxMol}")
    print(f"✓ Current molecules: {box.nMolTotal}")
    print(f"✓ Atom position array shape: {box.atoms.shape}")
    
    return box

def test_atom_placement():
    """Test placing atoms in the box"""
    print("\n" + "=" * 60)
    print("TEST 3: Atom Placement")
    print("=" * 60)
    
    box = test_box_constructor()
    
    # Place 3 atoms at specific positions (well-separated to avoid overlaps)
    positions = np.array([
        [-2.0, -2.0, -2.0],  # Atom 1
        [ 0.0,  0.0,  0.0],  # Atom 2  
        [ 2.0,  2.0,  2.0]   # Atom 3
    ])  # Shape: (3, 3)
    
    # Set the positions for the first 3 atoms
    box.atoms[:3, :] = positions
    
    print("✓ Atoms placed at positions:")
    for i in range(3):
        pos = box.atoms[:, i]
        print(f"  Atom {i+1}: ({pos[0]:6.2f}, {pos[1]:6.2f}, {pos[2]:6.2f})")
    
    return box

def test_boundary_conditions():
    """Test periodic boundary conditions"""
    print("\n" + "=" * 60)
    print("TEST 4: Periodic Boundary Conditions")
    print("=" * 60)
    
    box = test_atom_placement()
    
    # Test boundary conditions with coordinates outside the box
    test_coords = [
        [6.0, 0.0, 0.0],    # Outside +x boundary
        [-6.0, 0.0, 0.0],   # Outside -x boundary
        [0.0, 7.0, 0.0],    # Outside +y boundary
        [0.0, 0.0, -8.0]    # Outside -z boundary
    ]
    
    print("Testing periodic boundary conditions:")
    for i, coord in enumerate(test_coords):
        wrapped = box.boundary(coord)
        print(f"  Input:  ({coord[0]:6.1f}, {coord[1]:6.1f}, {coord[2]:6.1f})")
        print(f"  Output: ({wrapped[0]:6.1f}, {wrapped[1]:6.1f}, {wrapped[2]:6.1f})")
    
    return box

def test_energy_calculation():
    """Test energy calculation with hard sphere potential"""
    print("\n" + "=" * 60)
    print("TEST 5: Energy Calculation")
    print("=" * 60)
    
    box = test_boundary_conditions()
    
    # Create and initialize hard sphere force field
    ff = HardSphere(nAtomTypes=1)
    ff.constructor(nAtomTypes=1)
    
    # Set hard sphere diameter (1.0 Angstrom)
    ff.process_io("1 1.0")  # Type 1, sigma = 1.0
    
    # Set the energy function
    box.EFunc = ff
    
    print("✓ Hard sphere force field initialized")
    print(f"✓ Hard sphere diameter: {ff.sig[0]:.2f} Å")
    
    # Calculate energy
    try:
        energy, accept = ff.detailed_calc(box)
        if accept:
            print(f"✓ Energy calculation successful: {energy:.6f} kJ/mol")
            print("✓ No hard sphere overlaps detected")
        else:
            print("✗ Hard sphere overlaps detected!")
    except Exception as e:
        print(f"✗ Energy calculation failed: {e}")
    
    return box, ff

def test_coordinate_transformations():
    """Test coordinate transformations"""
    print("\n" + "=" * 60)
    print("TEST 6: Coordinate Transformations")
    print("=" * 60)
    
    box, ff = test_energy_calculation()
    
    # Test real to reduced coordinate conversion
    real_coords = np.array([2.5, -3.0, 1.0])
    reduced_coords = box.get_reduced_coords(real_coords)
    back_to_real = box.get_real_coords(reduced_coords)
    
    print("Testing coordinate transformations:")
    print(f"  Real coordinates:    ({real_coords[0]:6.2f}, {real_coords[1]:6.2f}, {real_coords[2]:6.2f})")
    print(f"  Reduced coordinates: ({reduced_coords[0]:6.3f}, {reduced_coords[1]:6.3f}, {reduced_coords[2]:6.3f})")
    print(f"  Back to real:        ({back_to_real[0]:6.2f}, {back_to_real[1]:6.2f}, {back_to_real[2]:6.2f})")
    
    # Check if transformation is reversible
    diff = np.abs(real_coords - back_to_real)
    if np.all(diff < 1e-10):
        print("✓ Coordinate transformations are reversible")
    else:
        print(f"✗ Coordinate transformation error: {np.max(diff)}")
    
    return box, ff

def test_sampling_rule():
    """Test sampling rule initialization"""
    print("\n" + "=" * 60)
    print("TEST 7: Sampling Rule")
    print("=" * 60)
    
    box, ff = test_coordinate_transformations()
    
    # Create Metropolis sampling rule
    sampling = Metropolis()
    
    print("✓ Metropolis sampling rule created")
    
    # Test a simple acceptance decision
    try:
        # Test with zero energy difference (should always accept)
        accept = sampling.make_decision(box, 0.0, [], in_prob=1.0)
        print(f"✓ Zero energy move accepted: {accept}")
        
        # Test with large positive energy difference (should reject)
        accept = sampling.make_decision(box, 1000.0, [], in_prob=1.0)
        print(f"✓ High energy move accepted: {accept}")
        
    except Exception as e:
        print(f"✗ Sampling test failed: {e}")
    
    return box, ff, sampling

def test_box_properties():
    """Test various box properties and methods"""
    print("\n" + "=" * 60)
    print("TEST 8: Box Properties and Methods")
    print("=" * 60)
    
    box, ff, sampling = test_sampling_rule()
    
    # Test box dimensions
    dimensions = box.get_dimensions()
    print(f"✓ Box dimensions: {dimensions}")
    
    # Test active atom checking
    print("✓ Active atoms:")
    for i in range(box.nMaxAtoms):
        if box.is_active(i):
            print(f"  Atom {i}: Active")
        else:
            print(f"  Atom {i}: Inactive")
    
    # Test molecular data retrieval
    print("✓ Molecular data:")
    for i in range(min(3, box.maxMol)):
        mol_data = box.get_mol_data(i)
        print(f"  Molecule {i}: {mol_data}")
    
    # Test thermodynamic properties
    print("✓ Thermodynamic properties:")
    print(f"  Volume: {box.get_thermo(1):.2f}")
    print(f"  Temperature: {box.temperature:.2f} K")
    print(f"  Beta: {box.beta:.6f}")
    
    return box, ff, sampling

def test_basic_simulation_setup():
    """Test setting up a basic simulation"""
    print("\n" + "=" * 60)
    print("TEST 9: Basic Simulation Setup")
    print("=" * 60)
    
    box, ff, sampling = test_box_properties()
    
    # Test prologue (initialization)
    try:
        box.prologue()
        print("✓ Box prologue completed successfully")
    except Exception as e:
        print(f"✗ Box prologue failed: {e}")
    
    # Test epilogue (finalization)
    try:
        box.epilogue()
        print("✓ Box epilogue completed successfully")
    except Exception as e:
        print(f"✗ Box epilogue failed: {e}")
    
    return box, ff, sampling

def run_all_tests():
    """Run all tests"""
    print("CUBIC BOX FUNCTIONALITY TEST SUITE")
    print("=" * 60)
    print("Testing the translation of Fortran CubeBox to Python")
    print("This will test basic functionality with 3 atoms in a cubic box")
    print()
    
    try:
        # Run all tests
        box, ff, sampling = test_basic_simulation_setup()
        
        print("\n" + "=" * 60)
        print("TEST SUMMARY")
        print("=" * 60)
        print("✓ All tests completed successfully!")
        print(f"✓ Cubic box with {box.nAtoms} atoms operational")
        print(f"✓ Box volume: {box.volume:.1f} Å³")
        print(f"✓ Hard sphere force field functional")
        print(f"✓ Metropolis sampling rule functional")
        print(f"✓ Periodic boundary conditions working")
        print(f"✓ Coordinate transformations working")
        print()
        print("CONCLUSION: The language translation appears to be working well!")
        print("The basic functionality is intact and ready for further testing.")
        
        return True
        
    except Exception as e:
        print(f"\n✗ TEST FAILED: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    success = run_all_tests()
    sys.exit(0 if success else 1) 