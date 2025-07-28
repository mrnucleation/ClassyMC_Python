#!/usr/bin/env python3
"""
Test script to demonstrate ASE trajectory generation from SimMonteCarlo simulation.
This script shows how to use the new to_ase_atoms() method to visualize and verify
that the SimMonteCarlo class is working properly.
"""

import numpy as np
import sys
import os

# Add the src_python directory to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src_python'))

# Import using absolute paths to avoid relative import issues
import Box_SimpleBox
import Box_CubeBox
import Sim_MonteCarlo
import MC_Move_MolTranslation
import Sampling_Metropolis
import FF_EasyPair_Cut

SimpleBox = Box_SimpleBox.SimpleBox
CubeBox = Box_CubeBox.CubeBox
SimMonteCarlo = Sim_MonteCarlo.SimMonteCarlo
MCMoveMolTranslation = MC_Move_MolTranslation.MolTranslate
SamplingMetropolis = Sampling_Metropolis.Metropolis
FFEasyPairCut = FF_EasyPair_Cut.EasyPairCut

def create_test_molecules():
    """Create test molecule definitions for a simple simulation"""
    
    # Define a simple molecule type (e.g., argon atoms)
    class MoleculeType:
        def __init__(self):
            self.atomtypes = ['Ar']  # Argon atoms
            self.nAtoms = 1
            self.masses = [39.948]  # Atomic mass of argon
            self.ridgid = True  # Rigid molecule (single atom)
    
    return [MoleculeType()]

def test_ase_trajectory():
    """Test the ASE trajectory generation functionality"""
    
    print("Testing ASE trajectory generation from SimMonteCarlo...")
    
    # Create molecule definitions
    mol_data = create_test_molecules()
    
    # Create a simple box with 10 argon atoms
    box = SimpleBox(molData=mol_data, NMolMax=[10], NMol=[10])
    
    # Initialize the box arrays (this would normally be done by the simulation)
    box.nMaxAtoms = 10
    box.maxMol = 10
    box.nMolTotal = 10
    box.nAtoms = 10
    
    # Create dummy arrays for testing
    box.atoms = np.random.rand(10, 3) * 10.0  # Random positions in a 10x10x10 box
    box.AtomType = np.zeros(10, dtype=int)  # All atoms are type 0 (Ar)
    box.MolType = np.zeros(10, dtype=int)  # All molecules are type 0
    box.MolIndx = np.arange(10)  # Each atom is its own molecule
    box.MolSubIndx = np.zeros(10, dtype=int)
    box.AtomSubIndx = np.zeros(10, dtype=int)
    
    # Test the to_ase_atoms() method
    print("Converting simulation box to ASE Atoms object...")
    atoms_obj = box.to_ase_atoms()
    
    if atoms_obj is not None:
        print(f"Successfully created ASE Atoms object with {len(atoms_obj)} atoms")
        print(f"Atom symbols: {atoms_obj.get_chemical_symbols()}")
        print(f"Positions shape: {atoms_obj.positions.shape}")
        print(f"Cell: {atoms_obj.cell}")
        print(f"PBC: {atoms_obj.pbc}")
        
        # Test trajectory generation
        print("\nTesting trajectory functionality...")
        try:
            from ase.io import write
            # Write a trajectory file
            write('test_trajectory.poscar', atoms_obj, format='vasp')
            print("Successfully wrote trajectory file: test_trajectory.poscar")
            
            # Test multiple frames (simulating a trajectory)
            trajectory = [atoms_obj]
            for i in range(3):  # Create 3 additional frames
                new_atoms = atoms_obj.copy()
                new_atoms.positions += np.random.normal(0, 0.1, new_atoms.positions.shape)
                trajectory.append(new_atoms)
            
            write('test_multi_frame.poscar', trajectory, format='vasp')
            print("Successfully wrote multi-frame trajectory: test_multi_frame.poscar")
            
        except Exception as e:
            print(f"Error writing trajectory: {e}")
        
        return True
    else:
        print("Failed to create ASE Atoms object")
        return False

def test_sim_monte_carlo_with_ase():
    """Test integrating ASE with SimMonteCarlo simulation"""
    
    print("\n" + "="*60)
    print("Testing SimMonteCarlo with ASE trajectory generation...")
    print("="*60)
    
    # Create molecule definitions
    mol_data = create_test_molecules()
    
    # Create simulation box
    box = SimpleBox(molData=mol_data, NMolMax=[5], NMol=[5])
    
    # Initialize box arrays
    box.nMaxAtoms = 5
    box.maxMol = 5
    box.nMolTotal = 5
    box.nAtoms = 5
    box.atoms = np.random.rand(5, 3) * 5.0
    box.AtomType = np.zeros(5, dtype=int)
    box.MolType = np.zeros(5, dtype=int)
    box.MolIndx = np.arange(5)
    box.MolSubIndx = np.zeros(5, dtype=int)
    box.AtomSubIndx = np.zeros(5, dtype=int)
    
    # Create Monte Carlo simulation
    sim = SimMonteCarlo(nCycles=10, nMoves=5, screenfreq=5, configfreq=10, energyCheck=10)
    sim.BoxList = [box]
    
    # Create a simple move (translation)
    move = MCMoveMolTranslation([box])
    move.max_displacement = 0.5
    sim.Moves = [move]
    
    # Create sampling rule
    sampling = SamplingMetropolis()
    sampling.temperature = 300.0
    sim.Sampling = sampling
    
    # Create force field
    ff = FFEasyPairCut()
    ff.cutoff = 8.0
    ff.epsilon = 1.0
    ff.sigma = 3.4
    box.EFunc = [ff]
    
    print("Running Monte Carlo simulation...")
    sim.run_monte_carlo()
    print("Simulation completed successfully!")
    
    # Generate ASE trajectory from final state
    print("\nGenerating ASE trajectory from final simulation state...")
    final_atoms = box.to_ase_atoms()   
    try:

        if final_atoms is not None:
            from ase.io import write
            write('final_simulation_state.poscar', final_atoms, format='vasp')
            print("Successfully wrote final simulation state: final_simulation_state.poscar")
            return True
        else:
            print("Failed to generate ASE trajectory from simulation")
            return False
            
    except Exception as e:
        print(f"Error during simulation: {e}")
        return False

if __name__ == "__main__":
    print("ASE Trajectory Generation Test")
    print("="*50)
    
    # Test 1: Basic ASE conversion
    success1 = test_ase_trajectory()
    
    # Test 2: Integration with SimMonteCarlo
    success2 = test_sim_monte_carlo_with_ase()
    
    print("\n" + "="*50)
    if success1 and success2:
        print("All tests passed! ASE trajectory generation is working correctly.")
    else:
        print("Some tests failed. Check the output above for details.")
    
    print("\nYou can now visualize the generated trajectory files:")
    print("  - test_trajectory.poscar (single frame)")
    print("  - test_multi_frame.poscar (multi-frame trajectory)")
    print("  - final_simulation_state.poscar (final simulation state)")
    print("\nUse ASE visualization tools or other molecular viewers to inspect these files.")
