#!/usr/bin/env python3
"""
Monte Carlo simulation script for 100 LJ atoms in vapor phase using NNTranslate move.

Simulation Conditions:
- 100 LJ atoms in a large cubic box (vapor phase)
- Lennard-Jones potential with cutoff = 2.5 sigma
- Temperature = 0.7 (reduced units, kT/epsilon = 0.7)
- Using neural network-guided translation moves (NNTranslate)
- Large box size to simulate vapor-like low density

This script is designed for cross-validation with LAMMPS simulations.
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import time
from src_python.Molecule_Definition import Molecule_Type
from src_python.Box_CubeBox import CubeBox
from src_python.FF_LJ_Cut import LJ_Cut
from src_python.MC_Move_NNTranslate import NNTranslate
from src_python.Sampling_Metropolis import Metropolis
from src_python.Sim_MonteCarlo import SimMonteCarlo
from src_python.MC_Move_MolTranslation import MolTranslate

from math import sqrt

# Set random seed for reproducibility
np.random.seed(12345)

def create_vapor_box():
    """
    Create a large cubic box with 100 LJ atoms at low density (vapor phase).
    """
    print("Creating vapor phase box with 100 LJ atoms...")
    
    natoms = 10
    
    # Define molecule topology (single LJ atom per molecule)
    mock_lj_data = {
        "atoms": [("Ar", "LJ")],  # Single LJ atom type
    }
    atomtypes = ["LJ"]
    LJ_type = Molecule_Type(mock_lj_data, atomtypes=atomtypes)
    
    # Set up box with 100 molecules
    NMol = [natoms]      # 100 LJ atoms
    NMolMin = [0]     # Minimum 0 (for GCMC, not used here)
    NMolMax = [natoms]   # Maximum 150 (allows for fluctuations)
    
    # Create cubic box
    box = CubeBox([LJ_type], NMolMin=NMolMin, NMolMax=NMolMax, NMol=NMol)
    
    # Set large box dimensions for vapor phase
    # For vapor phase, we want low density: rho* = rho * sigma^3 ~ 0.01-0.1
    # With 100 atoms and rho* = 0.05, box volume should be ~2000 sigma^3
    # So box length = (2000)^(1/3) ~ 12.6 sigma
    box_length = 8.0  # 15 sigma units (large vapor box)
    
    box.load_dimension([box_length])
    box.boxID = 1
    box.nDimension = 3
    
    # Set molecule count arrays
    box.NMol = np.array(NMol)
    box.NMolMin = np.array(NMolMin)
    box.NMolMax = np.array(NMolMax)
    
    # Initialize atom counts and arrays manually (following test pattern)
    box.nAtoms = sum(NMol)  # 100 atoms total
    box.nMaxAtoms = sum(NMolMax)  # 150 max atoms
    box.maxMol = len(NMol)  # 1 molecule type
    box.nMolTotal = sum(NMol)  # 100 molecules total
    
    # Allocate arrays
    box.atoms = np.zeros((box.nMaxAtoms, 3), dtype=np.float64)
    box.AtomType = np.zeros(box.nMaxAtoms, dtype=int)
    box.MolType = np.zeros(box.nMaxAtoms, dtype=int)
    box.MolIndx = np.zeros(box.nMaxAtoms, dtype=int)
    box.MolSubIndx = np.zeros(box.nMaxAtoms, dtype=int)
    box.AtomSubIndx = np.zeros(box.nMaxAtoms, dtype=int)
    #box.centerMass = np.zeros((box.nMolTotal, 3), dtype=np.float64)  # Center of mass for each molecule
    
    # Generate random initial positions for vapor phase
    boxL_val = box.boxL[0] if hasattr(box.boxL, '__len__') else box.boxL
    print(f"Box dimensions: {boxL_val:.3f} x {boxL_val:.3f} x {boxL_val:.3f}")
    
    # Calculate reduced density
    sigma = 1.0  # LJ sigma in reduced units
    reduced_density = 10 * (sigma**3) / box.volume
    
    # Place atoms randomly in the box - use first 100 slots in the allocated arrays
    box.atoms[:natoms] = np.random.uniform(-boxL_val/2, boxL_val/2, size=(natoms, 3))
    
    # Set up indexing arrays for first 100 atoms
    box.AtomType[:natoms] = 0                      # All atoms are type 0 (LJ)
    box.MolType[:natoms] = 0                       # All molecules are type 0
    box.MolIndx[:natoms] = np.arange(natoms)          # Each atom is its own molecule
    box.MolSubIndx[:natoms] = np.arange(natoms)       # Sub-index within molecule type
    box.AtomSubIndx[:natoms] = np.zeros(natoms, dtype=int)  # Sub-index within atom type (all 0 for single-atom molecules)
    
    print(f"Initial atom positions generated for {box.nAtoms} atoms")
    print(f"Position range: [{np.min(box.atoms[:natoms]):.3f}, {np.max(box.atoms[:natoms]):.3f}]")
    
    return box

def create_fcc_with_vacancies():
    """
    Create a small FCC (Face-Centered Cubic) lattice with two atoms missing.
    Useful for testing neural network on solid-state structures with defects.
    
    Returns:
        box: CubeBox with FCC lattice structure and two vacancies
    """
    print("Creating FCC lattice with two vacancies...")
    
    # FCC lattice parameters
    lattice_constant = sqrt(8)*2**(1/6)/2.0  # Distance between nearest neighbors (in LJ sigma units)
    n_cells = 2  # Number of unit cells in each direction (2x2x2 = 32 atoms - 2 vacancies = 30 atoms)
    
    # Define molecule topology (single LJ atom per molecule)
    mock_lj_data = {
        "atoms": [("Ar", "LJ")],  # Single LJ atom type
    }
    atomtypes = ["LJ"]
    LJ_type = Molecule_Type(mock_lj_data, atomtypes=atomtypes)
    
    # Calculate number of atoms in FCC lattice
    # FCC has 4 atoms per unit cell
    natoms_full = 4 * (n_cells ** 3)
    natoms_with_vacancies = natoms_full - 2  # Remove 2 atoms
    
    # Set up box
    NMol = [natoms_with_vacancies]
    NMolMin = [0]
    NMolMax = [natoms_with_vacancies]  # Match max atoms to actual atom count to avoid mismatch
    
    # Create cubic box
    box = CubeBox([LJ_type], NMolMin=NMolMin, NMolMax=NMolMax, NMol=NMol)
    
    # Set box dimensions to fit the FCC lattice
    box_length = n_cells * lattice_constant
    
    box.load_dimension([box_length])
    box.boxID = 1
    box.nDimension = 3
    
    # Set molecule count arrays
    box.NMol = np.array(NMol)
    box.NMolMin = np.array(NMolMin)
    box.NMolMax = np.array(NMolMax)
    
    # Initialize atom counts and arrays
    box.nAtoms = natoms_with_vacancies
    box.nMaxAtoms = natoms_with_vacancies
    box.maxMol = len(NMol)
    box.nMolTotal = natoms_with_vacancies
    
    # Allocate arrays
    box.atoms = np.zeros((box.nMaxAtoms, 3), dtype=np.float64)
    box.AtomType = np.zeros(box.nMaxAtoms, dtype=int)
    box.MolType = np.zeros(box.nMaxAtoms, dtype=int)
    box.MolIndx = np.zeros(box.nMaxAtoms, dtype=int)
    box.MolSubIndx = np.zeros(box.nMaxAtoms, dtype=int)
    box.AtomSubIndx = np.zeros(box.nMaxAtoms, dtype=int)
    
    # Generate FCC lattice positions
    # FCC unit cell basis (in units of lattice constant):
    # (0, 0, 0), (0.5, 0.5, 0), (0.5, 0, 0.5), (0, 0.5, 0.5)
    fcc_basis = np.array([
        [0.0, 0.0, 0.0],
        [0.5, 0.5, 0.0],
        [0.5, 0.0, 0.5],
        [0.0, 0.5, 0.5]
    ]) * lattice_constant
    
    # Generate all FCC positions
    positions = []
    for i in range(n_cells):
        for j in range(n_cells):
            for k in range(n_cells):
                cell_origin = np.array([i, j, k]) * lattice_constant
                for basis_pos in fcc_basis:
                    positions.append(cell_origin + basis_pos)
    
    positions = np.array(positions)
    
    # Center the lattice in the box (shift from corner to center)
    positions -= box_length / 2.0
    
    # Remove two atoms to create vacancies
    # Choose atoms near the center of the lattice for interesting dynamics
    center = np.array([0.0, 0.0, 0.0])
    distances = np.linalg.norm(positions - center, axis=1)
    
    # Remove the two atoms closest to the center
    vacancy_indices = np.argsort(distances)[:2]
    mask = np.ones(len(positions), dtype=bool)
    mask[vacancy_indices] = False
    positions = positions[mask]
    
    print(f"Created FCC lattice: {n_cells}x{n_cells}x{n_cells} unit cells")
    print(f"Total atoms in full lattice: {natoms_full}")
    print(f"Atoms after removing 2 vacancies: {natoms_with_vacancies}")
    print(f"Vacancy positions (removed):")
    for idx in vacancy_indices:
        print(f"  Vacancy {idx}: {positions[idx] if idx < len(positions) else 'removed'}")
    
    # Place atoms in box
    box.atoms[:natoms_with_vacancies] = positions
    
    # Set up indexing arrays
    box.AtomType[:natoms_with_vacancies] = 0  # All atoms are type 0 (LJ)
    box.MolType[:natoms_with_vacancies] = 0  # All molecules are type 0
    box.MolIndx[:natoms_with_vacancies] = np.arange(natoms_with_vacancies)  # Each atom is its own molecule
    box.MolSubIndx[:natoms_with_vacancies] = np.arange(natoms_with_vacancies)  # Sub-index within molecule type
    box.AtomSubIndx[:natoms_with_vacancies] = np.zeros(natoms_with_vacancies, dtype=int)  # Sub-index within atom type
    
    boxL_val = box.boxL[0] if hasattr(box.boxL, '__len__') else box.boxL
    print(f"Box dimensions: {boxL_val:.3f} x {boxL_val:.3f} x {boxL_val:.3f}")
    print(f"Lattice constant: {lattice_constant:.3f} sigma")
    print(f"Position range: [{np.min(box.atoms[:natoms_with_vacancies]):.3f}, {np.max(box.atoms[:natoms_with_vacancies]):.3f}]")
    
    # Calculate average nearest neighbor distance for verification
    if natoms_with_vacancies > 1:
        sample_pos = box.atoms[0]
        distances = []
        for i in range(1, min(natoms_with_vacancies, 10)):
            dist = np.linalg.norm(box.atoms[i] - sample_pos)
            distances.append(dist)
        if distances:
            print(f"Sample nearest neighbor distances: min={min(distances):.3f}, mean={np.mean(distances):.3f}")
    
    return box

def setup_lj_forcefield(box):
    """
    Set up LJ forcefield with specified parameters.
    """
    print("\nSetting up LJ forcefield...")
    
    # Create LJ forcefield for 1 atom type
    lj_ff = LJ_Cut(nAtomTypes=1)
    
    # Set LJ parameters (reduced units)
    # epsilon = 1.0 (energy unit)
    # sigma = 1.0 (length unit)
    lj_ff.epsilon[0] = 1.0
    lj_ff.sigma[0] = 1.0
    lj_ff.rMin[0] = 0.4  # rMin is typically sigma * 2^(-1/6) ≈ 0.89
    
    # Set cutoff distance
    cutoff = 2.5  # 2.5 sigma units
    lj_ff.rCut = cutoff
    lj_ff.rCutSq = cutoff**2
    
    # Update mixing tables (for single type, just copy parameters)
    lj_ff.epsTable[0, 0] = 4.0 * lj_ff.epsilon[0]  # 4*epsilon for LJ formula
    lj_ff.sigTable[0, 0] = lj_ff.sigma[0]**2        # sigma^2 for optimization
    lj_ff.rMinTable[0, 0] = lj_ff.rMin[0]**2        # rMin^2 for optimization
    
    # Add forcefield to box
    box.EFunc.append(lj_ff)
    
    print(f"LJ parameters:")
    print(f"  epsilon = {lj_ff.epsilon[0]:.3f}")
    print(f"  sigma = {lj_ff.sigma[0]:.3f}")
    print(f"  rCut = {lj_ff.rCut:.3f}")
    print(f"  rMin = {lj_ff.rMin[0]:.3f}")
    
    return lj_ff

def setup_temperature(box, temperature):
    """
    Set up temperature and related thermodynamic variables.
    """
    print(f"\nSetting temperature to {temperature:.3f} (reduced units)...")
    
    # In reduced units: kT/epsilon = temperature
    # Since epsilon = 1.0, kT = temperature
    box.temperature = temperature
    box.beta = 1.0 / temperature  # 1/(kT) in reduced units
    
    print(f"Temperature (T*): {box.temperature:.3f}")
    print(f"Beta (1/kT*): {box.beta:.6f}")

def setup_nntranslate_move(box):
    """
    Set up the neural network translation move.
    """
    print("\nSetting up NNTranslate move...")
    
    # Create NNTranslate move
    nn_move = NNTranslate([box])
    
    # Create and set a simple neural network
    # (In practice, this would be a trained model)
    nn_move.create_and_set_neural_network("model.pt")
    
    print(f"NNTranslate parameters:")
    print(f"  Neighbor cutoff: {nn_move.neighbor_cutoff:.3f}")
    print(f"  NN cutoff: {nn_move.cutoff_rc:.3f}")
    print(f"  Translation box: {nn_move.translation_box}")
    print(f"  Number of bins: {nn_move.nbins}")
    
    return nn_move

def compute_initial_energy(box):
    """
    Compute initial energy of the system.
    """
    print("\nComputing initial energy...")
    success = box.compute_energy()
    if not success:
        raise RuntimeError("Initial energy computation failed")
    
    print(f"Initial energy (E*): {box.ETotal:.6f}")
    print(f"Initial energy per particle: {box.ETotal/box.nAtoms:.6f}")
    
    return box.ETotal

def run_simulation(box, nn_move, n_cycles=50000, moves_per_cycle=30):
    """
    Run the Monte Carlo simulation with statistics tracking.
    """
    print(f"\nStarting simulation: {n_cycles} cycles, {moves_per_cycle} moves per cycle")
    print(f"Total moves: {n_cycles * moves_per_cycle}")
    
    # Set up Metropolis sampling
    sampling = Metropolis()
    
    # Set up Monte Carlo simulation
    mc = SimMonteCarlo(
        nCycles=n_cycles, 
        nMoves=moves_per_cycle, 
        screenfreq=100,      # Print every cycle
        configfreq=0,      # Don't save configurations
        energyCheck=5000      
    )
    
    mc.BoxList = [box]
    
    uniform_move = MolTranslate([box], limit=8.0, max_dist=0.1)
    uniform_move.maintFreq = 10
    mc.Moves = [uniform_move]
    mc.moveweights = [25.0]
    
    #Add NNTranslate move
    mc.Moves.append(nn_move)
    mc.moveweights.append(1.0)
    
    mc.Sampling = sampling
    
    # Track statistics
    initial_energy = box.ETotal
    
    print(f"\nStarting energy: {initial_energy:.6f}")
    print("=" * 60)
    
    start_time = time.time()
    
    # Run simulation
    mc.run_monte_carlo()
    
    end_time = time.time()
    simulation_time = end_time - start_time
    
    # Final statistics
    final_energy = box.ETotal
    total_attempts = nn_move.atmps
    total_accepts = nn_move.accpt
    acceptance_rate = total_accepts / total_attempts if total_attempts > 0 else 0.0
    
    print("=" * 60)
    print("SIMULATION SUMMARY")
    print("=" * 60)
    print(f"Total simulation time: {simulation_time:.2f} seconds")
    print(f"Time per move: {simulation_time/(n_cycles*moves_per_cycle)*1000:.3f} ms")
    print(f"Moves per second: {(n_cycles*moves_per_cycle)/simulation_time:.1f}")
    print()
    print(f"Total move attempts: {int(total_attempts)}")
    print(f"Total move accepts: {int(total_accepts)}")
    print(f"Overall acceptance rate: {acceptance_rate:.4f} ({acceptance_rate*100:.2f}%)")
    print()
    print(f"Initial energy: {initial_energy:.6f}")
    print(f"Final energy: {final_energy:.6f}")
    print(f"Energy change: {final_energy - initial_energy:.6f}")
    
    if box.nAtoms > 0:
        print(f"Final energy per particle: {final_energy/box.nAtoms:.6f}")
    else:
        print(f"Final energy per particle: N/A (no atoms)")
    
    # Detailed rejection statistics
    print()
    print("DETAILED STATISTICS")
    print("-" * 30)
    print(f"Overlap rejections: {nn_move.ovlaprej}")
    print(f"Constraint rejections: {nn_move.constrainrej}")
    print(f"Detailed balance rejections: {nn_move.detailedrej}")
    print(f"Neural network rejections: {nn_move.nnrej}")
    
    return {
        'acceptance_rate': acceptance_rate,
        'initial_energy': initial_energy,
        'final_energy': final_energy,
        'energy_per_particle': final_energy/box.nAtoms if box.nAtoms > 0 else float('nan'),
        'total_attempts': int(total_attempts),
        'total_accepts': int(total_accepts),
        'simulation_time': simulation_time,
        'rejection_stats': {
            'overlap': nn_move.ovlaprej,
            'constraint': nn_move.constrainrej,
            'detailed_balance': nn_move.detailedrej,
            'neural_network': nn_move.nnrej
        }
    }
#========================================================
def main():
    """
    Main simulation script.
    """
    print("="*80)
    print("NEURAL NETWORK TRANSLATION MOVE SIMULATION")
    print("100 LJ Atoms in Vapor Phase")
    print("="*80)
    
    # Simulation parameters
    temperature = 0.7  # Reduced temperature (kT/epsilon)
    
    try:
        # Step 1: Create vapor box
        #box = create_vapor_box()
        box = create_fcc_with_vacancies()
        
        # Step 2: Set up LJ forcefield
        lj_ff = setup_lj_forcefield(box)
        
        # Step 3: Set temperature
        setup_temperature(box, temperature)
        
        # Step 4: Set up NNTranslate move
        nn_move = setup_nntranslate_move(box)
        
        # Step 5: Compute initial energy
        initial_energy = compute_initial_energy(box)
        
        # Step 6: Run simulation
        results = run_simulation(box, nn_move)
        
        
        # Save final configuration
        boxL_val = box.boxL[0] if hasattr(box.boxL, '__len__') else box.boxL
        print(f"\nFinal configuration:")
        print(f"  Box dimensions: {boxL_val:.3f}^3")
        print(f"  Number of atoms: {box.nAtoms}")
        print(f"  Final energy: {results['final_energy']:.6f}")
        
        print("\nSimulation completed successfully!")
        return results
        
    except Exception as e:
        print(f"\nERROR: Simulation failed with exception: {e}")
        import traceback
        traceback.print_exc()
        return None

if __name__ == "__main__":
    results = main()
    if results is not None:
        print(f"\nKey Results:")
        print(f"  Acceptance Rate: {results['acceptance_rate']:.4f}")
        if not np.isnan(results['energy_per_particle']):
            print(f"  Final Energy/Particle: {results['energy_per_particle']:.6f}")
        else:
            print(f"  Final Energy/Particle: N/A (no atoms)")
        print(f"  Simulation Time: {results['simulation_time']:.2f} seconds")
