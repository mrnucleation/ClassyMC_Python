#!/usr/bin/env python3
"""
Production Monte Carlo simulation for FCC crystal structure with ~300 atoms.
This script creates an FCC lattice, sets up LJ forcefield, and runs
a long simulation to test translation move performance.
"""

import numpy as np
import math
from time import time

# Import the necessary modules from the simulation framework
from src_python.Script_LoadCoordinates import load_coords
from src_python.Molecule_Definition import Molecule_Type
from src_python.FF_LJ_Cut import LJ_Cut
from src_python.MC_Move_MolTranslation import MolTranslate
from src_python.Sampling_Metropolis import Metropolis
from src_python.Sim_MonteCarlo import SimMonteCarlo


def create_fcc_coordinates(n_unit_cells_per_side=4):
    """
    Create FCC crystal structure coordinates.

    Args:
        n_unit_cells_per_side: Number of unit cells along each dimension
                              Total atoms = 4 * (n_unit_cells_per_side)^3

    Returns:
        tuple: (coordinates, box_length)
    """
    # FCC unit cell lattice constant (typical for LJ systems)
    a = 5.0  # Angstroms

    # FCC basis atoms (relative to unit cell)
    basis = np.array([
        [0.0, 0.0, 0.0],
        [0.5, 0.5, 0.0],
        [0.5, 0.0, 0.5],
        [0.0, 0.5, 0.5]
    ])

    coordinates = []

    # Generate all unit cell positions
    for i in range(n_unit_cells_per_side):
        for j in range(n_unit_cells_per_side):
            for k in range(n_unit_cells_per_side):
                # Unit cell origin
                cell_origin = np.array([i, j, k]) * a

                # Add all 4 atoms in the basis
                for atom in basis:
                    pos = cell_origin + atom * a
                    coordinates.append(pos)

    coordinates = np.array(coordinates)
    box_length = n_unit_cells_per_side * a

    n_atoms = len(coordinates)
    print(f"Created FCC structure with {n_atoms} atoms")
    print(f"Box dimensions: {box_length:.2f} x {box_length:.2f} x {box_length:.2f} Å³")
    print(f"Number density: {n_atoms / (box_length**3):.4f} atoms/Å³")

    return coordinates, box_length


def create_clssy_file(coordinates, box_length, filename="FCC_LJ.clssy"):
    """
    Create a .clssy coordinate file for the simulation.

    Args:
        coordinates: Array of atomic coordinates
        box_length: Box dimension
        filename: Output filename
    """
    n_atoms = len(coordinates)

    with open(filename, 'w') as f:
        f.write(f"cube {box_length}\n")
        f.write(f"NMol {n_atoms}\n")  # Each atom is treated as a separate molecule
        f.write(f"NMax {n_atoms}\n")
        f.write("NMin 0\n")

        # Write each atom as a separate molecule
        for i, coord in enumerate(coordinates):
            x, y, z = coord
            f.write(f"0 {i} 0 {x:.6f} {y:.6f} {z:.6f}\n")

    print(f"Created coordinate file: {filename}")
    return filename


def setup_lj_system(coord_filename):
    """
    Set up the LJ system for simulation.

    Args:
        coord_filename: Path to coordinate file

    Returns:
        tuple: (box, molecule_type)
    """
    print("\n=== Setting up LJ System ===")

    # Define molecule type (single atom LJ)
    mock_lj_data = {
        "atoms": [("Ar", "LJ")],
    }
    atomtypes = ["LJ"]
    LJ_type = Molecule_Type(mock_lj_data, atomtypes=atomtypes)

    # Load coordinates
    box = load_coords(coord_filename, [LJ_type])
    box.temperature = 0.7
    box.beta = 1.0 / 0.7
    box.volume = box.boxL**3

    print(f"Box ID: {box.boxID}")
    print(f"Box Dimensions: {box.boxL}")
    print(f"Number of Molecules: {box.NMol}")
    print(f"Total Atoms: {box.nAtoms}")

    # Set up LJ forcefield
    lj_ff = LJ_Cut(nAtomTypes=1)
    lj_ff.epsilon = np.array([1.0])  # LJ epsilon parameter
    lj_ff.sigma = np.array([1.0])    # LJ sigma parameter
    lj_ff.rCut = 2.5                 # Cutoff distance
    lj_ff.rMin = np.zeros(1)
    lj_ff.rMinTable = np.zeros((1,1))

    box.EFunc.append(lj_ff)
    print("LJ Forcefield initialized and added to box.")

    # Compute initial energy
    success = box.compute_energy()
    if not success:
        raise RuntimeError("Initial energy computation failed")

    print(".6f")
    print(f"Initial Energy Breakdown: Inter={box.E_Inter:.6f}, Intra={box.E_Intra:.6f}")

    return box, LJ_type


def run_production_simulation(box, n_cycles=100, moves_per_cycle=300):
    """
    Run the production Monte Carlo simulation.

    Args:
        box: Simulation box
        n_cycles: Number of MC cycles
        moves_per_cycle: Number of moves per cycle
    """
    print(f"\n=== Running Production Simulation ===")
    print(f"Simulation Parameters:")
    print(f"  - Cycles: {n_cycles}")
    print(f"  - Moves per cycle: {moves_per_cycle}")
    print(f"  - Total moves: {n_cycles * moves_per_cycle}")

    # Set up Monte Carlo simulation
    sampling = Metropolis()

    # Create translation move
    trans_move = MolTranslate([box])
    trans_move.max_dist = 0.5  # Maximum translation distance in Angstroms
    trans_move.boxmax_dist[0] = 0.5  # Set for the first (and only) box

    # Configure simulation
    mc = SimMonteCarlo(
        nCycles=n_cycles,
        nMoves=moves_per_cycle,
        screenfreq=10,      # Print status every 1000 cycles
        configfreq=0,         # Don't save configurations
        energyCheck=1000      # Check energy conservation every 1000 cycles
    )

    mc.BoxList = [box]
    mc.Moves = [trans_move]
    mc.Sampling = sampling

    # Run simulation
    print("Starting simulation...")
    start_time = time()

    mc.run_monte_carlo()

    end_time = time()
    total_time = end_time - start_time

    # Print final statistics
    print("\n=== Simulation Complete ===")
    print(".3f")
    print(".1f")
    print(".6f")

    if hasattr(mc, 'acceptance_ratio'):
        print(".3f")

    return total_time


def main():
    """Main function to run the production simulation."""
    print("FCC Crystal Monte Carlo Production Run")
    print("=" * 50)

    # Create FCC structure with approximately 300 atoms
    # 4^3 = 64 unit cells * 4 atoms/unit = 256 atoms (close to 300)
    print("\n1. Generating FCC crystal structure...")
    coordinates, box_length = create_fcc_coordinates(n_unit_cells_per_side=4)

    # Create coordinate file
    coord_file = create_clssy_file(coordinates, box_length)

    # Set up the system
    print("\n2. Setting up simulation system...")
    box, lj_type = setup_lj_system(coord_file)

    # Run production simulation
    print("\n3. Running production simulation...")
    total_time = run_production_simulation(
        box,
        n_cycles=100,
        moves_per_cycle=300
    )

    print("\nProduction run completed successfully!")


if __name__ == "__main__":
    main()
