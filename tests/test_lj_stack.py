
from time import time
import math
from src_python.Script_LoadCoordinates import load_coords
from src_python.Molecule_Definition import Molecule_Type
from src_python.FF_LJ_Cut import LJ_Cut
from src_python.CoordinateTypes import Displacement, Addition, Deletion
from src_python.MC_Move_MolTranslation import MolTranslate
from src_python.CoordinateTypes import OrthoVolChange
from src_python.Sampling_Metropolis import Metropolis
from src_python.Sim_MonteCarlo import SimMonteCarlo

import numpy as np

#Seed the random number generator for reproducibility
np.random.seed(42)

def test_lj_stack():
    
    
    mock_lj_data = {
        "atoms": [("Ar", "LJ")],
    }
    atomtypes = ["LJ"]
    LJ_type = Molecule_Type(mock_lj_data, atomtypes=atomtypes)



    
    box = load_coords("SimpleLJ.clssy", [LJ_type])
    
    
    print(f"Box ID: {box.boxID}")
    print(f"Box Dimensions: {box.boxL}")
    print(f"Number of Molecules: {box.NMol}")
    print(f"Molecule Index: {box.MolIndx}")
    print(f"Molecule Sub Index: {box.MolSubIndx}")
    print(f"Atom Type: {box.AtomType}")
    print(f"Atom Sub Index: {box.AtomSubIndx}")
    print(f"Atom Positions: {box.atoms}")
    
    #Define LJ Forcefield
    lj_ff = LJ_Cut(nAtomTypes=1)
    lj_ff.rMin = np.zeros(1)
    lj_ff.rMinTable = np.zeros((1,1))
    print(f"Number of Atom Types: {lj_ff.nAtomTypes}")
    print(f"Epsilon values: {lj_ff.epsilon}")
    print(f"Sigma values: {lj_ff.sigma}")
    
    box.EFunc.append(lj_ff)
    
    print(box.EFunc)
    print("LJ Forcefield initialized and added to box.")
    
    #Compute the initial energy
    success = box.compute_energy()
    assert success, "Energy computation failed"
    print("Initial energy computed successfully.")
    print(f"Initial Energy: {box.E_Inter}")
    
    start_energy = box.ETotal
    
    #Set up a displacement object
    # Move atom 0 (first atom of first molecule)
    mol_type = 0  # First molecule type
    mol_index = 0  # First molecule
    atom_indices = np.array([0])  # First atom
    delta_x = np.array([1.12, 0.0, 0.0])
    new_positions = (box.atoms[0] + delta_x).reshape(1, -1)  # Reshape to (1, 3)

    disp = Displacement(mol_type, mol_index, atom_indices, new_positions)
    print(f"Displacement created: {disp}")
    print(f"Original position: {box.atoms[0]}")
    print(f"New position: {new_positions}")
    
    # Compute the energy after displacement
    E_Inter, E_Intra, accept = box.compute_energy_delta(disp)
    assert accept, "Energy computation after displacement failed"
    #Update the box with the new positions and energy
    
    print(f"Delta x: {E_Inter}")
    
    box.update_position(disp)
    #box.update_energy(E_Inter)
    
    new_running_energy = start_energy + E_Inter
    
    #Recompute the energy with the new positions
    success = box.compute_energy()
    assert success, "Energy computation failed"
    
    print(f"Energy difference: {box.ETotal - new_running_energy}")
    
    
    #Test index functions
    mol_indicies = box.get_molindicies()
    print(f"Molecule Indices: {mol_indicies}")
    
    sampling = Metropolis()
    
    # Create a displacement move for the LJ stack
    transMove = MolTranslate(
        [box], 
    )

    
    timestart = time()
    for imove in range(1000):
        accept = transMove.full_move(box, sampling)
    timeend = time()
    E_culmative = box.ETotal
    print(f"Time taken for 1000 moves: {timeend - timestart} seconds")
    #Recompute the energy after all moves
    success = box.compute_energy()
    assert success, "Energy computation failed after moves"
    assert np.isclose(box.ETotal, E_culmative), "Cumulative energy does not match computed energy"

def test_lj_stack_montecarlo():
    # Reuse the same LJ stack setup
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

    # Initial energy
    ok = box.compute_energy()
    assert ok, "Initial energy computation failed"

    # Moves and sampling
    sampling = Metropolis()
    trans_move = MolTranslate([box])

    # Configure and run Monte Carlo
    mc = SimMonteCarlo(nCycles=2, nMoves=50, screenfreq=1, configfreq=0, energyCheck=0)
    mc.BoxList = [box]
    mc.Moves = [trans_move]
    mc.Sampling = sampling

    mc.run_monte_carlo()

def test_lj_volume_changes():
    """
    Test volume change moves and related functionality for LJ system.
    Tests the various box and energy routines for volume changes.
    """
    print("\n=== Testing LJ Volume Changes ===\n")

    # Set up LJ system similar to existing tests
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

    # Initial energy computation
    success = box.compute_energy()
    assert success, "Initial energy computation failed"
    initial_energy = box.ETotal
    initial_volume = box.volume
    initial_box_dims = box.boxL.copy()

    print(f"Initial system setup:")
    print(f"  Box dimensions: {initial_box_dims}")
    print(f"  Volume: {initial_volume[0]:.6f}")
    print(f"  Energy: {initial_energy:.6f}")
    print(f"  Number of molecules: {box.NMol}")

    # Test 1: Box volume calculation and scaling
    print("\n--- Test 1: Box Volume Scaling ---")
    test_scale_factor = 1.1  # 10% increase
    new_volume = initial_volume * (test_scale_factor ** 3)
    new_box_dims = initial_box_dims * test_scale_factor

    print(f"  Original volume: {initial_volume[0]:.6f}")
    print(f"  Scale factor: {test_scale_factor}")
    print(f"  Expected new volume: {new_volume[0]:.6f}")
    print(f"  Expected new box dimensions: {new_box_dims}")

    # Test 2: OrthoVolChange displacement object
    print("\n--- Test 2: OrthoVolChange Object ---")
    scales = np.array([test_scale_factor, test_scale_factor, test_scale_factor])
    ortho_vol_disp = OrthoVolChange(scales=scales, volOld=initial_volume[0], volNew=new_volume[0])

    print(f"  OrthoVolChange created successfully")
    print(f"  Scales: {ortho_vol_disp.scales}")
    print(f"  Old volume: {ortho_vol_disp.volOld:.6f}")
    print(f"  New volume: {ortho_vol_disp.volNew:.6f}")

    # Test 3: Energy computation with volume change
    print("\n--- Test 3: Energy with Volume Change ---")

    # Make a temporary copy of atom positions for scaling
    original_positions = box.atoms.copy()

    # Scale the positions
    scaled_positions = original_positions * test_scale_factor
    box.atoms = scaled_positions
    box.boxL = new_box_dims

    # Update box volume
    box.volume = new_volume

    # Compute energy after volume change
    success = box.compute_energy()
    assert success, "Energy computation after volume change failed"
    scaled_energy = box.ETotal

    print(f"  Energy after scaling: {scaled_energy:.6f}")
    print(f"  Energy difference: {scaled_energy - initial_energy:.6f}")

    # Restore original state
    box.atoms = original_positions
    box.boxL = initial_box_dims
    box.volume = initial_volume
    box.compute_energy()  # Restore original energy

    # Test 4: Manual volume move simulation
    print("\n--- Test 4: Manual Volume Move Simulation ---")

    # Simulate what a volume move does without using the IsoVol class
    # This tests the core volume change logic

    # Propose volume change (log scale)
    maxDv = 0.01
    dV = maxDv * (2.0 * np.random.random() - 1.0)  # Random value between -maxDv and maxDv
    vol_new = initial_volume * math.exp(dV)

    # Calculate scale factors
    scale_factor = (vol_new / initial_volume) ** (1.0 / 3.0)
    new_box_dims_proposed = initial_box_dims * scale_factor
    new_positions_proposed = original_positions * scale_factor

    print(f"  Proposed volume change:")
    print(f"    dV: {dV:.6f}")
    print(f"    Old volume: {initial_volume[0]:.6f}")
    print(f"    New volume: {vol_new[0]:.6f}")
    print(f"    Scale factor: {scale_factor[0]:.6f}")

    # Apply the proposed change
    box.atoms = new_positions_proposed
    box.boxL = new_box_dims_proposed
    box.volume = vol_new

    # Compute energy with the proposed change
    success = box.compute_energy()
    assert success, "Energy computation with proposed volume change failed"
    proposed_energy = box.ETotal

    print(f"    Proposed energy: {proposed_energy:.6f}")
    print(f"    Energy change: {proposed_energy - initial_energy:.6f}")

    # Calculate acceptance probability (simplified, without Metropolis)
    n_mol = box.nMolTotal
    prob = (n_mol + 1) * math.log(vol_new[0] / initial_volume[0])
    print(f"    Log probability term: {prob:.6f}")

    # Test 5: Multiple volume change attempts
    print("\n--- Test 5: Multiple Volume Change Attempts ---")

    # Reset to original state
    box.atoms = original_positions
    box.boxL = initial_box_dims
    box.volume = initial_volume
    box.compute_energy()

    n_attempts = 5
    accepted_attempts = 0

    for i in range(n_attempts):
        # Store original state
        orig_energy = box.ETotal
        orig_volume = box.volume
        orig_positions = box.atoms.copy()
        orig_boxL = box.boxL.copy()

        # Propose volume change
        dV = maxDv * (2.0 * np.random.random() - 1.0)
        vol_new = orig_volume * math.exp(dV)
        scale_factor = (vol_new / orig_volume) ** (1.0 / 3.0)

        # Apply change
        box.atoms = orig_positions * scale_factor
        box.boxL = orig_boxL * scale_factor
        box.volume = vol_new

        # Compute new energy
        success = box.compute_energy()
        if success:
            new_energy = box.ETotal
            # Simple acceptance criterion (just check if energy decreased)
            if new_energy < orig_energy:
                accepted_attempts += 1
                print(f"    Attempt {i+1}: ACCEPTED (ΔE = {new_energy - orig_energy:.6f})")
            else:
                # Reject - restore original state
                box.atoms = orig_positions
                box.boxL = orig_boxL
                box.volume = orig_volume
                box.ETotal = orig_energy
                print(f"    Attempt {i+1}: REJECTED (ΔE = {new_energy - orig_energy:.6f})")
        else:
            # Reject due to energy computation failure
            box.atoms = orig_positions
            box.boxL = orig_boxL
            box.volume = orig_volume
            print(f"    Attempt {i+1}: REJECTED (energy computation failed)")

    acceptance_rate = (accepted_attempts / n_attempts) * 100
    print(f"  Manual volume move acceptance rate: {acceptance_rate:.1f}% ({accepted_attempts}/{n_attempts})")

    # Test 6: Volume scaling edge cases
    print("\n--- Test 6: Volume Scaling Edge Cases ---")

    # Reset to original state
    box.atoms = original_positions
    box.boxL = initial_box_dims
    box.volume = initial_volume
    box.compute_energy()

    # Test very small volume change
    small_scale = 1.001  # 0.1% change
    box.atoms = original_positions * small_scale
    box.boxL = initial_box_dims * small_scale
    box.volume = initial_volume * (small_scale ** 3)
    success = box.compute_energy()
    print(f"  Small volume change (0.1%): {'SUCCESS' if success else 'FAILED'}")

    # Reset
    box.atoms = original_positions
    box.boxL = initial_box_dims
    box.volume = initial_volume
    box.compute_energy()

    # Test larger volume change
    large_scale = 1.5  # 50% increase
    box.atoms = original_positions * large_scale
    box.boxL = initial_box_dims * large_scale
    box.volume = initial_volume * (large_scale ** 3)
    success = box.compute_energy()
    print(f"  Large volume change (50%): {'SUCCESS' if success else 'FAILED'}")

    # Reset to original state
    box.atoms = original_positions
    box.boxL = initial_box_dims
    box.volume = initial_volume
    box.compute_energy()

    print("\n=== LJ Volume Change Tests Completed ===\n")

def test_lj_addition_deletion():
    """
    Test addition and deletion moves for LJ system.
    This test verifies that energy calculations are consistent between
    delta calculations and full energy recomputations for addition and deletion operations.
    """
    print("\n=== Testing LJ Addition and Deletion Moves ===\n")

    # Set up LJ system similar to other tests
    mock_lj_data = {
        "atoms": [("Ar", "LJ")],
    }
    atomtypes = ["LJ"]
    LJ_type = Molecule_Type(mock_lj_data, atomtypes=atomtypes)

    box = load_coords("SimpleLJ.clssy", [LJ_type])

    # Define LJ Forcefield
    lj_ff = LJ_Cut(nAtomTypes=1)
    lj_ff.rMin = np.zeros(1)
    lj_ff.rMinTable = np.zeros((1,1))
    box.EFunc.append(lj_ff)

    # Initial energy computation
    success = box.compute_energy()
    assert success, "Initial energy computation failed"
    initial_energy = box.ETotal
    initial_nmol = box.NMol.copy()
    initial_natoms = box.nAtoms

    print(f"Initial system:")
    print(f"  Energy: {initial_energy:.6f}")
    print(f"  Number of molecules: {initial_nmol}")
    print(f"  Number of atoms: {initial_natoms}")

    # Step 1: Test Addition Move
    print("\n--- Step 1: Addition Move ---")

    # Store original state
    original_atoms = box.atoms.copy()
    original_energy = box.ETotal

    # Create Addition object
    # Add a new molecule at a random position
    mol_type = 0  # LJ molecule type
    mol_sub_index = initial_nmol[mol_type]  # Next available molecule slot

    # Generate random position for new atom
    new_position = np.array([[np.random.uniform(0, box.boxL[0]),
                             np.random.uniform(0, box.boxL[1]),
                             np.random.uniform(0, box.boxL[2])]])
    
    print(f"New position: {new_position}")
    print(f"New position shape: {new_position.shape}")

    # Create Addition object
    addition = Addition(molType=mol_type,
                       molIndx=box.atoms.shape[0],
                       atomTypes=np.array([0]),  # New atom index
                       newPositions=new_position)

    print(f"  Created Addition object:")
    print(f"    Molecule type: {addition.molType}")
    print(f"    Molecule sub-index: {addition.molIndx}")
    print(f"    Atom indices: {addition.atomTypes}")
    print(f"    New positions: {addition.X}")

    # Step 2: Compute energy delta using diffcalc
    print("\n--- Step 2: Compute Energy Delta ---")
    E_Inter_delta, accept = box.EFunc[0].diff_calc(box, addition, [], [])
    assert accept, "Energy delta computation failed for addition"

    print(f"  Energy delta from diffcalc:")
    print(f"    Intermolecular: {E_Inter_delta:.6f}")

    # Step 3: Apply the addition using add_mol
    print("\n--- Step 3: Apply Addition ---")

    # Use add_mol to add the molecule
    box.update_position(addition)

    print(f"  Added molecule successfully")
    print(f"  New number of molecules: {box.NMol}")
    print(f"  New number of atoms: {box.nAtoms}")

    # Step 4: Compute new energy with delta method
    new_energy_delta = original_energy + E_Inter_delta 
    print(f"  New energy (delta method): {new_energy_delta:.6f}")

    # Step 5: Recompute energy from scratch
    print("\n--- Step 5: Full Energy Recomputation ---")
    success = box.compute_energy()
    assert success, "Full energy recomputation failed"
    new_energy_full = box.ETotal

    print(f"  New energy (full recomputation): {new_energy_full:.6f}")

    # Step 6: Compare energies
    energy_diff = abs(new_energy_delta - new_energy_full)
    print(f"  Energy difference: {energy_diff:.10f}")

    if energy_diff > 1e-10:
        raise AssertionError(f"Energy mismatch after addition! Delta: {new_energy_delta:.10f}, Full: {new_energy_full:.10f}, Diff: {energy_diff:.10f}")

    print("  ✓ Addition energy calculation is consistent!")

    # Step 7: Test Deletion Move
    print("\n--- Step 7: Deletion Move ---")

    # Store current state before deletion
    current_energy = box.ETotal
    current_atoms = box.atoms.copy()

    # Choose a molecule to delete (not the original ones, delete the one we just added)
    mol_to_delete = box.MolIndx.max()  # The molecule we just added
    print(f"  Deleting molecule sub-index: {mol_to_delete}")

    # Create Deletion object
    # Find atoms belonging to this molecule
    mol_data = box.get_mol_data(mol_to_delete)
    
    atom_indices = mol_data['atomIndicies']

    deletion = Deletion(molType=mol_type,
                       molIndx=mol_to_delete,
                       atomTypes=np.array([0]),
                       atomIndicies=atom_indices)

    print(f"  Created Deletion object:")
    print(f"    Molecule type: {deletion.molType}")
    print(f"    Molecule index: {deletion.molIndx}")
    print(f"    Atom indices: {deletion.atomIndicies}")

    # Step 8: Compute energy delta for deletion
    print("\n--- Step 8: Compute Deletion Energy Delta ---")
    E_Inter_delta_del, accept = box.EFunc[0].diff_calc(box, deletion, [], [])
    assert accept, "Energy delta computation failed for deletion"

    print(f"  Energy delta from diffcalc:")
    print(f"    Intermolecular: {E_Inter_delta_del:.6f}")

    # Step 9: Apply the deletion
    print("\n--- Step 9: Apply Deletion ---")
    box.update_position(deletion)

    print(f"  Deleted molecule successfully")
    print(f"  New number of molecules: {box.NMol}")
    print(f"  New number of atoms: {box.nAtoms}")

    # Step 10: Compute new energy with delta method
    new_energy_delta_del = current_energy + E_Inter_delta_del
    print(f"  New energy (delta method): {new_energy_delta_del:.6f}")

    # Step 11: Recompute energy from scratch
    print("\n--- Step 11: Full Energy Recomputation After Deletion ---")
    success = box.compute_energy()
    assert success, "Full energy recomputation failed after deletion"
    new_energy_full_del = box.ETotal

    print(f"  New energy (full recomputation): {new_energy_full_del:.6f}")

    # Step 12: Compare energies
    energy_diff_del = abs(new_energy_delta_del - new_energy_full_del)
    print(f"  Energy difference: {energy_diff_del:.10f}")

    if energy_diff_del > 1e-10:
        raise AssertionError(f"Energy mismatch after deletion! Delta: {new_energy_delta_del:.10f}, Full: {new_energy_full_del:.10f}, Diff: {energy_diff_del:.10f}")

    print("  ✓ Deletion energy calculation is consistent!")

    # Step 13: Verify we're back to original state
    print("\n--- Step 13: Verify Original State ---")
    final_energy_diff = abs(new_energy_full_del - initial_energy)
    print(f"  Final energy: {new_energy_full_del:.6f}")
    print(f"  Original energy: {initial_energy:.6f}")
    print(f"  Energy difference from original: {final_energy_diff:.10f}")

    if final_energy_diff > 1e-10:
        raise AssertionError(f"Not back to original energy state! Original: {initial_energy:.10f}, Final: {new_energy_full_del:.10f}, Diff: {final_energy_diff:.10f}")

    print("  ✓ Successfully returned to original state!")

    print("\n=== LJ Addition/Deletion Tests Completed Successfully ===\n")

#-------------------------------

if __name__ == "__main__":
    print("Running LJ Stack tests...")

    # Run the LJ stack test
    #test_lj_stack()
    # Run the LJ stack Monte Carlo driver test
    #test_lj_stack_montecarlo()
    # Run the volume change test
    #test_lj_volume_changes()
    # Run the addition/deletion test
    test_lj_addition_deletion()

    print("LJ Stack tests completed.")
