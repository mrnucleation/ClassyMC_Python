from time import time
from src_python.Script_LoadCoordinates import load_coords
from src_python.Molecule_Definition import Molecule_Type
from src_python.FF_LJ_Cut import LJ_Cut
from src_python.CoordinateTypes import Displacement
from src_python.MC_Move_MolTranslation import MolTranslate
from src_python.Sampling_Metropolis import Metropolis
from src_python.Sim_MonteCarlo import SimMonteCarlo
from src_python.Box_SimpleBox import SimpleBox

import numpy as np

#Seed the random number generator for reproducibility
np.random.seed(42)

def test_sim_monte_carlo_comprehensive():
    """Comprehensive test of SimMonteCarlo class following TestLJStack.py flow"""
    
    print("=== Setting up simulation components ===")
    
    # Mock LJ data (same as TestLJStack.py)
    mock_lj_data = {
        "atoms": [("Ar", "LJ")],
    }
    atomtypes = ["LJ"]
    LJ_type = Molecule_Type(mock_lj_data, atomtypes=atomtypes)

    # Load coordinates from file
    box = load_coords("SimpleLJ.clssy", [LJ_type])
    
    print(f"Box ID: {box.boxID}")
    print(f"Box Dimensions: {box.boxL}")
    print(f"Number of Molecules: {box.NMol}")
    print(f"Molecule Index: {box.MolIndx}")
    print(f"Molecule Sub Index: {box.MolSubIndx}")
    print(f"Atom Type: {box.AtomType}")
    print(f"Atom Sub Index: {box.AtomSubIndx}")
    print(f"Atom Positions: {box.atoms}")
    
    # Define LJ Forcefield
    lj_ff = LJ_Cut(nAtomTypes=1)
    lj_ff.rMin = np.zeros(1)
    lj_ff.rMinTable = np.zeros((1,1))
    print(f"Number of Atom Types: {lj_ff.nAtomTypes}")
    print(f"Epsilon values: {lj_ff.epsilon}")
    print(f"Sigma values: {lj_ff.sigma}")
    
    box.EFunc.append(lj_ff)
    
    print(box.EFunc)
    print("LJ Forcefield initialized and added to box.")
    
    # Compute the initial energy
    success = box.compute_energy()
    assert success, "Energy computation failed"
    print("Initial energy computed successfully.")
    print(f"Initial Energy: {box.ETotal}")
    
    start_energy = box.ETotal
    
    # Test manual displacement (like in TestLJStack.py)
    print("\\n=== Testing manual displacement ===")
    disp = Displacement(LJ_type, 0, 0, box.atoms[0])
    print(f"Displacement created: {disp}")
    delta_x = np.array([1.12, 0.0, 0.0])
    disp.X = box.atoms[0] + delta_x
    
    # Compute the energy after displacement
    E_Inter, E_Intra, accept = box.compute_energy_delta(disp)
    assert accept, "Energy computation after displacement failed"
    print(f"Delta x energy: {E_Inter}")
    
    box.update_position(disp)
    new_running_energy = start_energy + E_Inter
    
    # Recompute the energy with the new positions
    success = box.compute_energy()
    assert success, "Energy computation failed"
    
    print(f"Energy difference: {box.ETotal - new_running_energy}")
    
    # Test index functions
    mol_indicies = box.get_molindicies()
    print(f"Molecule Indices: {mol_indicies}")
    
    print("\\n=== Setting up SimMonteCarlo ===")
    
    # Create SimMonteCarlo instance
    sim_mc = SimMonteCarlo(
        nCycles=100,          # Number of cycles
        nMoves=10,           # Moves per cycle
        screenfreq=25,       # Screen output frequency
        configfreq=50,       # Configuration output frequency
        energyCheck=25       # Energy check frequency
    )
    
    # Add the box to the simulation
    sim_mc.BoxList = [box]
    
    # Create Metropolis sampling rule
    sampling = Metropolis()
    sim_mc.Sampling = sampling
    
    # Create molecular translation move
    transMove = MolTranslate([box])
    sim_mc.Moves = [transMove]
    
    # Add molecular data
    sim_mc.MolData = [LJ_type]
    
    print("SimMonteCarlo setup complete.")
    print(f"Number of boxes: {len(sim_mc.BoxList)}")
    print(f"Number of moves: {len(sim_mc.Moves)}")
    print(f"Sampling rule: {type(sim_mc.Sampling).__name__}")
    
    # Run the Monte Carlo simulation
    print("\\n=== Running Monte Carlo simulation ===")
    start_time = time()
    
    try:
        sim_mc.run_monte_carlo()
        end_time = time()
        print(f"Simulation completed successfully in {end_time - start_time:.3f} seconds")
        print(f"Final energy: {box.ETotal}")
        print(f"Total moves attempted: {transMove.atmps}")
        print(f"Total moves accepted: {transMove.accpt}")
        print(f"Acceptance rate: {transMove.get_accept_rate():.2f}%")
        
        # Verify energy consistency
        success = box.compute_energy()
        assert success, "Final energy computation failed"
        print(f"Final energy verified: {box.ETotal}")
        
    except Exception as e:
        print(f"Error during simulation: {e}")
        import traceback
        traceback.print_exc()

def test_sim_monte_carlo_minimal():
    """Minimal test of SimMonteCarlo class"""
    
    print("\\n=== Minimal SimMonteCarlo Test ===")
    
    # Mock LJ data
    mock_lj_data = {
        "atoms": [("Ar", "LJ")],
    }
    atomtypes = ["LJ"]
    LJ_type = Molecule_Type(mock_lj_data, atomtypes=atomtypes)

    # Load coordinates
    box = load_coords("SimpleLJ.clssy", [LJ_type])
    box.temp = 0.6
    box.beta = 1.0 / box.temp
    
    # Define LJ Forcefield
    lj_ff = LJ_Cut(nAtomTypes=1)
    lj_ff.rMin = np.zeros(1)
    lj_ff.rMinTable = np.zeros((1,1))
    box.EFunc.append(lj_ff)
    
    # Compute initial energy
    box.compute_energy()
    
    # Create minimal SimMonteCarlo instance
    sim_mc = SimMonteCarlo(
        nCycles=100,           # Small number of cycles
        nMoves=10,            # Small number of moves per cycle
        screenfreq=1,        # Frequent screen output
        configfreq=5,        # Configuration output frequency
        energyCheck=1        # Frequent energy checks
    )
    
    # Setup required components
    sim_mc.BoxList = [box]
    sim_mc.Sampling = Metropolis()
    sim_mc.Moves = [MolTranslate([box])]
    sim_mc.MolData = [LJ_type]
    
    print("Minimal setup complete. Running simulation...")
    
    sim_mc.run_monte_carlo()
    try:
        print("Minimal simulation completed successfully!")
    except Exception as e:
        print(f"Error in minimal simulation: {e}")

#-------------------------------

if __name__ == "__main__":
    print("Running comprehensive SimMonteCarlo tests...")
    
    # Run the comprehensive test
    test_sim_monte_carlo_comprehensive()
    
    # Run the minimal test
    test_sim_monte_carlo_minimal()
    
    print("\\nAll SimMonteCarlo tests completed.")
