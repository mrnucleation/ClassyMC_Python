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

def test_sim_monte_carlo():
    """Test the SimMonteCarlo class with similar functionality to TestLJStack.py"""
    
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
    
    # Create SimMonteCarlo instance
    sim_mc = SimMonteCarlo(
        nCycles=10,           # Number of cycles
        nMoves=10,            # Moves per cycle
        screenfreq=5,         # Screen output frequency
        configfreq=10,        # Configuration output frequency
        energyCheck=5         # Energy check frequency
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
    print("Starting Monte Carlo simulation...")
    start_time = time()
    
    try:
        sim_mc.run_monte_carlo()
        end_time = time()
        print(f"Simulation completed successfully in {end_time - start_time:.3f} seconds")
        print(f"Final energy: {box.ETotal}")
        print(f"Total moves attempted: {transMove.atmps}")
        print(f"Total moves accepted: {transMove.accpt}")
        print(f"Acceptance rate: {transMove.get_accept_rate():.2f}%")
    except Exception as e:
        print(f"Error during simulation: {e}")
        import traceback
        traceback.print_exc()

#-------------------------------

if __name__ == "__main__":
    print("Running SimMonteCarlo tests...")
    
    # Run the SimMonteCarlo test
    test_sim_monte_carlo()
    
    print("SimMonteCarlo tests completed.")
