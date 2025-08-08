
from time import time
from src_python.Script_LoadCoordinates import load_coords
from src_python.Molecule_Definition import Molecule_Type
from src_python.FF_LJ_Cut import LJ_Cut
from src_python.CoordinateTypes import Displacement
from src_python.MC_Move_MolTranslation import MolTranslate
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
    disp = Displacement(LJ_type, 0, 0, box.atoms[0])
    print(f"Displacement created: {disp}")
    delta_x = np.array([1.12, 0.0, 0.0])
    disp.X = box.atoms[0] + delta_x
    
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

def test_sim_monte_carlo():
    """
    Test the Sim_MonteCarlo class with the same LJ stack setup
    """
    print("\n" + "="*50)
    print("Testing Sim_MonteCarlo with LJ Stack")
    print("="*50)
    
    # Set up the same LJ system as in test_lj_stack()
    mock_lj_data = {
        "atoms": [("Ar", "LJ")],
    }
    atomtypes = ["LJ"]
    LJ_type = Molecule_Type(mock_lj_data, atomtypes=atomtypes)
    
    # Load coordinates
    box = load_coords("SimpleLJ.clssy", [LJ_type])
    
    # Define LJ Forcefield
    lj_ff = LJ_Cut(nAtomTypes=1)
    lj_ff.rMin = np.zeros(1)
    lj_ff.rMinTable = np.zeros((1,1))
    box.EFunc.append(lj_ff)
    
    # Compute initial energy
    success = box.compute_energy()
    assert success, "Energy computation failed"
    print(f"Initial Energy: {box.ETotal}")
    
    # Create Monte Carlo simulation
    sim = SimMonteCarlo(
        nCycles=10,           # 10 cycles
        nMoves=100,           # 100 moves per cycle  
        screenfreq=5,         # Screen output every 5 cycles
        configfreq=10,        # Configuration output every 10 cycles
        energyCheck=10        # Energy check every 10 cycles
    )
    
    # Set up simulation components
    sim.BoxList = [box]  # Add our box to the simulation
    
    # Set up sampling method (Metropolis)
    sampling = Metropolis()
    sim.Sampling = sampling
    
    # Set up moves (molecular translation)
    transMove = MolTranslate([box])
    sim.Moves = [transMove]
    
    # Set move weights (equal probability for all moves)
    sim.moveweights = [1.0] * len(sim.Moves)
    
    try:
        # Run the simulation
        print("Starting Sim_MonteCarlo simulation...")
        sim.run_monte_carlo()
        print("Sim_MonteCarlo simulation completed successfully!")
        
        # Final energy check
        success = box.compute_energy()
        assert success, "Final energy computation failed"
        print(f"Final Energy: {box.ETotal}")
        
    except Exception as e:
        print(f"Error during Sim_MonteCarlo simulation: {e}")
        raise

#-------------------------------

if __name__ == "__main__":
    print("Running LJ Stack tests...")
    
    # Run the LJ stack test
    test_lj_stack()
    
    # Run the Sim_MonteCarlo test
    test_sim_monte_carlo()
    
    print("All tests completed.")
