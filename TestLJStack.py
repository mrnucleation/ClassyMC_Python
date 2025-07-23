

from src_python.Script_LoadCoordinates import load_coords
from src_python.Molecule_Definition import Molecule_Type
from src_python.FF_LJ_Cut import LJ_Cut
import numpy as np


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
    
    


if __name__ == "__main__":
    print("Running LJ Stack tests...")
    
    # Run the LJ stack test
    test_lj_stack()
    
    print("LJ Stack tests completed.")
