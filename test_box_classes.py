#!/usr/bin/env python3
"""
Test script for box classes - CubeBox and SimpleBox
Tests loading dimensions and coordinates functionality
"""

import numpy as np
import tempfile
import os
from src_python.Box_CubeBox import CubeBox
from src_python.Box_SimpleBox import SimpleBox
from src_python.Script_LoadCoordinates import load_coords
from src_python.Molecule_Definition import Molecule_Type


mock_lj_data = {
    "atoms": [("Ar", "LJ")],
}
atomtypes = ["LJ"]
LJ_type = Molecule_Type(mock_lj_data, atomtypes=atomtypes)



def test_cube_box_load_dimension():
    """Test CubeBox load_dimension method"""
    print("Testing CubeBox load_dimension...")


    
    # Initialize cube box
    cube_box = CubeBox([LJ_type])
    
    # Test valid dimension line
    result = cube_box.load_dimension("boxlength 10.0")
    assert result, f"Expected True, got {result}"
    assert cube_box.boxL == 10.0, f"Expected boxL=10.0, got {cube_box.boxL}"
    assert cube_box.boxL2 == 5.0, f"Expected boxL2=5.0, got {cube_box.boxL2}"
    
    # Test invalid dimension line
    result = cube_box.load_dimension("boxlength invalid")
    assert result == False, f"Expected False, got {result}"
    
    print("✓ CubeBox load_dimension test passed")

def test_simple_box_load_dimension():
    """Test SimpleBox load_dimension method"""
    print("Testing SimpleBox load_dimension...")
    
    # Create mock molecular data
    mol_data = [
        {'nAtoms': 2, 'masses': [1.0, 1.0], 'ridgid': False},
        {'nAtoms': 3, 'masses': [1.0, 1.0, 1.0], 'ridgid': False}
    ]
    
    # Initialize simple box
    simple_box = SimpleBox(mol_data, NMolMin=[0, 0], NMolMax=[10, 10], NMol=[5, 5])
    
    # Test dimension line (should return 0 for simple box)
    result = simple_box.load_dimension("any line")
    assert result, f"Expected True, got {result}"
    
    print("✓ SimpleBox load_dimension test passed")

def test_cube_box_coordinates():
    """Test loading coordinates into CubeBox"""
    print("Testing CubeBox coordinate loading...")
    
    # Create test coordinate file content
    coord_content = """cube 10.0
NMol 2 3
NMax 5 5
NMin 0 0
1 1 1 1.0 2.0 3.0
1 1 2 4.0 5.0 6.0
2 1 1 7.0 8.0 9.0
2 1 2 10.0 11.0 12.0
2 1 3 13.0 14.0 15.0
"""
    
    # Write to temporary file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
        f.write(coord_content)
        temp_file = f.name
    
    try:
        # Load coordinates
        box = load_coords(temp_file, [LJ_type])
        
        # Verify it's a CubeBox
        assert isinstance(box, CubeBox), f"Expected CubeBox, got {type(box)}"
        
        # Check dimensions
        assert box.boxL == 10.0, f"Expected boxL=10.0, got {box.boxL}"
        assert box.boxL2 == 5.0, f"Expected boxL2=5.0, got {box.boxL2}"
        
        # Check molecule counts
        assert np.array_equal(box.NMol, [2, 3]), f"Expected NMol=[2, 3], got {box.NMol}"
        assert np.array_equal(box.NMolMax, [5, 5]), f"Expected NMolMax=[5, 5], got {box.NMolMax}"
        assert np.array_equal(box.NMolMin, [0, 0]), f"Expected NMolMin=[0, 0], got {box.NMolMin}"
        
        # Check coordinates
        expected_coords = np.array([
            [1.0, 2.0, 3.0],
            [4.0, 5.0, 6.0],
            [7.0, 8.0, 9.0],
            [10.0, 11.0, 12.0],
            [13.0, 14.0, 15.0]
        ])
        
        assert np.allclose(box.atoms[:5], expected_coords), "Coordinates don't match expected values"
        
        print("✓ CubeBox coordinate loading test passed")
        
    finally:
        # Clean up
        os.unlink(temp_file)

def test_simple_box_coordinates():
    """Test loading coordinates into SimpleBox"""
    print("Testing SimpleBox coordinate loading...")
    
    # Create test coordinate file content
    coord_content = """nobox
NMol 2 3
NMax 5 5
NMin 0 0
1 1 1 1.0 2.0 3.0
1 1 2 4.0 5.0 6.0
2 1 1 7.0 8.0 9.0
2 1 2 10.0 11.0 12.0
2 1 3 13.0 14.0 15.0
"""
    
    # Write to temporary file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
        f.write(coord_content)
        temp_file = f.name
    
    try:
        # Load coordinates
        box = load_coords(temp_file)
        
        # Verify it's a SimpleBox
        assert isinstance(box, SimpleBox), f"Expected SimpleBox, got {type(box)}"
        
        # Check molecule counts
        assert np.array_equal(box.NMol, [2, 3]), f"Expected NMol=[2, 3], got {box.NMol}"
        assert np.array_equal(box.NMolMax, [5, 5]), f"Expected NMolMax=[5, 5], got {box.NMolMax}"
        assert np.array_equal(box.NMolMin, [0, 0]), f"Expected NMolMin=[0, 0], got {box.NMolMin}"
        
        # Check coordinates
        expected_coords = np.array([
            [1.0, 2.0, 3.0],
            [4.0, 5.0, 6.0],
            [7.0, 8.0, 9.0],
            [10.0, 11.0, 12.0],
            [13.0, 14.0, 15.0]
        ])
        
        assert np.allclose(box.atoms[:5], expected_coords), "Coordinates don't match expected values"
        
        print("✓ SimpleBox coordinate loading test passed")
        
    finally:
        # Clean up
        os.unlink(temp_file)

def test_boundary_conditions():
    """Test boundary conditions for CubeBox"""
    print("Testing CubeBox boundary conditions...")
    
    # Create mock molecular data
    mol_data = [{'nAtoms': 1, 'masses': [1.0], 'ridgid': False}]
    
    # Initialize cube box with length 10.0
    cube_box = CubeBox(mol_data, NMolMin=[0], NMolMax=[10], NMol=[1])
    cube_box.boxL = 10.0
    cube_box.boxL2 = 5.0
    
    # Test coordinates within bounds
    coords = np.array([3.0, 4.0, 5.0])
    result = cube_box.boundary(coords)
    assert np.allclose(result, coords), "Coordinates within bounds should not change"
    
    # Test coordinates outside bounds (positive)
    coords = np.array([6.0, 7.0, 8.0])
    result = cube_box.boundary(coords)
    expected = np.array([-4.0, -3.0, -2.0])  # 6-10, 7-10, 8-10
    assert np.allclose(result, expected), f"Expected {expected}, got {result}"
    
    # Test coordinates outside bounds (negative)
    coords = np.array([-6.0, -7.0, -8.0])
    result = cube_box.boundary(coords)
    expected = np.array([4.0, 3.0, 2.0])  # -6+10, -7+10, -8+10
    assert np.allclose(result, expected), f"Expected {expected}, got {result}"
    
    print("✓ CubeBox boundary conditions test passed")

def test_box_dimensions():
    """Test getting dimensions for CubeBox"""
    print("Testing CubeBox dimensions...")
    
    # Create mock molecular data
    mol_data = [{'nAtoms': 1, 'masses': [1.0], 'ridgid': False}]
    
    # Initialize cube box with length 10.0
    cube_box = CubeBox(mol_data, NMolMin=[0], NMolMax=[10], NMol=[1])
    cube_box.boxL = 10.0
    cube_box.boxL2 = 5.0
    
    # Test get_dimensions
    dimensions = cube_box.get_dimensions()
    expected = [[-5.0, 5.0], [-5.0, 5.0], [-5.0, 5.0]]
    assert len(dimensions) == 3, f"Expected 3 dimensions, got {len(dimensions)}"
    for i, dim in enumerate(dimensions):
        assert dim == expected[i], f"Dimension {i}: expected {expected[i]}, got {dim}"
    
    print("✓ CubeBox dimensions test passed")

def main():
    """Run all tests"""
    print("Running box class tests...\n")
    
    try:
        test_cube_box_load_dimension()
        test_simple_box_load_dimension()
        test_cube_box_coordinates()
        test_simple_box_coordinates()
        test_boundary_conditions()
        test_box_dimensions()
        
        print("\nAll tests passed!")
        
    except Exception as e:
        print(f"\nTest failed: {e}")
        raise

if __name__ == "__main__":
    main()
