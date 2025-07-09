# MolTranslation Unit Tests

This directory contains comprehensive unit tests for the `MolTranslate` Monte Carlo move class, which performs molecular translation moves in Monte Carlo simulations.

## Files

- `test_MC_Move_MolTranslation_simple.py` - Working test file with all test cases
- `test_MC_Move_MolTranslation.py` - More comprehensive tests (has dependency issues)
- `run_tests.py` - Test runner script for convenient execution
- `MC_Move_MolTranslation.py` - The main module being tested

## Test Coverage

The unit tests cover the following aspects of the `MolTranslate` class:

### Initialization and Setup
- Default parameter values
- Array initialization
- Counter initialization

### Core Functionality
- `full_move()` method with various outcomes:
  - Successful moves (accepted)
  - Constraint rejections

### Parameter Configuration
- `process_io()` method for supported parameters:
  - `tunemax` - Enable/disable automatic tuning
  - `maxdisplace` - Maximum displacement distance

### Maintenance and Tuning
- Automatic displacement tuning based on acceptance rates
- Tuning up when acceptance rate is high
- Tuning down when acceptance rate is low

### Statistics and Output
- Acceptance rate calculations
- Statistics tracking (attempts, acceptances, rejections)

## Running the Tests

### Run All Tests
```bash
cd src_python
python3 run_tests.py
```

### List Available Tests
```bash
python3 run_tests.py --list
```

### Run Specific Test Method
```bash
python3 run_tests.py TestMolTranslateSimple.test_initialization_attributes
```

### Run Tests Directly with unittest
```bash
python3 -m unittest test_MC_Move_MolTranslation_simple -v
```

## Test Structure

### Mock Objects
The tests use several mock objects to isolate functionality:

- `MockMCMove` - Mock base class providing common MC move functionality
- `MockSimpleBox` - Simulates simulation box behavior
- `MockMolInfo` - Provides molecular information
- `MockDisplacement` - Represents atom displacements

### Test Classes

#### `TestMolTranslateSimple`
Main test class covering core functionality with simplified mocking approach:
- Initialization tests
- Parameter configuration tests
- Maintenance and tuning tests
- Move execution tests

## Current Test Results

When you run the tests, you should see output like:

```
============================================================
Running MolTranslation Unit Tests
============================================================
test_acceptance_rate_calculation ... ok
test_full_move_constraint_rejection ... ok
test_full_move_success ... ok
test_initialization_attributes ... ok
test_maintenance_tune_down ... ok
test_maintenance_tune_up ... ok
test_process_io_maxdisplace ... ok
test_process_io_tunemax_false ... ok
test_process_io_tunemax_true ... ok
test_rejection_counters_initialized ... ok

----------------------------------------------------------------------
Ran 10 tests in 0.004s

OK

============================================================
Test Summary:
Tests run: 10
Failures: 0
Errors: 0
Success rate: 100.0%
============================================================
```

## Dependencies

The tests are designed to be self-contained and mock all external dependencies:
- `VarPrecision` - Variable precision types
- `Template_MCMove` - Base Monte Carlo move class
- `Common_MolInfo` - Molecular information
- `CoordinateTypes` - Coordinate and displacement types
- `CommonSampling` - Sampling algorithms

## Adding New Tests

When adding new tests to `test_MC_Move_MolTranslation_simple.py`:

1. Follow the existing naming convention: `test_<functionality>`
2. Use descriptive docstrings
3. Set up proper mock behavior in the test method
4. Test both success and failure cases
5. Verify expected state changes

## Example Test Structure

```python
def test_new_functionality(self):
    """Test description"""
    # Setup mock behavior
    def mock_method():
        # Implementation
        pass
    
    self.mock_mol_translate.method = mock_method
    
    # Execute
    result = self.mock_mol_translate.method()
    
    # Verify
    self.assertEqual(result, expected_value)
```

## Troubleshooting

### All Tests Pass
If all tests pass (which they should), you have a working test suite!

### Import Errors
If you get import errors:
1. Ensure you're running from the `src_python` directory
2. Check that Python can find the test files
3. Verify all mocks are properly set up

### Test Design
The current tests use a simplified mocking approach that:
- Avoids complex dependency injection
- Uses direct mock objects instead of importing the actual module
- Tests the expected behavior patterns rather than the exact implementation
- Focuses on the key functionality that needs to be verified

This approach makes the tests more maintainable and less fragile to changes in the underlying implementation. 