#!/usr/bin/env python3
"""
Test runner for MC_Move_MolTranslation unit tests

This script provides a convenient way to run the MolTranslation unit tests
with proper setup and error handling.
"""

import sys
import os
import unittest

# Add the current directory to Python path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

def run_mol_translation_tests():
    """Run all MolTranslation unit tests"""
    
    print("="*60)
    print("Running MolTranslation Unit Tests")
    print("="*60)
    
    try:
        # Import the test module (using simplified version that works)
        from test_MC_Move_MolTranslation_simple import TestMolTranslateSimple
        
        # Create test suite
        loader = unittest.TestLoader()
        suite = unittest.TestSuite()
        
        # Add test cases
        suite.addTests(loader.loadTestsFromTestCase(TestMolTranslateSimple))
        
        # Run tests
        runner = unittest.TextTestRunner(verbosity=2)
        result = runner.run(suite)
        
        # Print summary
        print("\n" + "="*60)
        print("Test Summary:")
        print(f"Tests run: {result.testsRun}")
        print(f"Failures: {len(result.failures)}")
        print(f"Errors: {len(result.errors)}")
        if result.testsRun > 0:
            success_rate = (result.testsRun - len(result.failures) - len(result.errors)) / result.testsRun * 100
            print(f"Success rate: {success_rate:.1f}%")
        print("="*60)
        
        # Return success status
        return len(result.failures) == 0 and len(result.errors) == 0
        
    except ImportError as e:
        print(f"Error importing test module: {e}")
        print("Make sure all dependencies are available")
        return False
    except Exception as e:
        print(f"Unexpected error: {e}")
        return False

def run_specific_test(test_name):
    """Run a specific test case"""
    
    print(f"Running specific test: {test_name}")
    print("="*40)
    
    try:
        # Import the test module
        import test_MC_Move_MolTranslation_simple
        
        # Create test suite for specific test
        suite = unittest.TestLoader().loadTestsFromName(test_name, test_MC_Move_MolTranslation_simple)
        
        # Run test
        runner = unittest.TextTestRunner(verbosity=2)
        result = runner.run(suite)
        
        return len(result.failures) == 0 and len(result.errors) == 0
        
    except Exception as e:
        print(f"Error running test {test_name}: {e}")
        return False

def list_available_tests():
    """List all available tests"""
    print("Available tests:")
    print("-" * 40)
    
    try:
        from test_MC_Move_MolTranslation_simple import TestMolTranslateSimple
        
        loader = unittest.TestLoader()
        test_names = loader.getTestCaseNames(TestMolTranslateSimple)
        
        for test_name in test_names:
            print(f"  TestMolTranslateSimple.{test_name}")
            
    except Exception as e:
        print(f"Error listing tests: {e}")

if __name__ == "__main__":
    # Parse command line arguments
    if len(sys.argv) > 1:
        arg = sys.argv[1]
        if arg == "--list" or arg == "-l":
            list_available_tests()
            success = True
        else:
            success = run_specific_test(arg)
    else:
        success = run_mol_translation_tests()
    
    # Exit with appropriate code
    sys.exit(0 if success else 1) 