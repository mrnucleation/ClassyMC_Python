#!/bin/bash

# Metropolis Sampling Test Runner
# Runs all Metropolis sampling tests for ClassyMC

echo "========================================================"
echo "         ClassyMC Metropolis Sampling Test Suite"
echo "========================================================"
echo ""

# Check if we're in the right directory
if [ ! -f "test_metropolis_simple.py" ]; then
    echo "❌ Error: test_metropolis_simple.py not found"
    echo "Please run this script from the test directory"
    exit 1
fi

# Check Python availability
if ! command -v python3 &> /dev/null; then
    echo "❌ Error: python3 not found"
    echo "Please install Python 3.x"
    exit 1
fi

# Track test results
TESTS_PASSED=0
TESTS_FAILED=0

echo "1. Running simple Metropolis tests..."
echo "----------------------------------------------------"
if python3 test_metropolis_simple.py; then
    echo "✅ Simple tests PASSED"
    TESTS_PASSED=$((TESTS_PASSED + 1))
else
    echo "❌ Simple tests FAILED"
    TESTS_FAILED=$((TESTS_FAILED + 1))
fi
echo ""

echo "2. Running comprehensive Metropolis tests..."
echo "----------------------------------------------------"
if python3 test_metropolis_sampling.py; then
    echo "✅ Comprehensive tests PASSED"
    TESTS_PASSED=$((TESTS_PASSED + 1))
else
    echo "❌ Comprehensive tests FAILED"
    TESTS_FAILED=$((TESTS_FAILED + 1))
fi
echo ""

echo "3. Running integration demonstration..."
echo "----------------------------------------------------"
if python3 demo_metropolis_integration.py; then
    echo "✅ Integration demo PASSED"
    TESTS_PASSED=$((TESTS_PASSED + 1))
else
    echo "❌ Integration demo FAILED"
    TESTS_FAILED=$((TESTS_FAILED + 1))
fi
echo ""

# Summary
echo "========================================================"
echo "                    Test Summary"
echo "========================================================"
echo "Tests passed: $TESTS_PASSED"
echo "Tests failed: $TESTS_FAILED"
echo ""

if [ $TESTS_FAILED -eq 0 ]; then
    echo "🎉 All Metropolis sampling tests completed successfully!"
    echo ""
    echo "Your Metropolis implementation appears to be working correctly."
    echo "The acceptance probabilities match theoretical predictions,"
    echo "detailed balance is satisfied, and integration with simulation"
    echo "components is functioning properly."
    exit 0
else
    echo "⚠️  Some tests failed. Please check the output above for details."
    echo ""
    echo "Common issues:"
    echo "- Import errors: Check that src_python/ modules are available"
    echo "- Statistical fluctuations: Rerun tests to see if errors persist"  
    echo "- Environment issues: Ensure numpy and other dependencies are installed"
    echo ""
    echo "See README_Metropolis_Tests.md for troubleshooting guidance."
    exit 1
fi 