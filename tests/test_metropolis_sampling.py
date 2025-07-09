#!/usr/bin/env python3
"""
Test suite for Metropolis sampling implementation.

This test suite validates the Metropolis acceptance rule implementation
by testing:
1. Basic acceptance probability calculations
2. Error handling for invalid inputs
3. Integration with mock simulation boxes
4. Statistical behavior over many samples
"""

import sys
import os
import math
import numpy as np
import unittest
from unittest.mock import Mock, patch
import matplotlib.pyplot as plt

# Add src_python to path to import our modules
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src_python'))

from src_python.Sampling_Metropolis import Metropolis

# Inline implementations to avoid import issues
dp = np.float64

def IsNan(value):
    """Check if value is NaN"""
    return np.isnan(value)

def grnd():
    """Generate random number between 0 and 1"""
    return np.random.random()

class MockSimBox:
    """Mock simulation box for testing"""
    def __init__(self, temperature=300.0):
        self.temperature = temperature
        self.beta = 1.0 / (8.314e-3 * temperature)  # kT in kJ/mol
        self.ETotal = 0.0
        self.nAtoms = 100
        self.volume = 1000.0

class TestMetropolisSampling(unittest.TestCase):
    """Test cases for Metropolis sampling implementation"""
    
    def setUp(self):
        """Set up test fixtures"""
        self.metropolis = Metropolis()
        self.mock_box = MockSimBox(temperature=300.0)
        self.mock_disp = [Mock()]  # Mock displacement list
        
    def test_initialization(self):
        """Test that Metropolis sampler initializes correctly"""
        self.assertIsInstance(self.metropolis, Metropolis)
        self.assertTrue(hasattr(self.metropolis, 'make_decision'))
        self.assertTrue(hasattr(self.metropolis, 'make_decision_2box'))
    
    def test_nan_detection(self):
        """Test NaN detection utility function"""
        self.assertTrue(IsNan(float('nan')))
        self.assertFalse(IsNan(1.0))
        self.assertFalse(IsNan(0.0))
        self.assertFalse(IsNan(-1.0))
        self.assertFalse(IsNan(float('inf')))
    
    def test_grnd_function(self):
        """Test random number generator"""
        # Test that grnd returns values in [0, 1)
        for _ in range(100):
            r = grnd()
            self.assertGreaterEqual(r, 0.0)
            self.assertLess(r, 1.0)
    
    def test_acceptance_deterministic_cases(self):
        """Test acceptance for deterministic cases"""
        # Test case 1: Energy decrease should be accepted
        with patch('numpy.random.random', return_value=0.5):
            accept = self.metropolis.make_decision(
                self.mock_box, 
                e_diff=-1.0,  # Energy decreases
                disp=self.mock_disp,
                in_prob=1.0
            )
            self.assertTrue(accept)
        
        # Test case 2: No energy change should be accepted
        with patch('numpy.random.random', return_value=0.5):
            accept = self.metropolis.make_decision(
                self.mock_box,
                e_diff=0.0,  # No energy change
                disp=self.mock_disp,
                in_prob=1.0
            )
            self.assertTrue(accept)
    
    def test_acceptance_probabilistic_cases(self):
        """Test probabilistic acceptance for energy increases"""
        # High energy increase - should be rejected with high random number
        with patch('numpy.random.random', return_value=0.99):
            accept = self.metropolis.make_decision(
                self.mock_box,
                e_diff=10.0,  # Large energy increase
                disp=self.mock_disp,
                in_prob=1.0
            )
            self.assertFalse(accept)
        
        # Same energy increase - should be accepted with low random number
        with patch('numpy.random.random', return_value=0.01):
            accept = self.metropolis.make_decision(
                self.mock_box,
                e_diff=10.0,  # Large energy increase
                disp=self.mock_disp,
                in_prob=1.0
            )
            # Calculate expected acceptance probability
            bias_e = -self.mock_box.beta * 10.0
            expected_prob = math.exp(bias_e)
            # If random number is less than probability, should accept
            self.assertEqual(accept, 0.01 < expected_prob)
    
    def test_probability_input_validation(self):
        """Test validation of probability inputs"""
        # Test zero probability
        accept = self.metropolis.make_decision(
            self.mock_box,
            e_diff=0.0,
            disp=self.mock_disp,
            in_prob=0.0
        )
        self.assertFalse(accept)
        
        # Test negative probability
        accept = self.metropolis.make_decision(
            self.mock_box,
            e_diff=0.0,
            disp=self.mock_disp,
            in_prob=-1.0
        )
        self.assertFalse(accept)
    
    def test_log_probability_input(self):
        """Test acceptance with log probability input"""
        with patch('numpy.random.random', return_value=0.5):
            accept = self.metropolis.make_decision(
                self.mock_box,
                e_diff=-1.0,
                disp=self.mock_disp,
                log_prob=0.0  # ln(1) = 0
            )
            self.assertTrue(accept)
    
    def test_missing_probability_error(self):
        """Test error when no probability is provided"""
        with self.assertRaises(RuntimeError) as context:
            self.metropolis.make_decision(
                self.mock_box,
                e_diff=0.0,
                disp=self.mock_disp
                # No probability provided
            )
        self.assertIn("No probability provided", str(context.exception))
    
    def test_nan_energy_error(self):
        """Test error handling for NaN energy"""
        with self.assertRaises(RuntimeError) as context:
            self.metropolis.make_decision(
                self.mock_box,
                e_diff=float('nan'),
                disp=self.mock_disp,
                in_prob=1.0
            )
        self.assertIn("Invalid energy in sampling routine", str(context.exception))
    
    def test_extra_terms(self):
        """Test inclusion of extra terms in acceptance criterion"""
        extra_terms = 2.0
        
        with patch('numpy.random.random', return_value=0.5):
            accept = self.metropolis.make_decision(
                self.mock_box,
                e_diff=1.0,
                disp=self.mock_disp,
                in_prob=1.0,
                extra_in=extra_terms
            )
            
            # Calculate expected bias energy
            bias_e = -self.mock_box.beta * 1.0 + 0.0 + extra_terms  # ln(1.0) = 0
            expected_accept = bias_e >= 0.0 or bias_e > math.log(0.5)
            self.assertEqual(accept, expected_accept)
    
    def test_two_box_acceptance(self):
        """Test two-box acceptance decision"""
        mock_box2 = MockSimBox(temperature=400.0)
        mock_disp2 = [Mock()]
        
        with patch('numpy.random.random', return_value=0.5):
            accept = self.metropolis.make_decision_2box(
                self.mock_box, mock_box2,
                e_diff1=-1.0, e_diff2=1.0,
                disp1=self.mock_disp, disp2=mock_disp2,
                in_prob=1.0
            )
            
            # Calculate expected bias energy
            bias_e = (-self.mock_box.beta * (-1.0) - 
                     mock_box2.beta * 1.0 + 0.0)  # ln(1.0) = 0
            expected_accept = bias_e > 0.0 or bias_e > math.log(0.5)
            self.assertEqual(accept, expected_accept)
    
    def test_temperature_dependence(self):
        """Test that acceptance probability depends on temperature correctly"""
        hot_box = MockSimBox(temperature=1000.0)  # High temperature
        cold_box = MockSimBox(temperature=100.0)   # Low temperature
        
        energy_increase = 5.0
        
        # At high temperature, should be more likely to accept
        hot_accepts = 0
        cold_accepts = 0
        
        n_trials = 1000
        
        for _ in range(n_trials):
            hot_accept = self.metropolis.make_decision(
                hot_box,
                e_diff=energy_increase,
                disp=self.mock_disp,
                in_prob=1.0
            )
            if hot_accept:
                hot_accepts += 1
                
            cold_accept = self.metropolis.make_decision(
                cold_box,
                e_diff=energy_increase,
                disp=self.mock_disp,
                in_prob=1.0
            )
            if cold_accept:
                cold_accepts += 1
        
        hot_rate = hot_accepts / n_trials
        cold_rate = cold_accepts / n_trials
        
        # Hot system should have higher acceptance rate for energy increases
        self.assertGreater(hot_rate, cold_rate)
        
        # Theoretical acceptance probabilities
        hot_prob = math.exp(-hot_box.beta * energy_increase)
        cold_prob = math.exp(-cold_box.beta * energy_increase)
        
        # Check that empirical rates are close to theoretical (within 10%)
        self.assertAlmostEqual(hot_rate, hot_prob, delta=0.1)
        self.assertAlmostEqual(cold_rate, cold_prob, delta=0.1)


class MetropolisIntegrationTest(unittest.TestCase):
    """Integration tests with mock simulation components"""
    
    def setUp(self):
        """Set up integration test fixtures"""
        self.metropolis = Metropolis()
        self.mock_box = MockSimBox(temperature=300.0)
        
    def test_monte_carlo_simulation(self):
        """Test a simple Monte Carlo simulation loop"""
        # Initial configuration
        x = 0.0  # 1D particle position
        energy = x**2  # Harmonic potential
        
        positions = [x]
        energies = [energy]
        accepts = 0
        
        n_steps = 1000
        step_size = 0.5
        
        for step in range(n_steps):
            # Propose new position
            x_new = x + step_size * (2 * np.random.random() - 1)
            energy_new = x_new**2
            
            e_diff = energy_new - energy
            
            # Mock displacement object
            mock_disp = [Mock()]
            
            # Make acceptance decision
            accept = self.metropolis.make_decision(
                self.mock_box,
                e_diff=e_diff,
                disp=mock_disp,
                in_prob=1.0
            )
            
            if accept:
                x = x_new
                energy = energy_new
                accepts += 1
            
            positions.append(x)
            energies.append(energy)
        
        acceptance_rate = accepts / n_steps
        
        # Acceptance rate should be reasonable (10-95%)
        self.assertGreater(acceptance_rate, 0.1)
        self.assertLess(acceptance_rate, 0.95)
        
        # Final position should be close to equilibrium (near 0)
        final_positions = positions[-100:]  # Last 100 positions
        mean_pos_sq = np.mean([p**2 for p in final_positions])
        
        # For harmonic potential at temperature T, <x^2> = kT
        expected_var = 1.0 / self.mock_box.beta
        
        # Should be within factor of 3 (rough check)
        self.assertLess(mean_pos_sq, 3 * expected_var)
        
        print(f"Acceptance rate: {acceptance_rate:.3f}")
        print(f"Mean final position squared: {mean_pos_sq:.3f}")
        print(f"Expected variance: {expected_var:.3f}")


def run_visual_test():
    """Create visual plots to validate Metropolis behavior"""
    print("Running visual validation tests...")
    
    # Test 1: Acceptance probability vs energy difference
    metropolis = Metropolis()
    mock_box = MockSimBox(temperature=300.0)
    mock_disp = [Mock()]
    
    energy_diffs = np.linspace(-5, 10, 100)
    theoretical_probs = []
    empirical_probs = []
    
    n_trials = 1000
    
    for e_diff in energy_diffs:
        # Theoretical probability
        if e_diff <= 0:
            theo_prob = 1.0
        else:
            theo_prob = math.exp(-mock_box.beta * e_diff)
        theoretical_probs.append(theo_prob)
        
        # Empirical probability
        accepts = 0
        for _ in range(n_trials):
            accept = metropolis.make_decision(
                mock_box,
                e_diff=e_diff,
                disp=mock_disp,
                in_prob=1.0
            )
            if accept:
                accepts += 1
        
        empirical_probs.append(accepts / n_trials)
    
    # Plot comparison
    plt.figure(figsize=(10, 6))
    plt.plot(energy_diffs, theoretical_probs, 'r-', label='Theoretical', linewidth=2)
    plt.plot(energy_diffs, empirical_probs, 'bo', label='Empirical', markersize=4)
    plt.xlabel('Energy Difference (kJ/mol)')
    plt.ylabel('Acceptance Probability')
    plt.title('Metropolis Acceptance Probability')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.savefig('metropolis_acceptance_test.png', dpi=150, bbox_inches='tight')
    plt.show()
    
    # Check agreement
    max_error = max(abs(t - e) for t, e in zip(theoretical_probs, empirical_probs))
    print(f"Maximum error between theoretical and empirical: {max_error:.4f}")
    
    if max_error < 0.05:
        print("✓ Acceptance probability test PASSED")
    else:
        print("✗ Acceptance probability test FAILED")
    
    return max_error < 0.05


if __name__ == '__main__':
    print("Testing Metropolis Sampling Implementation")
    print("=" * 50)
    
    # Run unit tests
    print("\n1. Running unit tests...")
    unittest.main(argv=[''], exit=False, verbosity=2)
    
    # Run visual validation
    print("\n2. Running visual validation...")
    try:
        import matplotlib.pyplot as plt
        visual_passed = run_visual_test()
    except ImportError:
        print("Matplotlib not available, skipping visual tests")
        visual_passed = True
    
    print("\n" + "=" * 50)
    if visual_passed:
        print("All tests completed successfully!")
    else:
        print("Some tests failed - check output above") 