# Metropolis Sampling Tests for ClassyMC

This directory contains comprehensive tests for the Metropolis sampling implementation in the ClassyMC Python codebase. The tests are designed to validate the correctness of the newly converted Monte Carlo acceptance rules.

## Test Files Overview

### 1. `test_metropolis_simple.py` - Self-Contained Unit Tests
**Purpose**: Basic functionality tests that run without external dependencies.

**What it tests**:
- ✅ Metropolis acceptance probability calculations
- ✅ Temperature dependence of acceptance rates  
- ✅ Error handling (NaN energies, missing probabilities)
- ✅ Two-box system acceptance decisions
- ✅ Detailed balance validation
- ✅ Integration with simple Monte Carlo simulation loop

**How to run**:
```bash
python test_metropolis_simple.py
```

**Expected output**: All unit tests should pass, plus detailed validation showing:
- Acceptance probabilities match theoretical predictions (< 5% error)
- Higher temperatures give higher acceptance rates for energy increases
- Detailed balance is satisfied (< 10% error)

### 2. `test_metropolis_sampling.py` - Comprehensive Test Suite  
**Purpose**: Advanced tests with visual validation and statistical analysis.

**What it tests**:
- 📊 Visual plots of acceptance probability vs energy difference
- 📈 Statistical validation over large sample sizes
- 🔬 Monte Carlo simulation with harmonic potential
- 🌡️ Temperature effects on equilibrium distributions

**Requirements**: 
- matplotlib (for plots)
- Access to `src_python/` modules

**How to run**:
```bash
python test_metropolis_sampling.py
```

### 3. `demo_metropolis_integration.py` - Integration Demo
**Purpose**: Demonstrate Metropolis sampling with actual ClassyMC components.

**What it shows**:
- 🔧 Integration with SimpleBox simulation boxes
- ⚛️ Hard sphere particle interactions
- 🎯 Realistic Monte Carlo simulation workflow
- 📊 Temperature scaling effects

**How to run**:
```bash
python demo_metropolis_integration.py
```

## Test Validation Criteria

### Core Metropolis Algorithm
The tests validate that the implementation correctly implements:

```
P(accept) = min(1, exp(-β * ΔE + ln(P_forward/P_backward) + extra_terms))
```

Where:
- `β = 1/(k_B * T)` is the inverse temperature
- `ΔE` is the energy difference
- `P_forward/P_backward` accounts for non-symmetric moves
- `extra_terms` include chemical potential, pressure terms, etc.

### Expected Behaviors

1. **Energy Decreases**: Always accepted (P = 1.0)
2. **Energy Increases**: Accepted with probability `exp(-β * ΔE)`
3. **Temperature Scaling**: Higher T → higher acceptance rates for ΔE > 0
4. **Detailed Balance**: `P(i→j)/P(j→i) = exp(-β(E_j - E_i))`

### Acceptance Rate Guidelines
- **Optimal range**: 20-80% for most move types
- **Too high (>90%)**: Moves too small, poor sampling efficiency
- **Too low (<10%)**: Moves too large, poor acceptance

## Common Issues and Troubleshooting

### Import Errors
If you see import errors for ClassyMC modules:
1. Ensure you're running from the correct directory
2. Check that `src_python/` contains the required modules
3. The tests include fallback mock implementations

### Statistical Fluctuations
- Acceptance rate tests use statistical sampling
- Small deviations (<5%) from theoretical values are normal
- Rerun tests if marginal failures occur

### Temperature Units
- Default temperature unit: Kelvin (K)
- Energy unit: kJ/mol (consistent with `β = 1/(R*T)` where R = 8.314 J/mol/K)

## Test Results Interpretation

### Successful Test Output
```
Testing Metropolis Sampling Implementation
==================================================

1. Running unit tests...
test_energy_decrease_accepted ... ok
test_temperature_effect ... ok
test_detailed_balance ... ok
...

2. Running detailed tests...
ΔE =  0.0: Empirical = 1.000, Theoretical = 1.000, Error = 0.000
ΔE =  1.0: Empirical = 0.698, Theoretical = 0.700, Error = 0.002
...
✓ Acceptance probability test PASSED
✓ Detailed balance test PASSED

🎉 All Metropolis sampling tests completed successfully!
```

### Failure Indicators
- ❌ **Large acceptance probability errors (>5%)**
  - Check random number generator
  - Verify temperature/energy units
  - Ensure sufficient statistical sampling

- ❌ **Detailed balance violations (>10%)**
  - Check implementation of probability terms
  - Verify energy difference calculations
  - Look for systematic biases

- ❌ **Runtime errors**
  - Check for NaN/infinity handling
  - Verify all required parameters are provided
  - Ensure mock objects have required attributes

## Performance Benchmarks

On a typical modern CPU, expect:
- **Simple tests**: ~5-10 seconds
- **Comprehensive tests**: ~30-60 seconds  
- **Integration demo**: ~10-20 seconds

## Extending the Tests

To add new validation tests:

1. **Add test method** to `TestMetropolis` class:
```python
def test_your_feature(self):
    """Test description"""
    # Setup
    # Test logic
    # Assertions
```

2. **Add edge cases** to cover:
   - Extreme temperatures (very hot/cold)
   - Large energy changes
   - Zero/negative probabilities
   - Special move types

3. **Add integration scenarios**:
   - Different force fields
   - Various box types
   - Multi-component systems

## References

- **Metropolis Algorithm**: Metropolis et al. (1953) J. Chem. Phys. 21, 1087
- **Monte Carlo Methods**: Frenkel & Smit, "Understanding Molecular Simulation"
- **Detailed Balance**: Ensures correct equilibrium sampling in MC simulations

## Support

If tests fail or you encounter issues:
1. Check this README for common solutions
2. Verify your Python environment (numpy, unittest)
3. Ensure ClassyMC modules are properly converted from Fortran
4. Run individual test functions to isolate problems

The tests are designed to be robust and provide clear error messages to help diagnose issues with the Metropolis sampling implementation. 