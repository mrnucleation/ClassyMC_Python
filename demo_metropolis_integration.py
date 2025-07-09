#!/usr/bin/env python3
"""
Demonstration of Metropolis sampling integration with ClassyMC components.

This script shows how to set up and run a basic Monte Carlo simulation
using the Metropolis acceptance rule with real simulation components:
- SimpleBox simulation box
- Hard sphere or LJ force field  
- Molecular translation moves
- Energy calculations
"""

import sys
import os
import numpy as np
import math

# Add src_python to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src_python'))

def create_simple_system():
    """
    Create a simple test system with a few particles in a box.
    Returns the simulation box and required components.
    """
    # Mock molecular data for single atom molecules
    mol_data = [{
        'nAtoms': 1,
        'atomType': [0],  # All atoms are type 0
        'masses': [1.0],
        'ridgid': True
    }]
    
    # Create and initialize simulation box
    try:
        from Box_SimpleBox import SimpleBox
        box = SimpleBox(mol_data)
        
        # Set molecular bounds
        NMolMin = [0]     # Minimum 0 molecules
        NMolMax = [10]    # Maximum 10 molecules  
        NMol = [4]        # Current 4 molecules
        
        box.constructor(NMolMin, NMolMax, NMol)
        
        # Set up basic 3D positions for 4 single-atom molecules
        box.nDimension = 3
        positions = np.array([
            [0.0, 0.0, 0.0],    # Molecule 1
            [1.5, 0.0, 0.0],    # Molecule 2  
            [0.0, 1.5, 0.0],    # Molecule 3
            [1.5, 1.5, 0.0]     # Molecule 4
        ]).T
        
        # Set positions in the box
        for i in range(4):
            box.atoms[:, i] = positions[:, i]
        
        return box, mol_data
        
    except ImportError as e:
        print(f"Could not import simulation components: {e}")
        print("Using mock components instead...")
        return create_mock_system()

def create_mock_system():
    """Create a mock system for testing when real components aren't available"""
    
    class MockBox:
        def __init__(self):
            self.temperature = 300.0
            self.beta = 1.0 / (8.314e-3 * self.temperature)
            self.nAtoms = 4
            self.nMolTotal = 4
            self.atoms = np.array([
                [0.0, 1.5, 0.0, 1.5],
                [0.0, 0.0, 1.5, 1.5], 
                [0.0, 0.0, 0.0, 0.0]
            ])
            self.ETotal = 0.0
            self.volume = 8.0  # 2x2x2 box
            
        def boundary(self, rx, ry, rz):
            """No boundary conditions for simple test"""
            return rx, ry, rz
            
        def compute_energy(self):
            """Simple hard sphere energy - infinite if overlapping"""
            sigma = 1.0  # Hard sphere diameter
            
            for i in range(self.nAtoms-1):
                for j in range(i+1, self.nAtoms):
                    dx = self.atoms[0,i] - self.atoms[0,j]
                    dy = self.atoms[1,i] - self.atoms[1,j] 
                    dz = self.atoms[2,i] - self.atoms[2,j]
                    r = math.sqrt(dx*dx + dy*dy + dz*dz)
                    
                    if r < sigma:
                        self.ETotal = float('inf')
                        return False
                        
            self.ETotal = 0.0
            return True
            
        def compute_energy_delta(self, atom_idx, new_pos):
            """Compute energy change for moving one atom"""
            old_pos = self.atoms[:, atom_idx].copy()
            
            # Temporarily move atom
            self.atoms[:, atom_idx] = new_pos
            
            # Check for overlaps
            sigma = 1.0
            for j in range(self.nAtoms):
                if j == atom_idx:
                    continue
                    
                dx = new_pos[0] - self.atoms[0,j]
                dy = new_pos[1] - self.atoms[1,j]
                dz = new_pos[2] - self.atoms[2,j]
                r = math.sqrt(dx*dx + dy*dy + dz*dz)
                
                if r < sigma:
                    # Restore position
                    self.atoms[:, atom_idx] = old_pos
                    return float('inf'), False
            
            # Restore position  
            self.atoms[:, atom_idx] = old_pos
            return 0.0, True  # No energy change if no overlap
    
    mol_data = [{'nAtoms': 1, 'atomType': [0]}]
    return MockBox(), mol_data

def create_metropolis_sampler():
    """Create and initialize Metropolis sampler"""
    
    try:
        from Sampling_Metropolis import Metropolis
        from CommonSampling import set_sampling
        
        metropolis = Metropolis()
        set_sampling(metropolis)
        return metropolis
        
    except ImportError:
        print("Could not import Metropolis sampler, using inline version...")
        
        # Inline Metropolis implementation
        class Metropolis:
            def __init__(self):
                self.accepts = 0
                self.attempts = 0
                
            def make_decision(self, trial_box, e_diff, disp, in_prob=1.0, log_prob=None, extra_in=None):
                self.attempts += 1
                
                if np.isnan(e_diff) or e_diff == float('inf'):
                    return False
                    
                if e_diff <= 0.0:
                    self.accepts += 1
                    return True
                    
                # Metropolis criterion for energy increase
                prob = math.exp(-trial_box.beta * e_diff)
                if np.random.random() < prob:
                    self.accepts += 1
                    return True
                    
                return False
                
            def get_acceptance_rate(self):
                return self.accepts / max(1, self.attempts)
        
        return Metropolis()

def run_monte_carlo_simulation(box, metropolis, n_steps=1000):
    """
    Run a simple Monte Carlo simulation with Metropolis sampling.
    
    Args:
        box: Simulation box
        metropolis: Metropolis sampler
        n_steps: Number of MC steps to run
    """
    print(f"\nRunning {n_steps} Monte Carlo steps...")
    print("Initial configuration:")
    print("Atom positions:")
    for i in range(box.nAtoms):
        pos = box.atoms[:, i]
        print(f"  Atom {i}: ({pos[0]:.2f}, {pos[1]:.2f}, {pos[2]:.2f})")
    
    # Initial energy
    box.compute_energy()
    print(f"Initial energy: {box.ETotal}")
    
    max_displacement = 0.1
    accepts = 0
    
    for step in range(n_steps):
        # Select random atom to move
        atom_idx = np.random.randint(box.nAtoms)
        
        # Propose random displacement
        displacement = max_displacement * (2 * np.random.random(3) - 1)
        old_pos = box.atoms[:, atom_idx].copy()
        new_pos = old_pos + displacement
        
        # Calculate energy change
        e_diff, accept_energy = box.compute_energy_delta(atom_idx, new_pos)
        
        if not accept_energy:
            continue  # Skip if energy calculation failed
            
        # Mock displacement object
        mock_disp = [type('MockDisp', (), {
            'atmIndx': atom_idx,
            'x_new': new_pos[0],
            'y_new': new_pos[1], 
            'z_new': new_pos[2]
        })()]
        
        # Make acceptance decision
        accept = metropolis.make_decision(box, e_diff, mock_disp)
        
        if accept:
            # Accept move
            box.atoms[:, atom_idx] = new_pos
            box.ETotal += e_diff
            accepts += 1
            
        # Print progress
        if (step + 1) % (n_steps // 10) == 0:
            acceptance_rate = accepts / (step + 1)
            print(f"Step {step+1:4d}: Acceptance rate = {acceptance_rate:.3f}")
    
    final_acceptance_rate = accepts / n_steps
    
    print(f"\nSimulation completed!")
    print(f"Final acceptance rate: {final_acceptance_rate:.3f}")
    print("Final configuration:")
    print("Atom positions:")
    for i in range(box.nAtoms):
        pos = box.atoms[:, i]
        print(f"  Atom {i}: ({pos[0]:.2f}, {pos[1]:.2f}, {pos[2]:.2f})")
    
    box.compute_energy()
    print(f"Final energy: {box.ETotal}")
    
    return final_acceptance_rate

def test_temperature_effects():
    """Test how temperature affects acceptance rates"""
    print("\n" + "="*60)
    print("Testing temperature effects on Metropolis sampling")
    print("="*60)
    
    temperatures = [100.0, 300.0, 500.0, 1000.0]
    n_steps = 500
    
    for temp in temperatures:
        print(f"\nTesting at T = {temp} K:")
        
        box, mol_data = create_mock_system()
        box.temperature = temp
        box.beta = 1.0 / (8.314e-3 * temp)
        
        metropolis = create_metropolis_sampler()
        acceptance_rate = run_monte_carlo_simulation(box, metropolis, n_steps)
        
        print(f"T = {temp:6.0f} K: Acceptance rate = {acceptance_rate:.3f}")

def validate_detailed_balance():
    """Validate that detailed balance is satisfied"""
    print("\n" + "="*60)
    print("Validating detailed balance")
    print("="*60)
    
    # Create simple 2-state system
    class TwoStateSystem:
        def __init__(self, temp=300.0):
            self.temperature = temp
            self.beta = 1.0 / (8.314e-3 * temp)
            self.state = 0  # 0 or 1
            self.energies = [0.0, 5.0]  # State 0 has energy 0, state 1 has energy 5
            
        def get_energy_diff(self, new_state):
            return self.energies[new_state] - self.energies[self.state]
            
        def set_state(self, new_state):
            self.state = new_state
    
    system = TwoStateSystem()
    metropolis = create_metropolis_sampler()
    
    # Test transitions 0→1 and 1→0
    n_trials = 10000
    
    # Transition 0→1
    system.state = 0
    accepts_01 = 0
    for _ in range(n_trials):
        e_diff = system.get_energy_diff(1)
        mock_disp = [type('MockDisp', (), {})()]
        accept = metropolis.make_decision(system, e_diff, mock_disp)
        if accept:
            accepts_01 += 1
    
    # Transition 1→0  
    system.state = 1
    accepts_10 = 0
    for _ in range(n_trials):
        e_diff = system.get_energy_diff(0)
        mock_disp = [type('MockDisp', (), {})()]
        accept = metropolis.make_decision(system, e_diff, mock_disp)
        if accept:
            accepts_10 += 1
    
    prob_01 = accepts_01 / n_trials
    prob_10 = accepts_10 / n_trials
    
    print(f"P(0→1) = {prob_01:.4f}")
    print(f"P(1→0) = {prob_10:.4f}")
    
    if prob_10 > 0:
        empirical_ratio = prob_01 / prob_10
        theoretical_ratio = math.exp(-system.beta * (system.energies[1] - system.energies[0]))
        
        print(f"Empirical ratio P(0→1)/P(1→0) = {empirical_ratio:.4f}")
        print(f"Theoretical ratio = exp(-βΔE) = {theoretical_ratio:.4f}")
        print(f"Relative error = {abs(empirical_ratio - theoretical_ratio)/theoretical_ratio:.4f}")
        
        if abs(empirical_ratio - theoretical_ratio)/theoretical_ratio < 0.1:
            print("✓ Detailed balance test PASSED")
        else:
            print("✗ Detailed balance test FAILED")
    else:
        print("Cannot test detailed balance (prob_10 = 0)")

def main():
    """Main demonstration function"""
    print("ClassyMC Metropolis Sampling Demonstration")
    print("="*60)
    
    try:
        # Create system
        print("\n1. Creating simulation system...")
        box, mol_data = create_simple_system()
        print(f"✓ Created system with {box.nAtoms} atoms")
        
        # Create Metropolis sampler
        print("\n2. Creating Metropolis sampler...")
        metropolis = create_metropolis_sampler()
        print("✓ Metropolis sampler created")
        
        # Run basic simulation
        print("\n3. Running basic Monte Carlo simulation...")
        acceptance_rate = run_monte_carlo_simulation(box, metropolis, n_steps=1000)
        
        if 0.2 <= acceptance_rate <= 0.8:
            print("✓ Acceptance rate is in reasonable range (20-80%)")
        else:
            print(f"⚠ Acceptance rate {acceptance_rate:.3f} is outside ideal range")
        
        # Test temperature effects
        test_temperature_effects()
        
        # Validate detailed balance
        validate_detailed_balance()
        
        print("\n" + "="*60)
        print("🎉 Metropolis sampling demonstration completed successfully!")
        print("="*60)
        
    except Exception as e:
        print(f"\n❌ Error during demonstration: {e}")
        import traceback
        traceback.print_exc()
        return False
    
    return True

if __name__ == '__main__':
    success = main()
    sys.exit(0 if success else 1) 