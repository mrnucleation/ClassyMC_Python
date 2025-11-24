#!/usr/bin/env python3
"""
Reinforcement Learning training script for NNTranslate move using PPO.

This script trains the neural network to improve acceptance rates by:
1. Running Monte Carlo simulations with the current policy
2. Collecting experiences (features, actions, rewards)
3. Optimizing the neural network using PPO
4. Saving improved models periodically

The reward is based on:
- Move acceptance (primary reward)
- Energy improvement (secondary reward)
- Move diversity (exploration bonus)
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from collections import deque
import time

from src_python.Molecule_Definition import Molecule_Type
from src_python.Box_CubeBox import CubeBox
from src_python.FF_LJ_Cut import LJ_Cut
from src_python.MC_Move_NNTranslate import NNTranslate, SimpleFeedforwardNet2Feature
from src_python.Sampling_Metropolis import Metropolis
from src_python.Sim_MonteCarlo import SimMonteCarlo
from src_python.MC_Move_MolTranslation import MolTranslate
from math import sqrt

# Set random seeds
np.random.seed(42)
torch.manual_seed(42)


class RLNNTranslate(NNTranslate):
    """
    Extended NNTranslate that collects experiences for RL training.
    """
    def __init__(self, BoxArray):
        super().__init__(BoxArray)
        self.experiences = []
        self.current_features = None
        self.current_log_probs = None
        self.current_bin_idx = None
        self.collect_experiences = False
        
    def full_move(self, trial_box, sampling):
        """Override to collect experiences during moves."""
        box_idx = trial_box.boxID - 1
        self.atmps += 1.0
        self.boxatmps[box_idx] += 1.0
        
        # Select random atom
        target_mol = trial_box.pick_random_molecule()
        atom_indices = target_mol['atomindicies'][0]
        target_atom_idx = atom_indices[0]
        target_mol_idx = target_mol['mol_index']
        target_pos = trial_box.atoms[target_atom_idx]
        mol_idx = target_mol_idx
        mol_type = target_mol['mol_type']
        
        # Store old energy
        old_energy = trial_box.ETotal
        
        # Gather neighbors and compute features
        neighbor_positions = self._gather_neighbors(trial_box, target_pos, target_mol_idx)
        features = self._compute_2channel_features(trial_box, target_pos, neighbor_positions)
        
        with torch.no_grad():
            features_tensor = features.unsqueeze(0)
            log_probs = self.nn_model(features_tensor).squeeze(0)
            log_probs_np = log_probs.cpu().numpy()
        
        # Sample bin
        selected_bin_idx, forward_prob = self._gumbel_sample(log_probs_np)
        
        # Generate position
        delta, bin_delta = self._generate_position_in_bin(selected_bin_idx)
        new_pos_raw = target_pos + delta
        new_pos = trial_box.boundary(new_pos_raw)
        
        # Create displacement
        from src_python.CoordinateTypes import Displacement
        disp = Displacement(
            molType=mol_type,
            molIndx=mol_idx,
            atmIndicies=np.array([target_atom_idx]),
            newPositions=new_pos.reshape(1, 3)
        )
        
        # Check constraints
        if not trial_box.check_constraint(disp):
            self.constrainrej += 1
            if self.collect_experiences:
                self._store_experience(features, log_probs_np, selected_bin_idx, 
                                     accepted=False, energy_change=0.0, 
                                     displacement_magnitude=np.linalg.norm(delta))
            return False
        
        # Compute energy delta
        e_inter, e_intra, accept_energy = trial_box.compute_energy_delta(
            disp, self.tempList, self.tempNnei, computeintra=False
        )
        if not accept_energy:
            self.ovlaprej += 1
            if self.collect_experiences:
                self._store_experience(features, log_probs_np, selected_bin_idx,
                                     accepted=False, energy_change=0.0,
                                     displacement_magnitude=np.linalg.norm(delta))
            return False
        
        e_diff = e_inter + e_intra
        
        if not trial_box.check_post_energy(disp, e_diff):
            self.constrainrej += 1
            if self.collect_experiences:
                self._store_experience(features, log_probs_np, selected_bin_idx,
                                     accepted=False, energy_change=e_diff,
                                     displacement_magnitude=np.linalg.norm(delta))
            return False
        
        # Calculate bias correction
        reverse_prob = self._calculate_reverse_probability(
            trial_box, new_pos_raw, target_pos, bin_delta, target_mol_idx
        )
        bias_correction = reverse_prob - forward_prob
        
        # Metropolis acceptance
        accept = sampling.make_decision(trial_box, e_diff, disp, log_prob=bias_correction)
        
        if accept:
            self.accpt += 1.0
            self.boxaccpt[box_idx] += 1.0
            trial_box.update_energy(e_diff, e_inter, e_intra)
            trial_box.update_position(disp)
            self.largest_displacement = max(self.largest_displacement, np.linalg.norm(delta))
        else:
            self.detailedrej += 1
        
        # Store experience
        if self.collect_experiences:
            self._store_experience(features, log_probs_np, selected_bin_idx,
                                 accepted=accept, energy_change=e_diff,
                                 displacement_magnitude=np.linalg.norm(delta))
        
        return accept
    
    def _store_experience(self, features, log_probs, bin_idx, accepted, 
                         energy_change, displacement_magnitude):
        """Store experience for later training."""
        experience = {
            'features': features.cpu().clone(),  # Store on CPU to save GPU memory
            'log_probs': log_probs.copy(),
            'bin_idx': bin_idx,
            'accepted': accepted,
            'energy_change': energy_change,
            'displacement': displacement_magnitude
        }
        self.experiences.append(experience)
    
    def clear_experiences(self):
        """Clear collected experiences."""
        self.experiences = []
    
    def get_experiences(self):
        """Return collected experiences."""
        return self.experiences


class PPOTrainer:
    """
    PPO trainer for the neural network translation move.
    """
    def __init__(self, model, lr=3e-4, gamma=0.99, clip_epsilon=0.2, 
                 value_coef=0.5, entropy_coef=0.01):
        self.model = model
        self.device = next(model.parameters()).device
        self.optimizer = optim.Adam(model.parameters(), lr=lr)
        
        # PPO hyperparameters
        self.gamma = gamma
        self.clip_epsilon = clip_epsilon
        self.value_coef = value_coef
        self.entropy_coef = entropy_coef
        
        # Value network (separate head for advantage estimation)
        self.value_net = nn.Sequential(
            nn.Linear(2, 32),
            nn.ReLU(),
            nn.Linear(32, 16),
            nn.ReLU(),
            nn.Linear(16, 1)
        ).to(self.device)
        self.value_optimizer = optim.Adam(self.value_net.parameters(), lr=lr)
        
        # Statistics tracking
        self.training_stats = []
        
    def compute_reward(self, experience):
        """
        Compute reward for a single experience.
        
        Reward components:
        1. Acceptance: +1.0 if accepted, -0.1 if rejected
        2. Energy improvement: +0.5 if energy decreases
        3. Displacement diversity: +0.2 if displacement > 0.1 sigma
        """
        reward = 0.0
        
        # Primary: Acceptance
        if experience['accepted']:
            reward += 1.0
        else:
            reward -= 0.1
        
        # Secondary: Energy improvement
        if experience['energy_change'] < 0:
            reward += 0.5 * min(abs(experience['energy_change']) / 10.0, 1.0)
        
        # Tertiary: Encourage exploration (reasonable displacement)
        if 0.1 < experience['displacement'] < 1.5:
            reward += 0.2
        
        return reward
    
    def compute_advantages(self, rewards, features_batch):
        """Compute advantages using GAE (Generalized Advantage Estimation)."""
        # Simplified: Use rewards directly as advantages for now
        # In full PPO, would compute value estimates and TD errors
        advantages = torch.tensor(rewards, dtype=torch.float32, device=self.device)
        
        # Normalize advantages
        advantages = (advantages - advantages.mean()) / (advantages.std() + 1e-8)
        
        return advantages
    
    def train_on_experiences(self, experiences, n_epochs=4):
        """
        Train the model using collected experiences with PPO.
        
        Args:
            experiences: List of experience dictionaries
            n_epochs: Number of training epochs per batch
        """
        if len(experiences) == 0:
            return None
        
        # Prepare batch data
        features_list = []
        old_log_probs_list = []
        actions_list = []
        rewards_list = []
        
        for exp in experiences:
            features_list.append(exp['features'])
            
            # Get old log prob for the action taken
            bin_idx = exp['bin_idx']
            old_log_prob = exp['log_probs'][bin_idx]
            old_log_probs_list.append(old_log_prob)
            
            # Action is the bin index (flattened)
            nbins = exp['log_probs'].shape
            action = bin_idx[0] * nbins[1] * nbins[2] + bin_idx[1] * nbins[2] + bin_idx[2]
            actions_list.append(action)
            
            # Compute reward
            reward = self.compute_reward(exp)
            rewards_list.append(reward)
        
        # Convert to tensors
        features_batch = torch.stack(features_list).to(self.device)
        old_log_probs = torch.tensor(old_log_probs_list, dtype=torch.float32, device=self.device)
        actions = torch.tensor(actions_list, dtype=torch.long, device=self.device)
        
        # Compute advantages
        advantages = self.compute_advantages(rewards_list, features_batch)
        
        # PPO training loop
        total_loss = 0.0
        total_policy_loss = 0.0
        total_entropy = 0.0
        
        for epoch in range(n_epochs):
            # Forward pass
            log_probs_3d = self.model(features_batch)  # (batch, nbins[0], nbins[1], nbins[2])
            
            # Flatten for action selection
            batch_size = log_probs_3d.shape[0]
            log_probs_flat = log_probs_3d.view(batch_size, -1)  # (batch, nbins_total)
            
            # Get log probs for taken actions
            new_log_probs = log_probs_flat.gather(1, actions.unsqueeze(1)).squeeze(1)
            
            # Compute ratio for PPO
            ratio = torch.exp(new_log_probs - old_log_probs)
            
            # Clipped surrogate objective
            surr1 = ratio * advantages
            surr2 = torch.clamp(ratio, 1.0 - self.clip_epsilon, 1.0 + self.clip_epsilon) * advantages
            policy_loss = -torch.min(surr1, surr2).mean()
            
            # Entropy bonus (encourage exploration)
            entropy = -(torch.exp(log_probs_flat) * log_probs_flat).sum(dim=1).mean()
            
            # Total loss
            loss = policy_loss - self.entropy_coef * entropy
            
            # Optimization step
            self.optimizer.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(self.model.parameters(), max_norm=0.5)
            self.optimizer.step()
            
            total_loss += loss.item()
            total_policy_loss += policy_loss.item()
            total_entropy += entropy.item()
        
        # Average over epochs
        avg_loss = total_loss / n_epochs
        avg_policy_loss = total_policy_loss / n_epochs
        avg_entropy = total_entropy / n_epochs
        avg_reward = np.mean(rewards_list)
        
        stats = {
            'loss': avg_loss,
            'policy_loss': avg_policy_loss,
            'entropy': avg_entropy,
            'avg_reward': avg_reward,
            'acceptance_rate': sum(exp['accepted'] for exp in experiences) / len(experiences)
        }
        
        self.training_stats.append(stats)
        return stats


def create_fcc_with_vacancies(natoms=10, nvacancies=5):
    """Create FCC lattice with vacancies for testing."""
    print("Creating FCC lattice with vacancies...")
    
    lattice_constant = sqrt(8)*2**(1/6)/2.0
    n_cells = 3
    
    mock_lj_data = {"atoms": [("Ar", "LJ")]}
    atomtypes = ["LJ"]
    LJ_type = Molecule_Type(mock_lj_data, atomtypes=atomtypes)
    
    natoms_full = 4 * (n_cells ** 3)
    natoms_with_vacancies = natoms_full - nvacancies
    
    NMol = [natoms_with_vacancies]
    NMolMin = [0]
    NMolMax = [natoms_with_vacancies]
    
    box = CubeBox([LJ_type], NMolMin=NMolMin, NMolMax=NMolMax, NMol=NMol)
    box_length = n_cells * lattice_constant
    
    box.load_dimension([box_length])
    box.boxID = 1
    box.nDimension = 3
    box.NMol = np.array(NMol)
    box.NMolMin = np.array(NMolMin)
    box.NMolMax = np.array(NMolMax)
    box.nAtoms = natoms_with_vacancies
    box.nMaxAtoms = natoms_with_vacancies
    box.maxMol = len(NMol)
    box.nMolTotal = natoms_with_vacancies
    
    box.atoms = np.zeros((box.nMaxAtoms, 3), dtype=np.float64)
    box.AtomType = np.zeros(box.nMaxAtoms, dtype=int)
    box.MolType = np.zeros(box.nMaxAtoms, dtype=int)
    box.MolIndx = np.zeros(box.nMaxAtoms, dtype=int)
    box.MolSubIndx = np.zeros(box.nMaxAtoms, dtype=int)
    box.AtomSubIndx = np.zeros(box.nMaxAtoms, dtype=int)
    
    # Generate FCC positions
    fcc_basis = np.array([
        [0.0, 0.0, 0.0],
        [0.5, 0.5, 0.0],
        [0.5, 0.0, 0.5],
        [0.0, 0.5, 0.5]
    ]) * lattice_constant
    
    positions = []
    for i in range(n_cells):
        for j in range(n_cells):
            for k in range(n_cells):
                cell_origin = np.array([i, j, k]) * lattice_constant
                for basis_pos in fcc_basis:
                    positions.append(cell_origin + basis_pos)
    
    positions = np.array(positions)
    positions -= box_length / 2.0
    
    # Remove atoms to create vacancies
    center = np.array([0.0, 0.0, 0.0])
    distances = np.linalg.norm(positions - center, axis=1)
    vacancy_indices = np.argsort(distances)[:nvacancies]
    mask = np.ones(len(positions), dtype=bool)
    mask[vacancy_indices] = False
    positions = positions[mask]
    
    box.atoms[:natoms_with_vacancies] = positions
    box.AtomType[:natoms_with_vacancies] = 0
    box.MolType[:natoms_with_vacancies] = 0
    box.MolIndx[:natoms_with_vacancies] = np.arange(natoms_with_vacancies)
    box.MolSubIndx[:natoms_with_vacancies] = np.arange(natoms_with_vacancies)
    box.AtomSubIndx[:natoms_with_vacancies] = np.zeros(natoms_with_vacancies, dtype=int)
    
    print(f"Created FCC lattice: {natoms_with_vacancies} atoms with {nvacancies} vacancies")
    return box


def setup_lj_forcefield(box):
    """Set up LJ forcefield."""
    lj_ff = LJ_Cut(nAtomTypes=1)
    lj_ff.epsilon[0] = 1.0
    lj_ff.sigma[0] = 1.0
    lj_ff.rMin[0] = 0.4
    cutoff = 2.5
    lj_ff.rCut = cutoff
    lj_ff.rCutSq = cutoff**2
    lj_ff.epsTable[0, 0] = 4.0 * lj_ff.epsilon[0]
    lj_ff.sigTable[0, 0] = lj_ff.sigma[0]**2
    lj_ff.rMinTable[0, 0] = lj_ff.rMin[0]**2
    box.EFunc.append(lj_ff)
    return lj_ff


def load_or_create_model(rl_model_path="model_01_rl.pt", base_model_path="model_01.pt"):
    """
    Load model, prioritizing RL-trained version.
    
    Args:
        rl_model_path: Path to RL-trained model
        base_model_path: Path to base model
    
    Returns:
        model: Loaded PyTorch model
        is_rl_model: Boolean indicating if RL model was loaded
    """
    model = SimpleFeedforwardNet2Feature()
    
    if os.path.exists(rl_model_path):
        print(f"Loading RL-trained model from {rl_model_path}")
        model.load_state_dict(torch.load(rl_model_path, map_location='cpu'))
        return model, True
    elif os.path.exists(base_model_path):
        print(f"Loading base model from {base_model_path}")
        model.load_state_dict(torch.load(base_model_path, map_location='cpu'))
        return model, False
    else:
        print(f"No existing model found, starting with random initialization")
        return model, False


def run_rl_training_episode(box, nn_move, n_cycles=50, moves_per_cycle=26):
    """
    Run a training episode to collect experiences.
    
    Uses 25:1 ratio of traditional to NN moves for efficiency.
    """
    sampling = Metropolis()
    
    mc = SimMonteCarlo(
        nCycles=n_cycles,
        nMoves=moves_per_cycle,
        screenfreq=10,
        configfreq=0,
        energyCheck=0  # Disable energy checks during training
    )
    
    mc.BoxList = [box]
    
    # Traditional translation move (25 weight)
    uniform_move = MolTranslate([box], limit=8.0, max_dist=0.1)
    uniform_move.maintFreq = 10
    
    # NN move (1 weight) - 25:1 ratio
    mc.Moves = [uniform_move, nn_move]
    mc.moveweights = [25.0, 1.0]
    mc.Sampling = sampling
    
    # Enable experience collection
    nn_move.collect_experiences = True
    nn_move.clear_experiences()
    
    # Run simulation
    initial_energy = box.ETotal
    mc.run_monte_carlo()
    final_energy = box.ETotal
    
    # Collect statistics
    experiences = nn_move.get_experiences()
    acceptance_rate = nn_move.accpt / nn_move.atmps if nn_move.atmps > 0 else 0.0
    
    return experiences, acceptance_rate, final_energy - initial_energy


def main():
    """Main RL training loop."""
    print("="*80)
    print("REINFORCEMENT LEARNING TRAINING FOR NNTRANSLATE")
    print("="*80)
    
    # Training hyperparameters
    n_training_iterations = 20  # Number of training iterations
    n_cycles_per_iteration = 50  # MC cycles per iteration
    moves_per_cycle = 26  # 25 traditional + 1 NN (on average)
    
    save_interval = 5  # Save model every N iterations
    temperature = 0.7
    
    # Load or create model
    model, is_rl_model = load_or_create_model()
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model.to(device)
    model.train()  # Set to training mode
    
    # Initialize PPO trainer
    trainer = PPOTrainer(model, lr=1e-4)
    
    # Track best model
    best_acceptance_rate = 0.0
    best_model_state = None
    
    print(f"\nStarting RL training:")
    print(f"  Training iterations: {n_training_iterations}")
    print(f"  MC cycles per iteration: {n_cycles_per_iteration}")
    print(f"  Moves per cycle: {moves_per_cycle}")
    print(f"  Temperature: {temperature}")
    print(f"  Device: {device}")
    print("="*80)
    
    start_time = time.time()
    
    for iteration in range(n_training_iterations):
        print(f"\n--- Training Iteration {iteration + 1}/{n_training_iterations} ---")
        
        # Create fresh system for this iteration
        box = create_fcc_with_vacancies()
        setup_lj_forcefield(box)
        box.temperature = temperature
        box.beta = 1.0 / temperature
        
        # Compute initial energy
        box.compute_energy()
        print(f"Initial energy: {box.ETotal:.6f}")
        
        # Create RL-enabled NN move
        nn_move = RLNNTranslate([box])
        nn_move.nn_model = model
        nn_move.nn_model.to(device)
        
        # Run training episode
        experiences, acceptance_rate, energy_change = run_rl_training_episode(
            box, nn_move, n_cycles=n_cycles_per_iteration, moves_per_cycle=moves_per_cycle
        )
        
        print(f"Collected {len(experiences)} experiences")
        print(f"Acceptance rate: {acceptance_rate:.4f} ({acceptance_rate*100:.2f}%)")
        print(f"Energy change: {energy_change:.6f}")
        
        # Train on experiences
        if len(experiences) > 0:
            stats = trainer.train_on_experiences(experiences, n_epochs=4)
            print(f"Training stats:")
            print(f"  Loss: {stats['loss']:.6f}")
            print(f"  Policy loss: {stats['policy_loss']:.6f}")
            print(f"  Entropy: {stats['entropy']:.6f}")
            print(f"  Avg reward: {stats['avg_reward']:.4f}")
        
        # Track best model
        if acceptance_rate > best_acceptance_rate:
            best_acceptance_rate = acceptance_rate
            best_model_state = model.state_dict().copy()
            print(f"  *** New best acceptance rate: {best_acceptance_rate:.4f} ***")
        
        # Save checkpoint
        if (iteration + 1) % save_interval == 0:
            if best_model_state is not None:
                torch.save(best_model_state, "model_01_rl.pt")
                print(f"  Saved best model to model_01_rl.pt")
    
    end_time = time.time()
    training_time = end_time - start_time
    
    # Final save
    if best_model_state is not None:
        torch.save(best_model_state, "model_01_rl.pt")
    
    print("\n" + "="*80)
    print("TRAINING COMPLETE")
    print("="*80)
    print(f"Total training time: {training_time:.2f} seconds")
    print(f"Best acceptance rate achieved: {best_acceptance_rate:.4f} ({best_acceptance_rate*100:.2f}%)")
    print(f"Final model saved to: model_01_rl.pt")
    
    # Print training history
    if len(trainer.training_stats) > 0:
        print("\nTraining History:")
        print(f"{'Iter':<6} {'Loss':<10} {'Reward':<10} {'Accept Rate':<12}")
        print("-" * 40)
        for i, stats in enumerate(trainer.training_stats):
            print(f"{i+1:<6} {stats['loss']:<10.6f} {stats['avg_reward']:<10.4f} {stats['acceptance_rate']:<12.4f}")
    
    return {
        'best_acceptance_rate': best_acceptance_rate,
        'training_time': training_time,
        'training_stats': trainer.training_stats
    }


if __name__ == "__main__":
    results = main()
    print(f"\nFinal Results:")
    print(f"  Best Acceptance Rate: {results['best_acceptance_rate']:.4f}")
    print(f"  Training Time: {results['training_time']:.2f} seconds")

