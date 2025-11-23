#!/usr/bin/env python3
"""
RL Optimization for NNTranslate move.
Reuses system setup from run_nn_translate_simulation.py.
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import time
import torch
import torch.optim as optim

# Import from source
from src_python.MC_Move_NNTranslate import NNTranslate, SimpleFeedforwardNet2Feature
from src_python.CoordinateTypes import Displacement
from src_python.Box_SimpleBox import SimpleBox
from src_python.Sampling_Metropolis import Metropolis
from src_python.Sim_MonteCarlo import SimMonteCarlo
import run_nn_translate_simulation as base_sim

class RLNNTranslate(NNTranslate):
    def __init__(self, BoxArray, learning_rate=1e-4, batch_size=32):
        super().__init__(BoxArray)
        self.learning_rate = learning_rate
        self.optimizer = None
        self.batch_size = batch_size
        self.saved_log_probs = []
        self.rewards = []
        self.rl_moves_count = 0
        self.training = True
        self.running_reward = 0
        
    def full_move(self, trial_box: SimpleBox, sampling):
        # Override full_move to enable gradients and RL training
        
        box_idx = trial_box.boxID - 1
        self.atmps += 1.0
        self.boxatmps[box_idx] += 1.0

        target_mol = trial_box.pick_random_molecule()
        atom_indices = target_mol['atomindicies'][0]
        target_atom_idx = atom_indices[0]
        target_mol_idx = target_mol['mol_index']
        target_pos = trial_box.atoms[target_atom_idx]
        mol_type = target_mol['mol_type']
        mol_idx = target_mol_idx

        neighbor_positions = self._gather_neighbors(trial_box, target_pos, target_mol_idx)
        
        features = self._compute_2channel_features(trial_box, target_pos, neighbor_positions)
        
        # Enable gradient for policy update
        self.nn_model.train() 
        
        features_tensor = features.unsqueeze(0)
        # Compute log probs (with gradient)
        log_probs_tensor = self.nn_model(features_tensor).squeeze(0)
        
        # Detach for sampling/Sim usage (expects numpy)
        self.log_probs = log_probs_tensor.detach().cpu().numpy()
        
        # Sample bin
        selected_bin_idx, forward_prob = self._gumbel_sample(self.log_probs)
        
        # Store log prob of selected action for REINFORCE
        # selected_bin_idx is tuple (i, j, k)
        selected_log_prob = log_probs_tensor[selected_bin_idx]
        self.saved_log_probs.append(selected_log_prob)
        
        # Generate position
        delta, bin_delta = self._generate_position_in_bin(selected_bin_idx)
        new_pos_raw = target_pos + delta
        new_pos = trial_box.boundary(new_pos_raw)

        disp = Displacement(
            molType=mol_type,
            molIndx=mol_idx,
            atmIndicies=np.array([target_atom_idx]),
            newPositions=new_pos.reshape(1, 3)
        )

        accepted = False
        
        # Check constraints and energy
        if not trial_box.check_constraint(disp):
            self.constrainrej += 1
            accepted = False
        else:
            e_inter, e_intra, potential_accept = trial_box.compute_energy_delta(
                disp, self.tempList, self.tempNnei, computeintra=False
            )
            if not potential_accept:
                self.ovlaprej += 1
                accepted = False
            else:
                e_diff = e_inter + e_intra
                if not trial_box.check_post_energy(disp, e_diff):
                    self.constrainrej += 1
                    accepted = False
                else:
                    # Bias correction
                    reverse_prob = self._calculate_reverse_probability(trial_box, new_pos_raw, target_pos, bin_delta, target_mol_idx)
                    bias_correction = reverse_prob - forward_prob
                    
                    accept = sampling.make_decision(trial_box, e_diff, disp, log_prob=bias_correction)
                    if accept:
                        self.accpt += 1.0
                        self.boxaccpt[box_idx] += 1.0
                        trial_box.update_energy(e_diff, e_inter, e_intra)
                        trial_box.update_position(disp)
                        self.largest_displacement = max(self.largest_displacement, np.linalg.norm(delta))
                        accepted = True
                    else:
                        self.detailedrej += 1
                        accepted = False
        
        # RL Reward
        # +1 for accept, -1 for reject to maximize acceptance
        reward = 1.0 if accepted else -1.0
        
        self.rewards.append(reward)
        self.rl_moves_count += 1
        
        if len(self.saved_log_probs) >= self.batch_size:
            self.update_policy()
            
        return accepted

    def update_policy(self):
        policy_loss = []
        
        # Basic REINFORCE
        for log_prob, reward in zip(self.saved_log_probs, self.rewards):
            policy_loss.append(-log_prob * reward)
            
        self.optimizer.zero_grad()
        
        if len(policy_loss) > 0:
            loss = torch.stack(policy_loss).sum()
            loss.backward()
            self.optimizer.step()
            
            self.running_reward += sum(self.rewards)
            
            # Log periodically
            if (self.rl_moves_count // self.batch_size) % 10 == 0:
                avg_reward = sum(self.rewards) / len(self.rewards)
                print(f"RL Update (Move {self.rl_moves_count}): Loss {loss.item():.4f}, Avg Reward {avg_reward:.4f}")
        
        # Clear memory
        self.saved_log_probs = []
        self.rewards = []
        
        self.nn_model.eval()

def setup_rl_nntranslate_move(box):
    print("\nSetting up RL NNTranslate move...")
    # Initialize RL move
    nn_move = RLNNTranslate([box], learning_rate=1e-4, batch_size=64)
    
    # Load initial model if available, else random init
    if os.path.exists("model_01.pt"):
        print("Loading model_01.pt...")
        nn_move.create_and_set_neural_network("model_01.pt")
    else:
        print("Warning: model_01.pt not found, initializing random model.")
        nn_move.create_and_set_neural_network(None)
        
    # Re-initialize optimizer since model parameters might have changed/loaded?
    # create_and_set_neural_network creates a NEW model instance.
    # We need to update the optimizer to track the new parameters.
    nn_move.optimizer = optim.Adam(nn_move.nn_model.parameters(), lr=1e-4)
    
    return nn_move

def main():
    print("="*80)
    print("RL OPTIMIZATION FOR NNTRANSLATE")
    print("="*80)
    
    temperature = 0.7
    
    try:
        # Step 1: Create system (Reuse existing function)
        box = base_sim.create_fcc_with_vacancies()
        
        # Step 2: Setup Forcefield
        lj_ff = base_sim.setup_lj_forcefield(box)
        
        # Step 3: Setup Temperature
        base_sim.setup_temperature(box, temperature)
        
        # Step 4: Setup RL Move
        nn_move = setup_rl_nntranslate_move(box)
        
        # Step 5: Initial Energy
        base_sim.compute_initial_energy(box)
        
        # Step 6: Run Simulation
        # 100 cycles * 500 moves = 50,000 moves for training
        n_cycles = 100 
        moves_per_cycle = 500
        
        print(f"\nStarting Training: {n_cycles} cycles, {moves_per_cycle} moves/cycle")
        
        sampling = Metropolis()
        mc = SimMonteCarlo(
            nCycles=n_cycles, 
            nMoves=moves_per_cycle, 
            screenfreq=10, 
            configfreq=0, 
            energyCheck=1000
        )
        
        mc.BoxList = [box]
        
        # Use ONLY the RL move for training efficiency
        mc.Moves = [nn_move]
        mc.moveweights = [1.0]
        
        mc.Sampling = sampling
        
        start_time = time.time()
        mc.run_monte_carlo()
        end_time = time.time()
        
        print(f"\nTraining Time: {end_time - start_time:.2f} seconds")
        
        # Save the optimized model
        output_model = "model_rf_01.pt"
        print(f"\nSaving optimized model to {output_model}...")
        torch.save(nn_move.nn_model.state_dict(), output_model)
        print("Save complete.")
        
        # Stats
        print("\nFinal Statistics:")
        print(f"Total Attempts: {nn_move.atmps}")
        print(f"Total Accepts: {nn_move.accpt}")
        if nn_move.atmps > 0:
            print(f"Acceptance Rate: {nn_move.accpt / nn_move.atmps:.4f}")
        
    except Exception as e:
        print(f"\nERROR: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()

