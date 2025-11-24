#!/usr/bin/env python3
"""
Multi-GPU Parallel RL training script for NNTranslate move using PPO.

This script parallelizes training across multiple GPUs by:
1. Running independent MC simulations on each GPU simultaneously
2. Collecting experiences from all GPUs
3. Aggregating experiences and training on the combined dataset
4. Synchronizing model weights across GPUs

This provides near-linear speedup with the number of GPUs.
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
import torch.multiprocessing as mp
from collections import deque
import time
from queue import Queue

from src_python.Molecule_Definition import Molecule_Type
from src_python.Box_CubeBox import CubeBox
from src_python.FF_LJ_Cut import LJ_Cut
from src_python.MC_Move_NNTranslate import NNTranslate, SimpleFeedforwardNet2Feature
from src_python.Sampling_Metropolis import Metropolis
from src_python.Sim_MonteCarlo import SimMonteCarlo
from src_python.MC_Move_MolTranslation import MolTranslate
from math import sqrt

# Import base classes from single-GPU version
from train_nn_with_rl import (
    RLNNTranslate, 
    PPOTrainer, 
    create_fcc_with_vacancies,
    setup_lj_forcefield
)


def check_gpu_availability():
    """Check how many GPUs are available."""
    if not torch.cuda.is_available():
        print("No CUDA GPUs available. Falling back to CPU.")
        return 0
    
    n_gpus = torch.cuda.device_count()
    print(f"Found {n_gpus} CUDA GPU(s):")
    for i in range(n_gpus):
        props = torch.cuda.get_device_properties(i)
        print(f"  GPU {i}: {props.name} ({props.total_memory / 1e9:.2f} GB)")
    
    return n_gpus


def run_parallel_episode(gpu_id, model_state_dict, n_cycles, moves_per_cycle, 
                        temperature, return_queue, seed_offset=0):
    """
    Run a single training episode on a specific GPU.
    
    Args:
        gpu_id: GPU device ID to use
        model_state_dict: Model state dictionary to load
        n_cycles: Number of MC cycles
        moves_per_cycle: Moves per cycle
        temperature: System temperature
        return_queue: Queue to return results
        seed_offset: Random seed offset for this process
    """
    try:
        # Set random seeds (different for each GPU)
        np.random.seed(42 + seed_offset + gpu_id * 1000)
        torch.manual_seed(42 + seed_offset + gpu_id * 1000)
        
        # Set device for this process
        device = torch.device(f'cuda:{gpu_id}')
        
        # Create system
        box = create_fcc_with_vacancies()
        setup_lj_forcefield(box)
        box.temperature = temperature
        box.beta = 1.0 / temperature
        box.compute_energy()
        
        # Create and configure NN move
        nn_move = RLNNTranslate([box])
        
        # Load model on this GPU
        model = SimpleFeedforwardNet2Feature()
        model.load_state_dict(model_state_dict)
        model.to(device)
        model.eval()
        
        nn_move.nn_model = model
        nn_move.device = device
        
        # Run episode
        sampling = Metropolis()
        
        mc = SimMonteCarlo(
            nCycles=n_cycles,
            nMoves=moves_per_cycle,
            screenfreq=0,  # Suppress output in parallel mode
            configfreq=0,
            energyCheck=0
        )
        
        mc.BoxList = [box]
        
        # Traditional and NN moves (25:1 ratio)
        uniform_move = MolTranslate([box], limit=8.0, max_dist=0.1)
        uniform_move.maintFreq = 10
        
        mc.Moves = [uniform_move, nn_move]
        mc.moveweights = [25.0, 1.0]
        mc.Sampling = sampling
        
        # Enable experience collection
        nn_move.collect_experiences = True
        nn_move.clear_experiences()
        
        initial_energy = box.ETotal
        
        # Run simulation
        mc.run_monte_carlo()
        
        final_energy = box.ETotal
        
        # Collect results
        experiences = nn_move.get_experiences()
        acceptance_rate = nn_move.accpt / nn_move.atmps if nn_move.atmps > 0 else 0.0
        
        # Move experiences to CPU to avoid GPU memory issues
        experiences_cpu = []
        for exp in experiences:
            exp_cpu = {
                'features': exp['features'].cpu().clone(),
                'log_probs': exp['log_probs'].copy(),
                'bin_idx': exp['bin_idx'],
                'accepted': exp['accepted'],
                'energy_change': exp['energy_change'],
                'displacement': exp['displacement']
            }
            experiences_cpu.append(exp_cpu)
        
        result = {
            'gpu_id': gpu_id,
            'experiences': experiences_cpu,
            'acceptance_rate': acceptance_rate,
            'energy_change': final_energy - initial_energy,
            'n_experiences': len(experiences_cpu),
            'n_attempts': int(nn_move.atmps),
            'n_accepts': int(nn_move.accpt)
        }
        
        # Return result through queue
        return_queue.put(result)
        
    except Exception as e:
        print(f"Error on GPU {gpu_id}: {e}")
        import traceback
        traceback.print_exc()
        return_queue.put({
            'gpu_id': gpu_id,
            'error': str(e),
            'experiences': [],
            'acceptance_rate': 0.0,
            'energy_change': 0.0,
            'n_experiences': 0
        })


def run_parallel_training_iteration(n_gpus, model, n_cycles, moves_per_cycle, 
                                   temperature, iteration):
    """
    Run training iteration in parallel across multiple GPUs.
    
    Args:
        n_gpus: Number of GPUs to use
        model: Current model
        n_cycles: MC cycles per GPU
        moves_per_cycle: Moves per cycle
        temperature: System temperature
        iteration: Current iteration number (for seeding)
    
    Returns:
        all_experiences: Combined experiences from all GPUs
        stats: Dictionary of statistics
    """
    # Get model state dict (on CPU)
    model_state_dict = {k: v.cpu() for k, v in model.state_dict().items()}
    
    # Create process queue for results
    ctx = mp.get_context('spawn')
    return_queue = ctx.Queue()
    
    # Launch parallel processes
    processes = []
    for gpu_id in range(n_gpus):
        p = ctx.Process(
            target=run_parallel_episode,
            args=(gpu_id, model_state_dict, n_cycles, moves_per_cycle, 
                  temperature, return_queue, iteration)
        )
        p.start()
        processes.append(p)
    
    # Collect results
    results = []
    for _ in range(n_gpus):
        result = return_queue.get()
        results.append(result)
    
    # Wait for all processes to complete
    for p in processes:
        p.join()
    
    # Aggregate experiences from all GPUs
    all_experiences = []
    total_attempts = 0
    total_accepts = 0
    total_energy_change = 0.0
    
    for result in results:
        if 'error' in result:
            print(f"  ⚠ GPU {result['gpu_id']} encountered error: {result['error']}")
            continue
        
        all_experiences.extend(result['experiences'])
        total_attempts += result['n_attempts']
        total_accepts += result['n_accepts']
        total_energy_change += result['energy_change']
        
        print(f"  GPU {result['gpu_id']}: {result['n_experiences']} experiences, "
              f"acceptance rate: {result['acceptance_rate']:.4f}")
    
    # Calculate aggregate statistics
    aggregate_acceptance = total_accepts / total_attempts if total_attempts > 0 else 0.0
    
    stats = {
        'n_experiences': len(all_experiences),
        'acceptance_rate': aggregate_acceptance,
        'total_attempts': total_attempts,
        'total_accepts': total_accepts,
        'avg_energy_change': total_energy_change / n_gpus
    }
    
    return all_experiences, stats


def load_or_create_model(rl_model_path="model_01_rl.pt", base_model_path="model_01.pt"):
    """Load model, prioritizing RL-trained version."""
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


def main():
    """Main parallel RL training loop."""
    print("="*80)
    print("  MULTI-GPU PARALLEL RL TRAINING FOR NNTRANSLATE")
    print("="*80)
    
    # Check GPU availability
    n_gpus = check_gpu_availability()
    
    if n_gpus == 0:
        print("\nNo GPUs available. Please use train_nn_with_rl.py for CPU training.")
        return None
    
    # Ask user how many GPUs to use
    print(f"\nHow many GPUs would you like to use? (1-{n_gpus})")
    try:
        user_input = input(f"Enter number [default: {n_gpus}]: ").strip()
        if user_input == "":
            n_gpus_to_use = n_gpus
        else:
            n_gpus_to_use = int(user_input)
            if n_gpus_to_use < 1 or n_gpus_to_use > n_gpus:
                print(f"Invalid number. Using all {n_gpus} GPUs.")
                n_gpus_to_use = n_gpus
    except (ValueError, EOFError, KeyboardInterrupt):
        print(f"Using all {n_gpus} GPUs.")
        n_gpus_to_use = n_gpus
    
    # Training hyperparameters
    n_training_iterations = 20
    n_cycles_per_gpu = 50  # Each GPU runs this many cycles
    moves_per_cycle = 26
    save_interval = 5
    temperature = 0.7
    
    print(f"\n" + "="*80)
    print(f"Training Configuration:")
    print(f"  GPUs to use: {n_gpus_to_use}")
    print(f"  Training iterations: {n_training_iterations}")
    print(f"  MC cycles per GPU: {n_cycles_per_gpu}")
    print(f"  Total cycles per iteration: {n_cycles_per_gpu * n_gpus_to_use}")
    print(f"  Moves per cycle: {moves_per_cycle}")
    print(f"  Expected experiences per iteration: ~{n_cycles_per_gpu * n_gpus_to_use}")
    print(f"  Temperature: {temperature}")
    print(f"  Speedup: ~{n_gpus_to_use}× faster than single GPU")
    print("="*80)
    
    # Load or create model
    model, is_rl_model = load_or_create_model()
    model.train()
    
    # Initialize PPO trainer (on GPU 0 for training)
    device = torch.device('cuda:0')
    model.to(device)
    trainer = PPOTrainer(model, lr=1e-4)
    
    # Track best model
    best_acceptance_rate = 0.0
    best_model_state = None
    
    print(f"\nStarting parallel RL training...")
    print("="*80)
    
    start_time = time.time()
    
    for iteration in range(n_training_iterations):
        iter_start = time.time()
        
        print(f"\n{'='*80}")
        print(f"  Training Iteration {iteration + 1}/{n_training_iterations}")
        print(f"{'='*80}")
        
        # Run parallel experience collection
        print(f"Running {n_gpus_to_use} parallel simulations...")
        experiences, stats = run_parallel_training_iteration(
            n_gpus_to_use, model, n_cycles_per_gpu, moves_per_cycle,
            temperature, iteration
        )
        
        print(f"\nAggregated results:")
        print(f"  Total experiences: {stats['n_experiences']}")
        print(f"  Total attempts: {stats['total_attempts']}")
        print(f"  Total accepts: {stats['total_accepts']}")
        print(f"  Aggregate acceptance rate: {stats['acceptance_rate']:.4f} "
              f"({stats['acceptance_rate']*100:.2f}%)")
        print(f"  Average energy change: {stats['avg_energy_change']:.6f}")
        
        # Train on aggregated experiences
        if len(experiences) > 0:
            print(f"\nTraining on {len(experiences)} experiences...")
            
            # Move experiences back to GPU 0 for training
            for exp in experiences:
                exp['features'] = exp['features'].to(device)
            
            train_stats = trainer.train_on_experiences(experiences, n_epochs=4)
            
            print(f"Training stats:")
            print(f"  Loss: {train_stats['loss']:.6f}")
            print(f"  Policy loss: {train_stats['policy_loss']:.6f}")
            print(f"  Entropy: {train_stats['entropy']:.6f}")
            print(f"  Avg reward: {train_stats['avg_reward']:.4f}")
        
        # Track best model
        if stats['acceptance_rate'] > best_acceptance_rate:
            best_acceptance_rate = stats['acceptance_rate']
            best_model_state = model.state_dict().copy()
            print(f"\n  *** New best acceptance rate: {best_acceptance_rate:.4f} ***")
        
        # Save checkpoint
        if (iteration + 1) % save_interval == 0:
            if best_model_state is not None:
                torch.save(best_model_state, "model_01_rl.pt")
                print(f"  Saved best model to model_01_rl.pt")
        
        iter_time = time.time() - iter_start
        print(f"\nIteration time: {iter_time:.2f} seconds")
        print(f"Estimated remaining time: {iter_time * (n_training_iterations - iteration - 1) / 60:.1f} minutes")
    
    end_time = time.time()
    training_time = end_time - start_time
    
    # Final save
    if best_model_state is not None:
        torch.save(best_model_state, "model_01_rl.pt")
    
    print("\n" + "="*80)
    print("  TRAINING COMPLETE")
    print("="*80)
    print(f"Total training time: {training_time:.2f} seconds ({training_time/60:.1f} minutes)")
    print(f"Average time per iteration: {training_time/n_training_iterations:.2f} seconds")
    print(f"Best acceptance rate achieved: {best_acceptance_rate:.4f} ({best_acceptance_rate*100:.2f}%)")
    print(f"Final model saved to: model_01_rl.pt")
    print(f"\nSpeedup vs single GPU: ~{n_gpus_to_use}×")
    print(f"Total experiences collected: {n_training_iterations * n_cycles_per_gpu * n_gpus_to_use}")
    
    # Print training history
    if len(trainer.training_stats) > 0:
        print("\nTraining History:")
        print(f"{'Iter':<6} {'Loss':<10} {'Reward':<10} {'Accept Rate':<12}")
        print("-" * 40)
        for i, stats in enumerate(trainer.training_stats):
            print(f"{i+1:<6} {stats['loss']:<10.6f} {stats['avg_reward']:<10.4f} {stats['acceptance_rate']:<12.4f}")
    
    print("="*80)
    
    return {
        'best_acceptance_rate': best_acceptance_rate,
        'training_time': training_time,
        'n_gpus_used': n_gpus_to_use,
        'training_stats': trainer.training_stats
    }


if __name__ == "__main__":
    # Set multiprocessing start method
    mp.set_start_method('spawn', force=True)
    
    results = main()
    
    if results is not None:
        print(f"\nFinal Results:")
        print(f"  Best Acceptance Rate: {results['best_acceptance_rate']:.4f}")
        print(f"  Training Time: {results['training_time']:.2f} seconds")
        print(f"  GPUs Used: {results['n_gpus_used']}")
        print(f"  Speedup: ~{results['n_gpus_used']}×")

