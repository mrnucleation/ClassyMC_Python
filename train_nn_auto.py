#!/usr/bin/env python3
"""
Automatic RL training launcher.

This script automatically detects available GPUs and chooses the best training mode:
- Multi-GPU parallel training (if multiple GPUs available)
- Single-GPU training (if one GPU available)
- CPU training (if no GPUs available)
"""

import sys
import os
import torch

def main():
    """Auto-detect hardware and launch appropriate training script."""
    
    print("="*80)
    print("  AUTOMATIC RL TRAINING LAUNCHER")
    print("="*80)
    
    # Check GPU availability
    if not torch.cuda.is_available():
        print("\nNo CUDA GPUs detected.")
        print("Launching CPU training mode...")
        print("="*80 + "\n")
        
        # Import and run CPU training
        import train_nn_with_rl
        results = train_nn_with_rl.main()
        return results
    
    n_gpus = torch.cuda.device_count()
    
    print(f"\nDetected {n_gpus} CUDA GPU(s):")
    for i in range(n_gpus):
        props = torch.cuda.get_device_properties(i)
        mem_gb = props.total_memory / 1e9
        print(f"  GPU {i}: {props.name} ({mem_gb:.1f} GB)")
    
    # Decide on training mode
    if n_gpus == 1:
        print("\n1 GPU detected.")
        print("Launching single-GPU training mode...")
        print("="*80 + "\n")
        
        # Import and run single-GPU training
        import train_nn_with_rl
        results = train_nn_with_rl.main()
        return results
    
    else:  # Multiple GPUs
        print(f"\n{n_gpus} GPUs detected.")
        print("Multi-GPU parallel training is available!")
        print("\nOptions:")
        print("  1. Multi-GPU parallel training (RECOMMENDED)")
        print(f"     - Uses all {n_gpus} GPUs simultaneously")
        print(f"     - Expected speedup: ~{n_gpus}×")
        print("  2. Single-GPU training")
        print("     - Uses only GPU 0")
        print("  3. Exit")
        
        try:
            choice = input("\nSelect option (1-3) [default: 1]: ").strip()
            if choice == "" or choice == "1":
                print("\nLaunching multi-GPU parallel training...")
                print("="*80 + "\n")
                
                # Import and run parallel training
                import train_nn_with_rl_parallel
                results = train_nn_with_rl_parallel.main()
                return results
            
            elif choice == "2":
                print("\nLaunching single-GPU training...")
                print("="*80 + "\n")
                
                # Import and run single-GPU training
                import train_nn_with_rl
                results = train_nn_with_rl.main()
                return results
            
            elif choice == "3":
                print("\nExiting...")
                return None
            
            else:
                print(f"\nInvalid choice: {choice}")
                print("Defaulting to multi-GPU training...")
                print("="*80 + "\n")
                
                import train_nn_with_rl_parallel
                results = train_nn_with_rl_parallel.main()
                return results
                
        except (EOFError, KeyboardInterrupt):
            print("\n\nExiting...")
            return None


if __name__ == "__main__":
    results = main()
    
    if results is not None:
        print("\n" + "="*80)
        print("  TRAINING SUMMARY")
        print("="*80)
        print(f"Best Acceptance Rate: {results['best_acceptance_rate']:.4f}")
        print(f"Training Time: {results['training_time']:.2f} seconds")
        
        if 'n_gpus_used' in results:
            print(f"GPUs Used: {results['n_gpus_used']}")
        
        print("\nModel saved to: model_01_rl.pt")
        print("Future simulations will automatically use this model.")
        print("="*80)

