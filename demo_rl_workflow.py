#!/usr/bin/env python3
"""
Demonstration of the RL training and inference workflow.

This script shows:
1. How to run RL training to improve the neural network
2. How to run inference with the trained model
3. How the model selection works (RL model vs base model)
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

import subprocess
import time


def print_section(title):
    """Print a formatted section header."""
    print("\n" + "="*80)
    print(f"  {title}")
    print("="*80 + "\n")


def check_model_files():
    """Check which model files exist."""
    print_section("Checking Model Files")
    
    base_model = "model_01.pt"
    rl_model = "model_01_rl.pt"
    
    print(f"Base model (model_01.pt): {'EXISTS' if os.path.exists(base_model) else 'NOT FOUND'}")
    print(f"RL model (model_01_rl.pt): {'EXISTS' if os.path.exists(rl_model) else 'NOT FOUND'}")
    
    return os.path.exists(base_model), os.path.exists(rl_model)


def run_rl_training():
    """Run the RL training script."""
    print_section("Running RL Training")
    print("This will train the neural network to improve acceptance rates...")
    print("Training will run for multiple iterations (~2-5 minutes)\n")
    
    start_time = time.time()
    
    # Run the training script
    result = subprocess.run(
        ["python3", "train_nn_with_rl.py"],
        capture_output=False,
        text=True
    )
    
    elapsed_time = time.time() - start_time
    
    if result.returncode == 0:
        print(f"\n✓ RL training completed successfully in {elapsed_time:.1f} seconds")
        return True
    else:
        print(f"\n✗ RL training failed with return code {result.returncode}")
        return False


def run_inference():
    """Run inference with the trained model."""
    print_section("Running Inference with Trained Model")
    print("Running simulation with the RL-optimized neural network...\n")
    
    start_time = time.time()
    
    # Run the inference script
    result = subprocess.run(
        ["python3", "run_nn_translate_simulation.py"],
        capture_output=False,
        text=True
    )
    
    elapsed_time = time.time() - start_time
    
    if result.returncode == 0:
        print(f"\n✓ Inference completed successfully in {elapsed_time:.1f} seconds")
        return True
    else:
        print(f"\n✗ Inference failed with return code {result.returncode}")
        return False


def main():
    """Main demonstration workflow."""
    print("="*80)
    print("  REINFORCEMENT LEARNING WORKFLOW DEMONSTRATION")
    print("="*80)
    
    # Check initial state
    has_base, has_rl = check_model_files()
    
    if not has_base:
        print("\n⚠ Warning: Base model (model_01.pt) not found!")
        print("The system will start with random initialization.")
    
    # Ask user what to do
    print("\nOptions:")
    print("  1. Run RL training (creates/updates model_01_rl.pt)")
    print("  2. Run inference only (uses existing model)")
    print("  3. Full workflow (training + inference)")
    print("  4. Exit")
    
    try:
        choice = input("\nSelect option (1-4): ").strip()
    except (EOFError, KeyboardInterrupt):
        print("\nExiting...")
        return
    
    if choice == "1":
        success = run_rl_training()
        if success:
            check_model_files()
    
    elif choice == "2":
        if not has_base and not has_rl:
            print("\n⚠ No model files found. Please run training first or provide model_01.pt")
            return
        success = run_inference()
    
    elif choice == "3":
        # Full workflow
        print_section("Full Workflow: Training + Inference")
        
        # Step 1: Training
        success = run_rl_training()
        if not success:
            print("\nTraining failed, skipping inference")
            return
        
        # Check models
        check_model_files()
        
        # Step 2: Inference
        time.sleep(1)  # Brief pause
        success = run_inference()
    
    elif choice == "4":
        print("\nExiting...")
        return
    
    else:
        print(f"\nInvalid choice: {choice}")
        return
    
    print_section("Workflow Complete")
    print("Summary:")
    print(f"  Base model: {'EXISTS' if os.path.exists('model_01.pt') else 'NOT FOUND'}")
    print(f"  RL model: {'EXISTS' if os.path.exists('model_01_rl.pt') else 'NOT FOUND'}")
    print("\nThe trained RL model will be automatically used in future simulations.")


if __name__ == "__main__":
    main()

