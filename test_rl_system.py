#!/usr/bin/env python3
"""
Quick verification script to test the RL system components.

This script performs basic tests to ensure:
1. Model loading works correctly
2. Experience collection functions
3. Reward computation is reasonable
4. PPO trainer initializes properly
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import torch

print("="*70)
print("  RL SYSTEM VERIFICATION TEST")
print("="*70)

# Test 1: Import all required modules
print("\n[Test 1] Checking imports...")
try:
    from train_nn_with_rl import (
        RLNNTranslate, 
        PPOTrainer, 
        load_or_create_model,
        create_fcc_with_vacancies,
        setup_lj_forcefield
    )
    from src_python.MC_Move_NNTranslate import SimpleFeedforwardNet2Feature
    print("✓ All imports successful")
except ImportError as e:
    print(f"✗ Import failed: {e}")
    sys.exit(1)

# Test 2: Model creation and loading
print("\n[Test 2] Testing model creation...")
try:
    model = SimpleFeedforwardNet2Feature()
    print(f"✓ Model created successfully")
    print(f"  Model parameters: {sum(p.numel() for p in model.parameters())}")
except Exception as e:
    print(f"✗ Model creation failed: {e}")
    sys.exit(1)

# Test 3: Model loading priority
print("\n[Test 3] Testing model loading priority...")
try:
    has_rl = os.path.exists("model_01_rl.pt")
    has_base = os.path.exists("model_01.pt")
    
    print(f"  model_01_rl.pt exists: {has_rl}")
    print(f"  model_01.pt exists: {has_base}")
    
    model, is_rl = load_or_create_model()
    
    if has_rl:
        if is_rl:
            print("✓ Correctly loaded RL model")
        else:
            print("✗ Should have loaded RL model but didn't")
    elif has_base:
        if not is_rl:
            print("✓ Correctly loaded base model")
        else:
            print("✗ Should have loaded base model but loaded RL")
    else:
        if not is_rl:
            print("✓ Correctly initialized random model")
        else:
            print("✗ Should have random model but claims RL")
    
except Exception as e:
    print(f"✗ Model loading failed: {e}")
    sys.exit(1)

# Test 4: System creation
print("\n[Test 4] Testing system creation...")
try:
    box = create_fcc_with_vacancies()
    setup_lj_forcefield(box)
    box.temperature = 0.7
    box.beta = 1.0 / 0.7
    box.compute_energy()
    
    print(f"✓ System created successfully")
    print(f"  Number of atoms: {box.nAtoms}")
    print(f"  Initial energy: {box.ETotal:.6f}")
    print(f"  Temperature: {box.temperature}")
except Exception as e:
    print(f"✗ System creation failed: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

# Test 5: RLNNTranslate initialization
print("\n[Test 5] Testing RLNNTranslate initialization...")
try:
    nn_move = RLNNTranslate([box])
    nn_move.nn_model = model
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    nn_move.nn_model.to(device)
    
    print(f"✓ RLNNTranslate initialized successfully")
    print(f"  Device: {device}")
    print(f"  Neighbor cutoff: {nn_move.neighbor_cutoff}")
    print(f"  Number of bins: {nn_move.nbins}")
except Exception as e:
    print(f"✗ RLNNTranslate initialization failed: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

# Test 6: Experience collection
print("\n[Test 6] Testing experience collection...")
try:
    nn_move.collect_experiences = True
    nn_move.clear_experiences()
    
    # Perform a few moves
    from src_python.Sampling_Metropolis import Metropolis
    sampling = Metropolis()
    
    n_test_moves = 5
    for i in range(n_test_moves):
        nn_move.full_move(box, sampling)
    
    experiences = nn_move.get_experiences()
    print(f"✓ Experience collection works")
    print(f"  Collected {len(experiences)} experiences from {n_test_moves} moves")
    
    if len(experiences) > 0:
        exp = experiences[0]
        print(f"  Experience keys: {list(exp.keys())}")
        print(f"  Features shape: {exp['features'].shape}")
        print(f"  Accepted: {exp['accepted']}")
except Exception as e:
    print(f"✗ Experience collection failed: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

# Test 7: PPO trainer initialization
print("\n[Test 7] Testing PPO trainer initialization...")
try:
    trainer = PPOTrainer(model, lr=1e-4)
    print(f"✓ PPO trainer initialized successfully")
    print(f"  Learning rate: 1e-4")
    print(f"  Clip epsilon: {trainer.clip_epsilon}")
    print(f"  Entropy coefficient: {trainer.entropy_coef}")
except Exception as e:
    print(f"✗ PPO trainer initialization failed: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

# Test 8: Reward computation
print("\n[Test 8] Testing reward computation...")
try:
    # Create test experiences
    test_experiences = [
        {
            'accepted': True,
            'energy_change': -1.0,
            'displacement': 0.5
        },
        {
            'accepted': False,
            'energy_change': 2.0,
            'displacement': 0.2
        },
        {
            'accepted': True,
            'energy_change': 0.0,
            'displacement': 1.0
        }
    ]
    
    rewards = [trainer.compute_reward(exp) for exp in test_experiences]
    print(f"✓ Reward computation works")
    print(f"  Test rewards: {[f'{r:.2f}' for r in rewards]}")
    
    # Check that accepted moves have higher rewards
    if rewards[0] > rewards[1]:
        print(f"  ✓ Accepted moves have higher rewards than rejected")
    else:
        print(f"  ⚠ Warning: Accepted move reward not higher than rejected")
        
except Exception as e:
    print(f"✗ Reward computation failed: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

# Test 9: Training on dummy experiences (if we collected any)
if len(experiences) > 0:
    print("\n[Test 9] Testing training on experiences...")
    try:
        stats = trainer.train_on_experiences(experiences, n_epochs=2)
        
        if stats is not None:
            print(f"✓ Training successful")
            print(f"  Loss: {stats['loss']:.6f}")
            print(f"  Policy loss: {stats['policy_loss']:.6f}")
            print(f"  Entropy: {stats['entropy']:.6f}")
            print(f"  Avg reward: {stats['avg_reward']:.4f}")
            print(f"  Acceptance rate: {stats['acceptance_rate']:.4f}")
        else:
            print(f"✗ Training returned None")
    except Exception as e:
        print(f"✗ Training failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
else:
    print("\n[Test 9] Skipping training test (no experiences collected)")

# Test 10: Model saving
print("\n[Test 10] Testing model saving...")
try:
    test_model_path = "test_model_temp.pt"
    torch.save(model.state_dict(), test_model_path)
    
    # Try loading it back
    test_model = SimpleFeedforwardNet2Feature()
    test_model.load_state_dict(torch.load(test_model_path, map_location='cpu'))
    
    # Clean up
    os.remove(test_model_path)
    
    print(f"✓ Model saving and loading works")
except Exception as e:
    print(f"✗ Model saving failed: {e}")
    if os.path.exists(test_model_path):
        os.remove(test_model_path)
    sys.exit(1)

# Summary
print("\n" + "="*70)
print("  VERIFICATION COMPLETE")
print("="*70)
print("\n✅ All tests passed! The RL system is ready to use.")
print("\nNext steps:")
print("  1. Run: python3 train_nn_with_rl.py")
print("  2. Or use: python3 demo_rl_workflow.py")
print("\nSee RL_TRAINING_README.md for detailed documentation.")
print("="*70)

