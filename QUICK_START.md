# Quick Start Guide - RL Training System

## 🚀 Fastest Start (Recommended)

```bash
# 1. Run the verification test (optional but recommended)
python3 test_rl_system.py

# 2. Use the interactive demo
python3 demo_rl_workflow.py
```

Follow the on-screen menu to train and run simulations.

## 📝 Manual Workflow

### Step 1: Train the Model (First Time)

```bash
python3 train_nn_with_rl.py
```

**What happens:**
- Trains neural network using PPO for ~20 iterations
- Takes 2-5 minutes on CPU, <1 min on GPU
- Creates `model_01_rl.pt` (RL-optimized model)
- Saves best model automatically

**Expected output:**
```
Best acceptance rate achieved: 0.4286 (42.86%)
Final model saved to: model_01_rl.pt
```

### Step 2: Run Simulations (Always)

```bash
python3 run_nn_translate_simulation.py
```

**What happens:**
- Automatically uses `model_01_rl.pt` if it exists
- Falls back to `model_01.pt` if no RL model
- Runs Monte Carlo simulation with optimized NN moves
- Reports acceptance rates and statistics

## 🎯 Key Points

### Model Priority (Automatic)
1. `model_01_rl.pt` ← **Used first** (RL-trained)
2. `model_01.pt` ← Fallback (base model)
3. Random initialization ← Last resort

### The 25:1 Ratio
- **25 traditional moves** : **1 NN move**
- Why? Traditional moves are fast, NN moves are slower
- Still collects ~1000 NN experiences per training run
- Optimal balance of efficiency and data collection

### Reward Function
The NN learns to maximize:
- ✅ **Acceptance** (+1.0 if accepted, -0.1 if rejected)
- ✅ **Energy improvement** (up to +0.5 for lower energy)
- ✅ **Move diversity** (+0.2 for reasonable displacements)

## 🔧 Common Commands

### Train with default settings
```bash
python3 train_nn_with_rl.py
```

### Verify system works
```bash
python3 test_rl_system.py
```

### Run simulation
```bash
python3 run_nn_translate_simulation.py
```

### Check which models exist
```bash
ls -lh model_01*.pt
```

### Visualize training progress (future feature)
```bash
python3 plot_rl_progress.py --show
```

## 📊 What to Expect

### Before Training
```
Acceptance Rate: 10-30%
Energy sampling: Moderate
Exploration: Random-ish
```

### After Training
```
Acceptance Rate: 30-50%+ (typical 2× improvement)
Energy sampling: Better
Exploration: More physical
```

### Training Output Example
```
--- Training Iteration 5/20 ---
Collected 42 experiences
Acceptance rate: 0.4286 (42.86%)
Energy change: -2.453821
Training stats:
  Loss: 0.234567
  Avg reward: 0.8234
  *** New best acceptance rate: 0.4286 ***
  Saved best model to model_01_rl.pt
```

## ⚙️ Customization

### Quick Tweaks (edit `train_nn_with_rl.py`)

**Faster training (less accurate):**
```python
n_training_iterations = 10      # Less iterations
n_cycles_per_iteration = 25     # Fewer cycles
```

**Better training (slower):**
```python
n_training_iterations = 50      # More iterations
n_cycles_per_iteration = 100    # More cycles
```

**More efficient (fewer NN moves):**
```python
mc.moveweights = [50.0, 1.0]    # 50:1 ratio
```

**More NN training (less efficient):**
```python
mc.moveweights = [10.0, 1.0]    # 10:1 ratio
```

### Change Reward Emphasis

In `PPOTrainer.compute_reward()`:
```python
# Emphasize acceptance more
if accepted:
    reward += 2.0  # Was 1.0

# Care less about rejection
else:
    reward -= 0.05  # Was -0.1
```

## 🐛 Troubleshooting

### "No model file found"
**Solution:** Normal on first run. System will use random initialization.

### Low acceptance rate after training (<20%)
**Solution:** 
```python
# Increase training
n_training_iterations = 50

# Check base model exists
ls model_01.pt
```

### Training is slow
**Solution:**
```python
# Reduce iterations
n_training_iterations = 10

# Use more traditional moves
mc.moveweights = [50.0, 1.0]
```

### "CUDA out of memory"
**Solution:**
```python
# Force CPU usage
device = torch.device("cpu")

# Or reduce batch size (not directly configurable in current version)
```

## 📚 Documentation

- **This file**: Quick start
- **`RL_TRAINING_README.md`**: Detailed guide
- **`IMPLEMENTATION_SUMMARY.md`**: Technical details

## ✨ Tips

1. **Run verification first**: `python3 test_rl_system.py`
2. **Use demo for learning**: `python3 demo_rl_workflow.py`
3. **Train once, use forever**: The RL model persists
4. **Retrain periodically**: As you change systems/parameters
5. **Check model files**: `ls -lh model_01*.pt`

## 🎓 Understanding Output

### During Training
```
Acceptance rate: 0.4286 (42.86%)  ← Higher is better (aim for 40%+)
Loss: 0.234567                    ← Should decrease over time
Avg reward: 0.8234                ← Should increase over time
```

### During Simulation
```
Neural Network Translation Acceptance Rate: 0.4532  ← Your NN performance
Overall acceptance rate: 0.6234                     ← Combined (trad + NN)
```

## 💡 Pro Tips

1. **Baseline first**: Run simulation before training to see baseline
2. **Compare**: Run before/after training to see improvement
3. **Iterate**: If not satisfied, train longer or adjust rewards
4. **Monitor**: Check if acceptance rate improves during training
5. **Save old models**: Rename `model_01_rl.pt` before retraining

## 🚦 Status Indicators

✅ **Good signs during training:**
- Acceptance rate increasing
- Loss decreasing
- Entropy positive (>0)
- Rewards increasing

⚠️ **Warning signs:**
- Acceptance rate stuck or decreasing
- Loss not changing
- Entropy near zero (no exploration)

## Summary

```bash
# Complete workflow (first time)
python3 train_nn_with_rl.py          # Train (2-5 min)
python3 run_nn_translate_simulation.py # Run with RL model

# Or use interactive mode
python3 demo_rl_workflow.py          # Guided workflow
```

**That's it!** The system handles everything else automatically.

---

Need more details? See `RL_TRAINING_README.md`

