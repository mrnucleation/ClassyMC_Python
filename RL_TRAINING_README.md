# Reinforcement Learning for Neural Network Translation Moves

This document describes the reinforcement learning (RL) system for optimizing the neural network used in `NNTranslate` Monte Carlo moves.

## Overview

The RL training system uses **Proximal Policy Optimization (PPO)** to improve the acceptance rate of neural network-guided translation moves. The neural network learns to propose better atom positions by receiving rewards based on:

1. **Move acceptance** (primary reward: +1.0 for accepted, -0.1 for rejected)
2. **Energy improvement** (secondary reward: up to +0.5 for energy decrease)
3. **Move diversity** (exploration bonus: +0.2 for reasonable displacements)

## Files

- **`train_nn_with_rl.py`**: Main RL training script
- **`run_nn_translate_simulation.py`**: Inference script (uses trained model)
- **`demo_rl_workflow.py`**: Interactive demonstration of the workflow
- **`src_python/MC_Move_NNTranslate.py`**: Core NNTranslate implementation (updated to auto-detect RL model)

## Model Priority

The system automatically prioritizes models in this order:

1. **`model_01_rl.pt`** - RL-trained model (created by training script)
2. **`model_01.pt`** - Base/pretrained model (fallback)
3. **Random initialization** - Used if no models exist

This happens automatically when you call `create_and_set_neural_network()` without arguments.

## Quick Start

### Option 1: Interactive Demo

```bash
python3 demo_rl_workflow.py
```

This provides an interactive menu to:
- Run RL training
- Run inference with trained model
- Execute full workflow (training + inference)

### Option 2: Direct Training

```bash
# Train the model (creates model_01_rl.pt)
python3 train_nn_with_rl.py

# Run simulation with trained model
python3 run_nn_translate_simulation.py
```

## RL Training Details

### Training Configuration

The default configuration in `train_nn_with_rl.py`:

```python
n_training_iterations = 20      # Number of training iterations
n_cycles_per_iteration = 50     # MC cycles per iteration
moves_per_cycle = 26            # Total moves (25 traditional + 1 NN on average)
temperature = 0.7               # Reduced temperature
save_interval = 5               # Save checkpoint every N iterations
```

This provides:
- **~1000 MC cycles total** (20 iterations × 50 cycles)
- **~26,000 total moves** (1000 cycles × 26 moves)
- **~1000 NN moves** (1/26 of total moves due to 25:1 ratio)
- **Training time**: ~2-5 minutes on CPU, <1 minute on GPU

### Why 25:1 Ratio?

The 25:1 ratio of traditional to NN moves provides:

1. **Efficiency**: Traditional moves are fast; NN moves are slower
2. **Experience quality**: System equilibrates between NN moves
3. **Stability**: Prevents NN from dominating before it's trained
4. **Data collection**: Still collects sufficient experiences (~1000 per training run)

### PPO Hyperparameters

```python
lr = 1e-4                    # Learning rate
gamma = 0.99                 # Discount factor
clip_epsilon = 0.2           # PPO clipping parameter
entropy_coef = 0.01          # Entropy bonus (exploration)
n_epochs = 4                 # Training epochs per batch
```

## Reward Function

The reward for each move is computed as:

```python
reward = 0.0

# Primary: Acceptance
if accepted:
    reward += 1.0
else:
    reward -= 0.1

# Secondary: Energy improvement
if energy_change < 0:
    reward += 0.5 * min(abs(energy_change) / 10.0, 1.0)

# Tertiary: Exploration bonus
if 0.1 < displacement < 1.5:
    reward += 0.2
```

This encourages:
- High acceptance rates (primary goal)
- Energy minimization (thermodynamic correctness)
- Diverse moves (avoid getting stuck)

## Architecture

### Neural Network

The neural network (`SimpleFeedforwardNet2Feature`) takes 2-channel gaussian field features and outputs log probabilities for each bin:

```
Input: (2, 201, 201, 201) - 2 channels, 201³ bins
       Channel 0: Repulsive field (max of Gaussian)
       Channel 1: Attractive field (sum of modified Gaussians)

Network: 
  - Linear(2 → 15) + LeakyReLU(0.02)
  - Linear(15 → 8) + LeakyReLU(0.02)
  - Linear(8 → 6) + LeakyReLU(0.03)
  - Linear(6 → 1) + LeakyReLU(0.03)

Output: (201, 201, 201) - Log probabilities for each bin
```

### RL Components

1. **Experience Collection** (`RLNNTranslate`):
   - Extends `NNTranslate` to record experiences during moves
   - Stores: features, log_probs, actions, acceptance, energy_change

2. **PPO Trainer**:
   - Computes advantages using rewards
   - Updates policy using clipped surrogate objective
   - Adds entropy bonus for exploration

## Customization

### Adjust Training Duration

For faster training (less accurate):
```python
n_training_iterations = 10      # Fewer iterations
n_cycles_per_iteration = 25     # Fewer cycles
```

For better training (slower):
```python
n_training_iterations = 50      # More iterations
n_cycles_per_iteration = 100    # More cycles
```

### Modify Reward Function

Edit the `compute_reward` method in `PPOTrainer`:

```python
def compute_reward(self, experience):
    reward = 0.0
    
    # Your custom reward logic here
    if experience['accepted']:
        reward += 2.0  # Increase acceptance weight
    
    # Add custom criteria
    if experience['displacement'] > 0.5:
        reward += 0.5  # Reward larger moves
    
    return reward
```

### Change Move Ratio

In `run_rl_training_episode`:

```python
# More NN moves (10:1 ratio)
mc.moveweights = [10.0, 1.0]

# Equal ratio (1:1)
mc.moveweights = [1.0, 1.0]

# More traditional moves (50:1 ratio)
mc.moveweights = [50.0, 1.0]
```

## Monitoring Training

The training script prints:

```
--- Training Iteration 5/20 ---
Collected 42 experiences
Acceptance rate: 0.4286 (42.86%)
Energy change: -2.453821
Training stats:
  Loss: 0.234567
  Policy loss: 0.223456
  Entropy: 0.011234
  Avg reward: 0.8234
  *** New best acceptance rate: 0.4286 ***
  Saved best model to model_01_rl.pt
```

Key metrics:
- **Acceptance rate**: Higher is better (aim for >40%)
- **Loss**: Should decrease over time
- **Entropy**: Should be positive (ensures exploration)
- **Avg reward**: Should increase over time

## Troubleshooting

### Low Acceptance Rate (<20%)

- Increase training iterations
- Adjust reward function to emphasize acceptance more
- Check if base model exists and is reasonable
- Verify energy function is correct

### Training Too Slow

- Reduce `n_training_iterations` or `n_cycles_per_iteration`
- Increase move ratio (more traditional moves)
- Use GPU if available

### Model Not Loading

Check file locations:
```bash
ls -lh model_01*.pt
```

Should see:
```
model_01.pt      # Base model (optional)
model_01_rl.pt   # RL-trained model (created by training)
```

### Acceptance Rate Decreasing

- RL might be overfitting to specific configurations
- Increase entropy coefficient: `entropy_coef = 0.02`
- Use more diverse training configurations
- Add regularization to reward function

## Advanced Usage

### Batch Training on Multiple Systems

```python
# Create different systems for training
systems = [
    create_fcc_with_vacancies(nvacancies=3),
    create_fcc_with_vacancies(nvacancies=5),
    create_vapor_box()
]

# Train on each system in rotation
for box in systems:
    experiences, _, _ = run_rl_training_episode(box, nn_move)
    trainer.train_on_experiences(experiences)
```

### Transfer Learning

```python
# Start from a partially trained model
model, is_rl = load_or_create_model(
    rl_model_path="model_checkpoint_v2.pt",
    base_model_path="model_01.pt"
)
```

### Curriculum Learning

```python
# Start with easy systems, progress to harder ones
temperatures = [1.0, 0.85, 0.7, 0.55]  # Easier to harder

for temp in temperatures:
    box.temperature = temp
    box.beta = 1.0 / temp
    experiences, _, _ = run_rl_training_episode(box, nn_move)
    trainer.train_on_experiences(experiences)
```

## Performance Tips

1. **Use GPU**: Training is 5-10× faster on GPU
   ```python
   device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
   ```

2. **Batch experiences**: Collect more experiences before training
   ```python
   n_cycles_per_iteration = 100  # Collect more data
   ```

3. **Adjust learning rate**: Start higher, decay over time
   ```python
   lr_schedule = [1e-3, 5e-4, 1e-4, 5e-5]
   for i, lr in enumerate(lr_schedule):
       trainer.optimizer.lr = lr
   ```

## Citation

If you use this RL training system in your research, please cite:

```
[Your paper/repository citation here]
```

## License

[Your license here]

## Support

For questions or issues:
- Open an issue on GitHub
- Email: [your-email]
- Documentation: [link]

