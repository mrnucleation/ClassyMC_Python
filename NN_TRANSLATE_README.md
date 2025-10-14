# Neural Network-Based Translation Move (NNTranslate)

This document describes the new `NNTranslate` class, which extends the standard molecule translation move in ClassyMC to use neural network predictions for generating new atom positions.

## Overview

The `NNTranslate` move is an extension of the existing `MolTranslate` class that uses a neural network model to propose new positions for atoms during Monte Carlo simulations. The key innovation is that instead of using random displacements, the move:

1. **Gathers neighboring atoms** within a specified cutoff distance (default: 10 Å)
2. **Creates a binned 3D space** around the target atom
3. **Computes gaussian fields** centered at each neighboring atom
4. **Uses a neural network** to predict log-probabilities for each bin
5. **Samples a bin** using Gumbel's method
6. **Generates a uniform random position** within the selected bin
7. **Accounts for bias** in the acceptance probability calculation

## Files Added

- `src_python/MC_Move_NNTranslate.py` - Main implementation of the NNTranslate class
- `src_python/NN_Model_Interface.py` - Abstract base class and example for neural network models
- `src_python/Example_NN_Usage.py` - Usage examples and integration guide

## Key Features

### Periodic Boundary Conditions
The move properly handles periodic boundary conditions when gathering neighbors and computing distances.

### Configurable Parameters
- `neighbor_cutoff`: Distance to search for neighboring atoms (default: 10.0 Å)
- `bin_size`: Size of each bin in the 3D grid (default: 0.5 Å)
- `n_bins_per_dim`: Number of bins per dimension (default: 20, giving 8000 total bins)
- `gaussian_sigma`: Standard deviation for gaussian fields (default: 2.0 Å)

### Bias Correction
Since this is a biased move (non-uniform proposal distribution), the acceptance probability correctly accounts for the forward and reverse proposal probabilities.

### Fallback Mechanism
If the neural network fails or no neighbors are found, the move falls back to a simple random displacement.

## Usage

### Basic Setup

```python
from src_python.MC_Move_NNTranslate import NNTranslate
from src_python.NN_Model_Interface import NNPositionPredictor

# Create your neural network model
class MyModel(NNPositionPredictor):
    def predict_log_probs(self, features):
        # Your neural network prediction logic here
        return log_probabilities

# Set up the move
nn_move = NNTranslate(box_array)
nn_move.set_neural_network(MyModel())

# Configure parameters
nn_move.neighbor_cutoff = 10.0
nn_move.bin_size = 0.5
nn_move.n_bins_per_dim = 20
```

### Integration with Existing Code

The `NNTranslate` class follows the same interface as other MC moves in ClassyMC:

```python
# In your simulation loop
accepted = nn_move.full_move(trial_box, sampling_method)

# At the beginning of simulation blocks
nn_move.prologue()

# At the end of simulation blocks
nn_move.epilogue()
```

### Input File Configuration

You can configure the move parameters through input files using the standard ClassyMC format:

```
mc move nntranslate neighbor_cutoff 10.0
mc move nntranslate bin_size 0.5
mc move nntranslate n_bins_per_dim 20
mc move nntranslate gaussian_sigma 2.0
```

## Neural Network Model Interface

Your neural network model should inherit from `NNPositionPredictor` and implement the `predict_log_probs` method:

```python
from src_python.NN_Model_Interface import NNPositionPredictor
import numpy as np

class MyNeuralNetwork(NNPositionPredictor):
    def predict_log_probs(self, features: np.ndarray) -> np.ndarray:
        """
        Args:
            features: Shape (n_bins, n_features)
                     Features include gaussian fields, relative positions, and atom types

        Returns:
            log_probs: Shape (n_bins,) - Unnormalized log probabilities
        """
        # Your model prediction code here
        return self.model.predict(features)
```

## Feature Representation

The features passed to the neural network include:

1. **Gaussian fields**: For each bin, the gaussian contribution from each neighboring atom
2. **Relative positions**: Position vectors from the target atom to each neighbor
3. **Atom types**: Type information for each neighboring atom

These features are concatenated and passed to the neural network as a 2D array of shape `(n_bins, n_features)`.

## Performance Considerations

- **Memory usage**: The binning approach creates a 3D grid (default 20×20×20 = 8000 bins)
- **Computation time**: Gaussian field computation scales with the number of neighbors
- **Neural network**: The NN prediction time depends on your model complexity

## Differences from Standard Translation

1. **Proposal distribution**: Non-uniform (learned from data) vs uniform random
2. **Computational cost**: Higher due to neighbor search and NN evaluation
3. **Acceptance probability**: Includes bias correction terms
4. **Applicability**: Best suited for systems where local structure matters

## Example Output

```
(Neural Network Translation)
Neighbor cutoff:   10.000 Å
Bin size:    0.500 Å
Bins per dimension: 20
Gaussian sigma:    2.000 Å

Neural Network Translation Moves Accepted:        850
Neural Network Translation Moves Attempted:      1000
Neural Network Translation Acceptance Rate:   85.000
Neural Network Translation, Rejections due to overlap:        120
Neural Network Translation, Rejections due to constraint:       15
Neural Network Translation, Rejections due to detailed balance: 10
Neural Network Translation, Rejections due to NN issues:         5
```

## Troubleshooting

### Common Issues

1. **No neighbors found**: The move will fall back to random displacement
2. **Neural network errors**: Check your model's input/output format
3. **Memory issues**: Reduce `n_bins_per_dim` for smaller grids
4. **Poor acceptance**: Adjust `bin_size` or neural network training

### Debug Tips

- Enable verbose output to see detailed statistics
- Use the example model (`create_example_model()`) for testing
- Check that your neural network returns the correct shape of log probabilities

## Future Enhancements

Potential improvements to consider:

1. **Adaptive binning**: Dynamically adjust bin size based on local density
2. **Multiple atom moves**: Extend to move entire molecules
3. **Parallel evaluation**: Batch neural network predictions
4. **Training integration**: Built-in training loop for the neural network
5. **Feature engineering**: Additional structural descriptors

## References

This implementation is based on the existing ClassyMC framework and extends the standard molecule translation move with neural network-guided position proposals, similar to approaches used in enhanced sampling methods and machine learning-accelerated molecular dynamics.
