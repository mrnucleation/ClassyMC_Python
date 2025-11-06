from src_python.Template_MCMove import MCMove
import numpy as np
from random import random, choice
from src_python.VarPrecision import dp
from src_python.Box_SimpleBox import SimpleBox
from src_python.CoordinateTypes import Displacement
import math
from scipy.special import logsumexp
import torch
import torch.nn as nn
import time
import matplotlib.pyplot as plt
import logging


#Log to a file
logging.basicConfig(level=logging.INFO, filename='NNTranslate.log')
logger = logging.getLogger(__name__)


def time_function(func):
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        print(f"Time taken for {func.__name__}: {end_time - start_time} seconds")
        return result
    
    return wrapper


#=======================================================================
'''
class SimpleFeedforwardNet2Feature(nn.Module):
    def __init__(self):
        super(SimpleFeedforwardNet2Feature, self).__init__()
        self.fc1 = nn.Linear(3, 30)
        self.act1 = nn.LeakyReLU(0.02)
        self.act2 = nn.LeakyReLU(0.02)
        self.fc2 = nn.Linear(30, 20)
        self.fc3 = nn.Linear(20, 10)
        self.fc4 = nn.Linear(10, 1)

    def forward(self, x):
        # Input shape is expected to be (batch_size, 3, D1, D2, D3)
        in_shape = x.shape
        
        # Permute to (batch_size, D1, D2, D3, 3)
        x = x.permute(0, 2, 3, 4, 1)
        
        # Reshape to (batch_size*D1*D2*D3, 3)
        xn = x.reshape(-1, x.shape[-1])

        x1 = self.fc1(xn)
        x1 = self.act1(x1)

        x2 = self.fc2(x1)
        x2 = self.act1(x2)
        
        x3 = self.fc3(x2)
        x3 = self.act2(x3)
        
        out = self.fc4(x3)
        out = self.act2(out)
        
        # Reshape back to (batch_size, D1, D2, D3)
        out = out.reshape(in_shape[0], in_shape[2], in_shape[3], in_shape[4])
        
        # Normalize the output log probability
        out = out - torch.logsumexp(out, dim=(1,2,3), keepdim=True)
        
        return out
'''

import torch
import torch.nn as nn


#=======================================================================
class SimpleFeedforwardNet2Feature(nn.Module):
    def __init__(self):
        super(SimpleFeedforwardNet2Feature, self).__init__()
        #self.fc1 = nn.Linear(2, 35)
        #self.act1 = nn.LeakyReLU(0.02)
        #self.act2 = nn.LeakyReLU(0.02)
        #self.fc2 = nn.Linear(35, 23)
        #self.fc3 = nn.Linear(23, 12)
        #self.fc4 = nn.Linear(12, 1)
        
        self.fc1 = nn.Linear(2, 15)
        self.act1 = nn.LeakyReLU(0.02)
        self.act2 = nn.LeakyReLU(0.03)
        self.fc2 = nn.Linear(15, 8)
        self.fc3 = nn.Linear(8, 6)
        self.fc4 = nn.Linear(6, 1)

    def forward(self, x):
        # Input shape is expected to be (batch_size, 2, D1, D2, D3)
        in_shape = x.shape
        
        # Permute to (batch_size, D1, D2, D3, 2)
        x = x.permute(0, 2, 3, 4, 1)
        
        # Reshape to (batch_size*D1*D2*D3, 2)
        xn = x.reshape(-1, x.shape[-1])

        x1 = self.fc1(xn)
        x1 = self.act1(x1)

        x2 = self.fc2(x1)
        x2 = self.act1(x2)
        
        x3 = self.fc3(x2)
        x3 = self.act2(x3)
        
        out = self.fc4(x3)
        out = self.act2(out)
        
        # Reshape back to (batch_size, D1, D2, D3)
        out = out.reshape(in_shape[0], in_shape[2], in_shape[3], in_shape[4])
        
        # Normalize the output log probability
        out = out - torch.logsumexp(out, dim=(1,2,3), keepdim=True)
        
        return out
#=======================================================================
class NNTranslate(MCMove):
    """
    Performs a Monte Carlo translation move on a randomly selected atom using
    a neural network model to generate position proposals.

    This class extends the standard molecule translation move by:
    1. Gathering neighboring atoms within a cutoff distance (with PBC)
    2. Creating a binned space around the target atom
    3. Computing gaussian fields centered at neighboring atoms
    4. Using a neural network to predict log probabilities for each bin
    5. Sampling a bin using Gumbel's method
    6. Generating a uniform random position within the selected bin
    7. Accounting for bias in the acceptance probability

    The neural network model should be provided separately and implement
    a predict_log_probs method that takes gaussian field features as input.
    """
    #----------------------------------------------------------------------------
    def __init__(self, BoxArray):
        super().__init__()
        # --- Default attributes from Fortran type definition ---
        self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        self.verbose = True
        self.proportional = True

        # --- Neural Network Parameters ---
        # Simulation box definition: 5x5x5 Angstroms
        self.nbins = [201, 201, 201]  # 50x50x50 bins
        self.nbins = torch.tensor(self.nbins, device=self.device)
        
        #Assert the bins are odd
        assert self.nbins[0] % 2 == 1, "Number of bins must be odd"
        assert self.nbins[1] % 2 == 1, "Number of bins must be odd"
        assert self.nbins[2] % 2 == 1, "Number of bins must be odd"
        
        self.middle_bin_index = torch.tensor([self.nbins[0] // 2, self.nbins[1] // 2, self.nbins[2] // 2], device=self.device)
        self.binsize = [0.02, 0.02, 0.02]
        self.binsize_half = [self.binsize[0]/2, self.binsize[1]/2, self.binsize[2]/2]
        self.binsize = torch.tensor(self.binsize, device=self.device)
        self.binsize_half = torch.tensor(self.binsize_half, device=self.device)
 
 
        self.largest_displacement = 0.0       # Largest displacement created by the move

        #compute the translation box
        self.translation_box = torch.tensor([
            (-self.binsize[0] * (self.nbins[0] // 2), self.binsize[0] * (self.nbins[0] // 2)),
            (-self.binsize[1] * (self.nbins[1] // 2), self.binsize[1] * (self.nbins[1] // 2)),
            (-self.binsize[2] * (self.nbins[2] // 2), self.binsize[2] * (self.nbins[2] // 2))
        ], device=self.device)
        print(f"Translation box: {self.translation_box}")
        self.trans_box_low = torch.tensor([
            -self.binsize[0] * (self.nbins[0] // 2),
            -self.binsize[1] * (self.nbins[1] // 2),
            -self.binsize[2] * (self.nbins[2] // 2)
        ], device=self.device)
        self.trans_box_high = torch.tensor([
            self.binsize[0] * (self.nbins[0] // 2),
            self.binsize[1] * (self.nbins[1] // 2),
            self.binsize[2] * (self.nbins[2] // 2)
        ], device=self.device)


        self.trans_box_volume = (self.translation_box[0][1] - self.translation_box[0][0]) * (self.translation_box[1][1] - self.translation_box[1][0]) * (self.translation_box[2][1] - self.translation_box[2][0])
        self.trans_box_volume = self.trans_box_volume.cpu().numpy()
        self.bin_volume = self.binsize[0] * self.binsize[1] * self.binsize[2]
        self.bin_volume = self.bin_volume.cpu().numpy()
        
        self.middle_bin_index = torch.tensor([self.nbins[0] // 2, self.nbins[1] // 2, self.nbins[2] // 2], device=self.device)

        self.neighbor_cutoff = 4.5   # Cutoff for neighbor detection
        self.gaussian_sigma = 1.0    # Sigma for gaussian fields (both channels use same sigma)
        self.gaussian_sigma2 = 0.15   # Sigma for second gaussian field
        self.gaussian_sigma3 = 0.32   # Sigma for third gaussian field

        self.cutoff_rc = 2.5         # Cutoff radius for Behler-Parrinello function
        self.cutoff_d = 0.1          # Cutoff transition width
        self.nn_model = None         # Neural network model (to be set externally)

        # --- Rejection Counters ---
        self.ovlaprej = 0
        self.constrainrej = 0
        self.detailedrej = 0
        self.nnrej = 0  # Rejection due to NN model issues


        n_boxes = len(BoxArray)

        # Initialize box probabilities
        if len(self.boxProb) == 0:
            if n_boxes > 0:
                self.boxProb = np.full(n_boxes, 1.0 / n_boxes, dtype=dp)
            else:
                self.boxProb = np.array([], dtype=dp)

        # Initialize box-specific arrays
        self.boxatmps = np.full(n_boxes, 1e-50, dtype=dp)
        self.boxaccpt = np.zeros(n_boxes, dtype=dp)

        # Initialize arrays for binning and gaussian fields
        self.bin_centers = None
        self.gaussian_fields = None
        self.log_probs = None
        
        # Pre-compute bin centers for the simulation box
        self._initialize_bin_centers()

    #----------------------------------------------------------------------------
    #@time_function
    def _initialize_bin_centers(self):
        """Initialize the bin centers for the simulation box."""
        # Vectorized version: create grid indices using meshgrid
        i_indices = np.arange(self.nbins[0].cpu().numpy())
        j_indices = np.arange(self.nbins[1].cpu().numpy())
        k_indices = np.arange(self.nbins[2].cpu().numpy())
        
        # Create 3D grids of indices
        i_grid, j_grid, k_grid = np.meshgrid(i_indices, j_indices, k_indices, indexing='ij')
        
        # Compute all bin centers at once using broadcasting
        bin_centers_np = np.stack([
            self.translation_box[0][0].cpu().numpy() + self.binsize[0].cpu().numpy() * (i_grid + 0.5),
            self.translation_box[1][0].cpu().numpy() + self.binsize[1].cpu().numpy() * (j_grid + 0.5),
            self.translation_box[2][0].cpu().numpy() + self.binsize[2].cpu().numpy() * (k_grid + 0.5)
        ], axis=-1)

        # Convert to PyTorch tensor and move to device
        self.bin_centers = torch.tensor(bin_centers_np, dtype=torch.float32, device=self.device)
        #print(f"Bin centers: {self.bin_centers}")
    #----------------------------------------------------------------------------
    #@time_function
    def _cutoff_function(self, r, rc=None, d=None):
        """Behler-Parrinello cutoff function for numpy arrays"""
        if rc is None:
            rc = self.cutoff_rc
        if d is None:
            d = self.cutoff_d
            
        result = np.ones_like(r)
        
        # For r >= (rc + d), set to 0
        mask_high = r >= (rc + d)
        result[mask_high] = 0.0
        
        # For (rc - d) < r < (rc + d), apply cosine cutoff
        mask_mid = (r > (rc - d)) & (r < (rc + d))
        result[mask_mid] = 0.5 * (np.cos(np.pi * (r[mask_mid] - rc + d) / (2 * d)) + 1.0)
        
        # For r <= (rc - d), already set to 1 by default

        return result

    #----------------------------------------------------------------------------
    def _cutoff_function_torch(self, r, rc=None, d=None):
        """Behler-Parrinello cutoff function for PyTorch tensors"""
        if rc is None:
            rc = self.cutoff_rc
        if d is None:
            d = self.cutoff_d

        result = torch.ones_like(r, device=self.device)

        # For r >= (rc + d), set to 0
        mask_high = r >= (rc + d)
        result[mask_high] = 0.0

        # For (rc - d) < r < (rc + d), apply cosine cutoff
        mask_mid = (r > (rc - d)) & (r < (rc + d))
        result[mask_mid] = 0.5 * (torch.cos(torch.pi * (r[mask_mid] - rc + d) / (2 * d)) + 1.0)

        # For r <= (rc - d), already set to 1 by default

        return result

    #----------------------------------------------------------------------------
    def set_neural_network(self, nn_model):
        """
        Set the neural network model for position prediction.

        Args:
            nn_model: PyTorch neural network model
        """
        self.nn_model = nn_model
        self.nn_model.to(self.device)
    #----------------------------------------------------------------------------
    def create_and_set_neural_network(self, model_path=None):
        """
        Create and set an instance of SimpleFeedforwardNet2Feature.
        
        Args:
            model_path: Optional path to load pretrained model weights
        """
        self.nn_model = SimpleFeedforwardNet2Feature()
        if model_path is not None:
            self.nn_model.load_state_dict(torch.load(model_path, map_location='cpu'))
        self.nn_model.eval()  # Set to evaluation mode
        self.nn_model.to(self.device)

    #----------------------------------------------------------------------------
    #@time_function
    def full_move(self, trial_box: SimpleBox, sampling):
        """
        Performs one neural network-based atom translation MC move.
        Returns a boolean indicating if the move was accepted.
        """
        box_idx = trial_box.boxID - 1  # Convert 1-based ID to 0-based index

        self.atmps += 1.0
        self.boxatmps[box_idx] += 1.0



        # Select a random atom
        
        target_mol = trial_box.pick_random_molecule()
        
        atom_indices = target_mol['atomindicies'][0]  # Extract the array of atom indices
        target_atom_idx = atom_indices[0]  # Get the first atom index
        target_mol_idx = target_mol['mol_index']
        target_pos = trial_box.atoms[target_atom_idx]

        # Get molecule information for the target atom
        mol_idx = target_mol_idx
        mol_type = target_mol['mol_type']

        # --- Gather neighboring atoms within cutoff ---
        neighbor_positions = self._gather_neighbors(
            trial_box, target_pos, target_mol_idx
        )
        

        #if len(neighbor_positions) == 0:
        #    # Generate a random uniform position within the translation box
        #    delta = np.random.uniform(self.trans_box_low, self.trans_box_high, size=3)
        #    forward_prob = -np.log(self.trans_box_volume)
        #else:
        features = self._compute_2channel_features(trial_box, target_pos, neighbor_positions)
        with torch.no_grad():
            features_tensor = features.unsqueeze(0)  # Add batch dimension (features is already a tensor)
            self.log_probs = self.nn_model(features_tensor).squeeze(0)  # Remove batch dimension
            self.log_probs = self.log_probs.cpu().numpy()
            
        
        # --- Sample bin using Gumbel's method ---
        selected_bin_idx, forward_prob = self._gumbel_sample(self.log_probs)
        #forward_prob = forward_prob - np.log(self.bin_volume)

        # --- Generate position within selected bin ---
        delta, bin_delta = self._generate_position_in_bin(selected_bin_idx)
        new_pos_raw = target_pos + delta
        new_pos = trial_box.boundary(new_pos_raw)


        # Create displacement object
        disp = Displacement(
            molType=mol_type,
            molIndx=mol_idx,
            atmIndicies=np.array([target_atom_idx]),
            newPositions=new_pos.reshape(1, 3)  # Single atom move
        )

        # --- Check constraints and calculate energy ---
        if not trial_box.check_constraint(disp):
            self.constrainrej += 1
            return False

        e_inter, e_intra, accept = trial_box.compute_energy_delta(
            disp,
            self.tempList,
            self.tempNnei,
            computeintra=False
        )
        if not accept:
            self.ovlaprej += 1
            return False

        e_diff = e_inter + e_intra

        if not trial_box.check_post_energy(disp, e_diff):
            self.constrainrej += 1
            return False

        # --- Calculate forward and reverse probabilities for bias ---
        reverse_prob = self._calculate_reverse_probability(trial_box, new_pos_raw, target_pos, bin_delta, target_mol_idx)
        
       
        # Include bias in the acceptance probability
        bias_correction = reverse_prob - forward_prob

        # --- Accept/Reject ---
        accept = sampling.make_decision(trial_box, e_diff, disp, log_prob=bias_correction)
        if accept:
            self.accpt += 1.0
            self.boxaccpt[box_idx] += 1.0
            trial_box.update_energy(e_diff, e_inter, e_intra)
            trial_box.update_position(disp)
            
            self.largest_displacement = max(self.largest_displacement, np.linalg.norm(delta))
        else:
            self.detailedrej += 1
            features_np = features.cpu().numpy()
            '''
            logger.info(f"Number of neighboring atoms: {len(neighbor_positions)}")
            logger.info(f"Selected bin index: {selected_bin_idx}")
            logger.info(f"Gaussian Features shape: {features_np.shape}")
            logger.info(f"Gaussian Features min: {np.min(features_np)}")
            logger.info(f"Gaussian Features max: {np.max(features_np)}")
            logger.info(f"Gaussian Features mean: {np.mean(features_np)}")
            logger.info(f"Gaussian Features std: {np.std(features_np)}")           
            logger.info(f"Log probabilities shape: {self.log_probs.shape}")
            logger.info(f"Log probabilities type: {self.log_probs.dtype}")
            logger.info(f"Log probabilities min: {np.min(self.log_probs)}")
            logger.info(f"Log probabilities max: {np.max(self.log_probs)}")
            logger.info(f"Log probabilities mean: {np.mean(self.log_probs)}")
            logger.info(f"Log probabilities std: {np.std(self.log_probs)}")           
            logger.info(f"Forward probability: {forward_prob:8.5f}")
            logger.info(f"Reverse probability: {reverse_prob:8.5f}")
            logger.info(f"Bias correction: {bias_correction:8.5f}")
            logger.info(f"e_diff: {e_diff:8.5f}")
            logger.info(f"Accepted: {accept}")
            '''
        return accept

    #----------------------------------------------------------------------------
    #@time_function
    def _gather_neighbors(self, trial_box, target_pos, target_mol_idx):
        """
        Gather neighboring atoms within cutoff distance, handling periodic boundaries.
        Optimized to minimize CPU/GPU transfers.
        """
        atoms = trial_box.atoms[np.where(trial_box.MolIndx != target_mol_idx)]
        # Calculate displacement vectors from target to all atoms
        displacements = atoms[:, :3] - target_pos.reshape(1, 3)

        # Apply periodic boundary conditions using minimum image convention
        displacements = trial_box.boundary(displacements)
        
        # Convert to torch tensor once (no redundant .to(device) call)
        displacements = torch.tensor(displacements, dtype=torch.float32, device=self.device)
        
        # Calculate distances and filter by neighbor_cutoff
        distances = torch.sqrt(torch.sum(displacements**2, dim=1))
        mask = distances < self.neighbor_cutoff
        
        displacements = displacements[mask]

        return displacements
    #----------------------------------------------------------------------------
    def _compute_2channel_features(self, trial_box, target_pos, neighbor_positions):
        """
        Compute 2-channel gaussian features with optimized memory usage.
        Returns features for neural network input with shape (2, nbins[0], nbins[1], nbins[2]).
        
        Channel 0: Max of exp(-rsq/(2*sigma2^2)) * cutoff
        Channel 1: Sum of (exp(-(r-lj_sigma)^2/(2*sigma2^2)) - exp(-rsq/(2*sigma2^2))) * cutoff
        
        Optimized by:
        - Processing atoms in chunks to reduce peak memory
        - Avoiding redundant cutoff function calls
        - Reusing computed distances
        """
        atoms = neighbor_positions

        # Check if there are no neighboring atoms
        if len(atoms) == 0:
            logger.info(f"No neighboring atoms found")
            return torch.zeros((2, self.nbins[0], self.nbins[1], self.nbins[2]), 
                             dtype=torch.float32, device=self.device)

        atoms_tensor = atoms
        n_atoms = atoms_tensor.shape[0]
        
        # Initialize 2-channel feature array
        gauss_features = torch.zeros((2, self.nbins[0], self.nbins[1], self.nbins[2]), 
                                    dtype=torch.float32, device=self.device)
        
        # Process in chunks to reduce memory usage
        chunk_size = min(10, n_atoms)  # Process up to 10 atoms at a time
        
        # Precompute constants
        inv_2sigma2_sq = 1.0 / (2.0 * self.gaussian_sigma2 * self.gaussian_sigma2)
        lj_mu = 2**(1.0/6.0)
        
        for i in range(0, n_atoms, chunk_size):
            end_idx = min(i + chunk_size, n_atoms)
            atoms_chunk = atoms_tensor[i:end_idx]
            
            # Expand dimensions for broadcasting
            bin_centers_expanded = self.bin_centers.unsqueeze(0)  # (1, nbins[0], nbins[1], nbins[2], 3)
            atoms_expanded = atoms_chunk.unsqueeze(1).unsqueeze(1).unsqueeze(1)  # (chunk, 1, 1, 1, 3)

            # Calculate distances
            diff = bin_centers_expanded - atoms_expanded  # (chunk, nbins[0], nbins[1], nbins[2], 3)
            rsq = torch.sum(diff**2, dim=-1)  # (chunk, nbins[0], nbins[1], nbins[2])
            r = torch.sqrt(rsq + 1e-10)  # Add small epsilon for numerical stability

            # Calculate cutoff function once
            cutoff = self._cutoff_function_torch(r, rc=self.cutoff_rc, d=self.cutoff_d)

            # Calculate exponentials with cutoff functions
            e1 = torch.exp(-rsq * inv_2sigma2_sq) * cutoff
            e2 = torch.exp(-((r - lj_mu)**2) * inv_2sigma2_sq) * cutoff

            e1 = torch.exp(-rsq/(2.0*self.gaussian_sigma*self.gaussian_sigma)) * cutoff  # (natoms-1, nbins[0], nbins[1], nbins[2])
            e2 = torch.exp(-rsq/(2.0*self.gaussian_sigma2*self.gaussian_sigma2)) * cutoff  # (natoms-1, nbins[0], nbins[1], nbins[2])
            e3 = torch.exp(-(r-2**(1.0/6.0))**2/(2.0*self.gaussian_sigma3*self.gaussian_sigma3)) * cutoff  # (natoms-1, nbins[0], nbins[1], nbins[2])
            
            # Max over all atoms for channel 0
            gauss_features[0, :, :, :] = torch.maximum(gauss_features[0, :, :, :],torch.max(e1, dim=0)[0]) 
            # Sum over all atoms for channel 1
            gauss_features[1, :, :, :] = torch.sum(e2-2.0*e3, dim=0)

            # Accumulate features
            # Channel 0: Max of e1
            #gauss_features[0, :, :, :] = torch.maximum(gauss_features[0, :, :, :],torch.max(e1, dim=0)[0])
            # Channel 1: Sum of (e2 - e1)
            #gauss_features[1, :, :, :] += torch.sum(e2 - e1, dim=0)

        return gauss_features

    #----------------------------------------------------------------------------
    #@time_function
    def _gumbel_sample(self, log_probs):
        """
        Sample from categorical distribution using Gumbel's method.
        Optimized to work directly with numpy arrays efficiently.
        log_probs should be a 3D array with shape (nbins[0], nbins[1], nbins[2]).
        """
        

        # Add Gumbel noise to log probabilities
        gumbel_noise = np.random.gumbel(size=log_probs.shape)
        perturbed_logits = log_probs + gumbel_noise

        # Return 3D index of maximum as tuple (i, j, k)
        flat_idx = np.argmax(perturbed_logits)
        bin_idx = np.unravel_index(flat_idx, log_probs.shape)
        probval = log_probs[bin_idx]
        
        return bin_idx, probval

    #----------------------------------------------------------------------------
   
    #@time_function
    def _generate_position_in_bin(self, bin_idx):
        """
        Generate a uniform random position within the selected bin.
        bin_idx is a tuple (i, j, k) representing the 3D bin indices.
        Optimized to minimize GPU/CPU transfers.
        """
        i, j, k = bin_idx
        bin_center = self.bin_centers[i, j, k]
        
        # Generate random position within bin on GPU, then convert once
        random_offset = (torch.rand(3, device=self.device) - 0.5) * self.binsize
        
        bin_delta = torch.tensor(bin_idx, device=self.device) - self.middle_bin_index
        new_pos = bin_center + random_offset
        
        
        # Single CPU transfer at the end
        return new_pos.cpu().numpy(), bin_delta

    #----------------------------------------------------------------------------
    #@time_function
    def _calculate_forward_probability(self, selected_bin_idx):
        """
        Calculate the forward probability for the bias correction.
        This is the probability of selecting the chosen bin.
        selected_bin_idx is a tuple (i, j, k).
        """
        i, j, k = selected_bin_idx
        #bin_prob = np.exp(self.log_probs[i, j, k])


        # Multiply by probability of selecting position within bin
        # (uniform within bin, so 1/volume of bin)
        #bin_volume = self.binsize[0] * self.binsize[1] * self.binsize[2]
        #position_prob = 1.0 / bin_volume

        return self.log_probs[i, j, k] - np.log(self.bin_volume)

    #----------------------------------------------------------------------------
    #@time_function
    def _calculate_reverse_probability(self, trial_box, new_pos, old_pos, bin_delta, target_mol_idx):
        """
        Calculate the reverse probability for the bias correction.
        This computes the probability of proposing old_pos when centered at new_pos.
        Optimized to minimize CPU/GPU transfers.
        """
        # Gather neighbors around the new position
        neighbor_positions = self._gather_neighbors(
            trial_box, new_pos, target_mol_idx
        )
        
        # Compute features centered at new position
        reverse_features = self._compute_2channel_features(trial_box, new_pos, neighbor_positions)
        
        # Get neural network predictions for reverse direction
        with torch.no_grad():
            features_tensor = reverse_features.unsqueeze(0)  # Add batch dimension
            reverse_log_probs = self.nn_model(features_tensor).squeeze(0)
            reverse_log_probs = reverse_log_probs.cpu().numpy()
        
        # Find which bin the old position would fall into
        # Calculate relative position from new_pos to old_pos
        
        rel_pos = old_pos - new_pos
        #rel_pos = trial_box.boundary(rel_pos)
        
        # Find which bin this relative position corresponds to (vectorized)
        # Convert to torch tensor for vectorized operations
        #rel_pos_tensor = torch.tensor(rel_pos, dtype=torch.float32, device=self.device)
        
        # Calculate bin indices (vectorized)
        #bin_coords = (rel_pos_tensor - self.trans_box_low) / self.binsize
        #bin_indices = torch.floor(bin_coords).long()
        
        bin_indices = -bin_delta + self.middle_bin_index
        
        
        #Check if the bin indices are within the bounds of the bins, raise an error if not
        if torch.any(bin_indices < 0) or torch.any(bin_indices > self.nbins):
            print(f"New position: {new_pos}")
            print(f"Old position: {old_pos}")
            print(f"Relative position: {rel_pos}")
            print(f"Bin indices: {bin_indices}")
            raise ValueError("Bin indices are out of bounds")
        
        # Convert to CPU numpy for indexing (single transfer)
        bin_indices = bin_indices.cpu().numpy()
        
        return reverse_log_probs[bin_indices[0], bin_indices[1], bin_indices[2]]# - np.log(self.bin_volume)
    #----------------------------------------------------------------------------
    def screenout(self):
        """Prints information about the move."""
        print(f"Neural Network Translation Acceptance Rate: {self.get_accept_rate():15.8f}")   
        print(f"Largest displacement: {self.largest_displacement:8.3f} Å")
    #----------------------------------------------------------------------------
    def maintenance(self):
        """
        No maintenance needed for NN-based moves.
        """
        pass

    #----------------------------------------------------------------------------
    def prologue(self):
        """Prints information at the start of a simulation block."""
        print(f"(Neural Network Translation)")
        print(f"Number of bins: {self.nbins}")
        print(f"Bin sizes: {self.binsize}")
        print(f"Neighbor cutoff: {self.neighbor_cutoff:8.3f} Å")
        print(f"Gaussian sigma 1: {self.gaussian_sigma:8.3f} Å")
        print(f"Gaussian sigma 2: {self.gaussian_sigma2:8.3f} Å")
        print(f"Cutoff radius: {self.cutoff_rc:8.3f} Å")
        print(f"Translation box: {self.translation_box}")
        if self.nn_model is None:
            print("Warning: No neural network model set!")

    #----------------------------------------------------------------------------
    def epilogue(self):
        """Prints summary statistics at the end of a simulation block."""
        print()
        print(f"Neural Network Translation Moves Accepted:  {round(self.accpt):15d}")
        print(f"Neural Network Translation Moves Attempted: {round(self.atmps):15d}")

        accpt_rate = self.get_accept_rate()
        print(f"Neural Network Translation Acceptance Rate: {accpt_rate:15.8f}")

        if self.verbose:
            print(f"Neural Network Translation, Rejections due to overlap:         {self.ovlaprej:15d}")
            print(f"Neural Network Translation, Rejections due to constraint:      {self.constrainrej:15d}")
            print(f"Neural Network Translation, Rejections due to detailed balance:{self.detailedrej:15d}")
            print(f"Neural Network Translation, Rejections due to NN issues:       {self.nnrej:15d}")

    #----------------------------------------------------------------------------
    def process_io(self, line: str):
        """
        Parses a line from an input file to set parameters.
        """
        parts = line.strip().lower().split()
        if len(parts) < 5:
            print(f"Invalid input line: {line.strip()}")
            raise IOError("Invalid input line format")

        command, value_str = parts[3], parts[4]
        try:
            if command == "neighbor_cutoff":
                self.neighbor_cutoff = float(value_str)
            elif command == "gaussian_sigma":
                self.gaussian_sigma = float(value_str)
            elif command == "gaussian_sigma2":
                self.gaussian_sigma2 = float(value_str)
            elif command == "cutoff_rc":
                self.cutoff_rc = float(value_str)
            elif command == "cutoff_d":
                self.cutoff_d = float(value_str)
            elif command == "proportional":
                self.proportional = value_str.lower() in ['true', '.true.', 't', '1']
            elif command == "updatefreq":
                self.maintFreq = int(value_str)
            else:
                return -1  # Command not recognized
        except (ValueError, IndexError):
            return -1  # Error parsing value
        return 0
#=======================================================================
