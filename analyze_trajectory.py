#!/usr/bin/env python3
"""
Trajectory Analysis Script for Monte Carlo Simulations

This script computes the Root Mean Squared Displacement (RMSD) of each atom
in the simulation as a function of time step and plots the results.

The script reads trajectory data from a Gromacs .gro file and computes:
- RMSD for each individual atom
- Average RMSD across all atoms
- Standard deviation of RMSD

Usage:
    python analyze_trajectory.py [trajectory_file]
    
Default trajectory file: trajectory.gro
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from ase.io import read

def read_gromacs_custom(filename):
    """
    Custom Gromacs .gro file reader that handles various formatting issues.
    
    Parameters:
    -----------
    filename : str
        Path to the Gromacs trajectory file
        
    Returns:
    --------
    list of tuples
        Each tuple contains (positions, box_vectors) for a frame
        positions: numpy array of shape (n_atoms, 3)
        box_vectors: numpy array of shape (3,) for cubic box or (9,) for general box
    """
    frames = []
    
    with open(filename, 'r') as f:
        while True:
            # Read title line
            title = f.readline()
            if not title:
                break
            
            # Read number of atoms
            n_atoms_line = f.readline()
            if not n_atoms_line:
                break
            
            try:
                n_atoms = int(n_atoms_line.strip())
            except ValueError:
                print(f"Warning: Could not parse number of atoms: {n_atoms_line}")
                break
            
            # Read atom positions
            positions = []
            for i in range(n_atoms):
                line = f.readline()
                if not line:
                    break
                
                # Parse Gromacs format (flexible parsing)
                # Format: resnum resname atomname atomnum x y z [vx vy vz]
                # Positions are in nm, we'll convert to Angstrom
                try:
                    parts = line.split()
                    if len(parts) >= 6:
                        # Take the last 3 or 6 values (positions, and optionally velocities)
                        # The positions are typically at indices 3, 4, 5 (after resnum, resname, atomname, atomnum)
                        # But with flexible parsing, they're the 4th, 5th, and 6th columns (0-indexed: 3, 4, 5)
                        x = float(parts[-6]) * 10.0  # Convert nm to Angstrom
                        y = float(parts[-5]) * 10.0
                        z = float(parts[-4]) * 10.0
                        positions.append([x, y, z])
                    else:
                        print(f"Warning: Line has insufficient columns: {line}")
                        continue
                except (ValueError, IndexError) as e:
                    print(f"Warning: Could not parse atom line: {line.strip()} - {e}")
                    continue
            
            # Read box vectors
            box_line = f.readline()
            if box_line:
                try:
                    box_parts = [float(x) for x in box_line.split()]
                    box_vectors = np.array(box_parts) * 10.0  # Convert nm to Angstrom
                except ValueError:
                    box_vectors = np.array([0.0, 0.0, 0.0])
            else:
                box_vectors = np.array([0.0, 0.0, 0.0])
            
            if len(positions) == n_atoms:
                frames.append((np.array(positions), box_vectors))
            else:
                print(f"Warning: Expected {n_atoms} atoms but read {len(positions)}")
    
    return frames

def read_trajectory(filename='trajectory.gro'):
    """
    Read trajectory from a Gromacs .gro file or XYZ file.
    
    Parameters:
    -----------
    filename : str
        Path to the trajectory file
        
    Returns:
    --------
    list of positions
        List of numpy arrays with atomic positions at each time step
    """
    print(f"Reading trajectory from {filename}...")
    try:
        if filename.endswith('.gro'):
            # Use custom Gromacs reader
            frames = read_gromacs_custom(filename)
            print(f"Successfully read {len(frames)} frames from Gromacs file")
            # Convert to list of position arrays
            trajectory = [frame[0] for frame in frames]
            return trajectory
        else:
            # Try ASE for other formats (XYZ, etc.)
            trajectory_ase = read(filename, index=':', format=None)
            print(f"Successfully read {len(trajectory_ase)} frames")
            # Convert to list of position arrays
            trajectory = [atoms.get_positions() for atoms in trajectory_ase]
            return trajectory
    except FileNotFoundError:
        print(f"ERROR: Trajectory file '{filename}' not found!")
        sys.exit(1)
    except Exception as e:
        print(f"ERROR: Failed to read trajectory: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

def compute_rmsd(trajectory):
    """
    Compute Root Mean Squared Displacement (RMSD) for each atom as a function of time.
    
    RMSD is calculated relative to the initial configuration (frame 0).
    
    Parameters:
    -----------
    trajectory : list of numpy.ndarray
        List of atomic position arrays at each time step, each shape (n_atoms, 3)
        
    Returns:
    --------
    rmsd_per_atom : numpy.ndarray
        RMSD for each atom at each time step, shape (n_frames, n_atoms)
    rmsd_avg : numpy.ndarray
        Average RMSD across all atoms at each time step, shape (n_frames,)
    rmsd_std : numpy.ndarray
        Standard deviation of RMSD at each time step, shape (n_frames,)
    """
    n_frames = len(trajectory)
    n_atoms = len(trajectory[0])
    
    print(f"Computing RMSD for {n_atoms} atoms over {n_frames} frames...")
    
    # Get initial positions (reference configuration)
    initial_positions = trajectory[0]
    
    # Initialize RMSD arrays
    rmsd_per_atom = np.zeros((n_frames, n_atoms))
    
    # Compute RMSD for each frame
    for frame_idx, current_positions in enumerate(trajectory):
        # Compute displacement for each atom
        displacement = current_positions - initial_positions
        
        # Compute squared displacement for each atom
        squared_displacement = np.sum(displacement**2, axis=1)
        
        # RMSD for each atom at this frame
        rmsd_per_atom[frame_idx, :] = np.sqrt(squared_displacement)
    
    # Compute average and standard deviation across all atoms
    rmsd_avg = np.mean(rmsd_per_atom, axis=1)
    rmsd_std = np.std(rmsd_per_atom, axis=1)
    
    print(f"RMSD computation complete!")
    print(f"Final average RMSD: {rmsd_avg[-1]:.6f}")
    print(f"Final RMSD std dev: {rmsd_std[-1]:.6f}")
    
    return rmsd_per_atom, rmsd_avg, rmsd_std

def compute_msd(trajectory):
    """
    Compute Mean Squared Displacement (MSD) for each atom as a function of time.
    
    MSD is calculated relative to the initial configuration (frame 0).
    
    Parameters:
    -----------
    trajectory : list of numpy.ndarray
        List of atomic position arrays at each time step, each shape (n_atoms, 3)
        
    Returns:
    --------
    msd_per_atom : numpy.ndarray
        MSD for each atom at each time step, shape (n_frames, n_atoms)
    msd_avg : numpy.ndarray
        Average MSD across all atoms at each time step, shape (n_frames,)
    msd_std : numpy.ndarray
        Standard deviation of MSD at each time step, shape (n_frames,)
    """
    n_frames = len(trajectory)
    n_atoms = len(trajectory[0])
    
    print(f"Computing MSD for {n_atoms} atoms over {n_frames} frames...")
    
    # Get initial positions (reference configuration)
    initial_positions = trajectory[0]
    
    # Initialize MSD arrays
    msd_per_atom = np.zeros((n_frames, n_atoms))
    
    # Compute MSD for each frame
    for frame_idx, current_positions in enumerate(trajectory):
        # Compute displacement for each atom
        displacement = current_positions - initial_positions
        
        # Compute squared displacement for each atom
        squared_displacement = np.sum(displacement**2, axis=1)
        
        # MSD for each atom at this frame
        msd_per_atom[frame_idx, :] = squared_displacement
    
    # Compute average and standard deviation across all atoms
    msd_avg = np.mean(msd_per_atom, axis=1)
    msd_std = np.std(msd_per_atom, axis=1)
    
    print(f"MSD computation complete!")
    print(f"Final average MSD: {msd_avg[-1]:.6f}")
    print(f"Final MSD std dev: {msd_std[-1]:.6f}")
    
    return msd_per_atom, msd_avg, msd_std

def plot_rmsd(rmsd_per_atom, rmsd_avg, rmsd_std, output_file='rmsd_analysis.png'):
    """
    Create comprehensive RMSD plots.
    
    Parameters:
    -----------
    rmsd_per_atom : numpy.ndarray
        RMSD for each atom at each time step
    rmsd_avg : numpy.ndarray
        Average RMSD across all atoms at each time step
    rmsd_std : numpy.ndarray
        Standard deviation of RMSD at each time step
    output_file : str
        Output filename for the plot
    """
    print(f"Creating RMSD plots...")
    
    n_frames, n_atoms = rmsd_per_atom.shape
    time_steps = np.arange(n_frames)
    
    # Set up the plotting style
    sns.set_style("whitegrid")
    sns.set_context("paper", font_scale=1.2)
    
    # Create figure with subplots
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('Root Mean Squared Displacement (RMSD) Analysis', fontsize=16, fontweight='bold')
    
    # Plot 1: Individual atom RMSD (top left)
    ax1 = axes[0, 0]
    # Plot a subset of atoms if there are too many (for readability)
    n_plot = min(10, n_atoms)
    colors = plt.cm.viridis(np.linspace(0, 1, n_plot))
    for i in range(n_plot):
        ax1.plot(time_steps, rmsd_per_atom[:, i], alpha=0.7, linewidth=1.5, 
                color=colors[i], label=f'Atom {i+1}')
    ax1.set_xlabel('Time Step')
    ax1.set_ylabel('RMSD (Å)')
    ax1.set_title(f'RMSD for Individual Atoms (showing {n_plot}/{n_atoms})')
    if n_plot <= 10:
        ax1.legend(loc='best', fontsize=8)
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Average RMSD with error bars (top right)
    ax2 = axes[0, 1]
    ax2.plot(time_steps, rmsd_avg, 'b-', linewidth=2, label='Mean RMSD')
    ax2.fill_between(time_steps, rmsd_avg - rmsd_std, rmsd_avg + rmsd_std, 
                     alpha=0.3, color='blue', label='±1 Std Dev')
    
    #Restrict y-axis to 0-10
    ax2.set_ylim(0, 5)
    
    ax2.set_xlabel('Time Step')
    ax2.set_ylabel('RMSD (Å)')
    ax2.set_title('Average RMSD Over All Atoms')
    ax2.legend(loc='best')
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: Heatmap of RMSD for all atoms (bottom left)
    ax3 = axes[1, 0]
    # Downsample if there are too many atoms for visualization
    if n_atoms > 50:
        step = n_atoms // 50
        rmsd_plot = rmsd_per_atom[:, ::step].T
        atom_labels = np.arange(0, n_atoms, step)
    else:
        rmsd_plot = rmsd_per_atom.T
        atom_labels = np.arange(n_atoms)
    
    im = ax3.imshow(rmsd_plot, aspect='auto', cmap='hot', interpolation='nearest',
                    extent=[0, n_frames, len(atom_labels)-0.5, -0.5])
    ax3.set_xlabel('Time Step')
    ax3.set_ylabel('Atom Index')
    ax3.set_title('RMSD Heatmap (All Atoms)')
    cbar = plt.colorbar(im, ax=ax3)
    cbar.set_label('RMSD (Å)', rotation=270, labelpad=20)
    
    # Plot 4: Distribution of RMSD at final time step (bottom right)
    ax4 = axes[1, 1]
    final_rmsd = rmsd_per_atom[-1, :]
    ax4.hist(final_rmsd, bins=30, color='steelblue', alpha=0.7, edgecolor='black')
    ax4.axvline(np.mean(final_rmsd), color='red', linestyle='--', linewidth=2, 
                label=f'Mean: {np.mean(final_rmsd):.3f}')
    ax4.axvline(np.median(final_rmsd), color='orange', linestyle='--', linewidth=2,
                label=f'Median: {np.median(final_rmsd):.3f}')
    ax4.set_xlabel('RMSD (Å)')
    ax4.set_ylabel('Number of Atoms')
    ax4.set_title('Distribution of RMSD at Final Time Step')
    ax4.legend(loc='best')
    ax4.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"RMSD plot saved to {output_file}")
    plt.close()

def plot_msd(msd_per_atom, msd_avg, msd_std, output_file='msd_analysis.png'):
    """
    Create comprehensive MSD plots.
    
    Parameters:
    -----------
    msd_per_atom : numpy.ndarray
        MSD for each atom at each time step
    msd_avg : numpy.ndarray
        Average MSD across all atoms at each time step
    msd_std : numpy.ndarray
        Standard deviation of MSD at each time step
    output_file : str
        Output filename for the plot
    """
    print(f"Creating MSD plots...")
    
    n_frames, n_atoms = msd_per_atom.shape
    time_steps = np.arange(n_frames)
    
    # Set up the plotting style
    sns.set_style("whitegrid")
    sns.set_context("paper", font_scale=1.2)
    
    # Create figure with subplots
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('Mean Squared Displacement (MSD) Analysis', fontsize=16, fontweight='bold')
    
    # Plot 1: Individual atom MSD (top left)
    ax1 = axes[0, 0]
    # Plot a subset of atoms if there are too many (for readability)
    n_plot = min(10, n_atoms)
    colors = plt.cm.viridis(np.linspace(0, 1, n_plot))
    for i in range(n_plot):
        ax1.plot(time_steps, msd_per_atom[:, i], alpha=0.7, linewidth=1.5, 
                color=colors[i], label=f'Atom {i+1}')
    ax1.set_xlabel('Time Step')
    ax1.set_ylabel('MSD (Ų)')
    ax1.set_title(f'MSD for Individual Atoms (showing {n_plot}/{n_atoms})')
    if n_plot <= 10:
        ax1.legend(loc='best', fontsize=8)
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Average MSD with error bars (top right)
    ax2 = axes[0, 1]
    ax2.plot(time_steps, msd_avg, 'b-', linewidth=2, label='Mean MSD')
    ax2.fill_between(time_steps, msd_avg - msd_std, msd_avg + msd_std, 
                     alpha=0.3, color='blue', label='±1 Std Dev')
    ax2.set_xlabel('Time Step')
    ax2.set_ylabel('MSD (Ų)')
    ax2.set_title('Average MSD Over All Atoms')
    ax2.legend(loc='best')
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: Heatmap of MSD for all atoms (bottom left)
    ax3 = axes[1, 0]
    # Downsample if there are too many atoms for visualization
    if n_atoms > 50:
        step = n_atoms // 50
        msd_plot = msd_per_atom[:, ::step].T
        atom_labels = np.arange(0, n_atoms, step)
    else:
        msd_plot = msd_per_atom.T
        atom_labels = np.arange(n_atoms)
    
    im = ax3.imshow(msd_plot, aspect='auto', cmap='hot', interpolation='nearest',
                    extent=[0, n_frames, len(atom_labels)-0.5, -0.5])
    ax3.set_xlabel('Time Step')
    ax3.set_ylabel('Atom Index')
    ax3.set_title('MSD Heatmap (All Atoms)')
    cbar = plt.colorbar(im, ax=ax3)
    cbar.set_label('MSD (Ų)', rotation=270, labelpad=20)
    
    # Plot 4: Distribution of MSD at final time step (bottom right)
    ax4 = axes[1, 1]
    final_msd = msd_per_atom[-1, :]
    ax4.hist(final_msd, bins=30, color='steelblue', alpha=0.7, edgecolor='black')
    ax4.axvline(np.mean(final_msd), color='red', linestyle='--', linewidth=2, 
                label=f'Mean: {np.mean(final_msd):.3f}')
    ax4.axvline(np.median(final_msd), color='orange', linestyle='--', linewidth=2,
                label=f'Median: {np.median(final_msd):.3f}')
    ax4.set_xlabel('MSD (Ų)')
    ax4.set_ylabel('Number of Atoms')
    ax4.set_title('Distribution of MSD at Final Time Step')
    ax4.legend(loc='best')
    ax4.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"MSD plot saved to {output_file}")
    plt.close()

def save_results(rmsd_per_atom, rmsd_avg, rmsd_std, msd_per_atom, msd_avg, msd_std, 
                output_file='displacement_data.npz'):
    """
    Save RMSD and MSD results to a compressed numpy file for later analysis.
    
    Parameters:
    -----------
    rmsd_per_atom : numpy.ndarray
        RMSD for each atom at each time step
    rmsd_avg : numpy.ndarray
        Average RMSD across all atoms at each time step
    rmsd_std : numpy.ndarray
        Standard deviation of RMSD at each time step
    msd_per_atom : numpy.ndarray
        MSD for each atom at each time step
    msd_avg : numpy.ndarray
        Average MSD across all atoms at each time step
    msd_std : numpy.ndarray
        Standard deviation of MSD at each time step
    output_file : str
        Output filename for the numpy data
    """
    print(f"Saving results to {output_file}...")
    np.savez_compressed(output_file,
                       rmsd_per_atom=rmsd_per_atom,
                       rmsd_avg=rmsd_avg,
                       rmsd_std=rmsd_std,
                       msd_per_atom=msd_per_atom,
                       msd_avg=msd_avg,
                       msd_std=msd_std)
    print(f"Results saved successfully!")

def main():
    """Main analysis workflow."""
    print("="*80)
    print("MONTE CARLO TRAJECTORY ANALYSIS")
    print("Root Mean Squared Displacement (RMSD) Analysis")
    print("="*80)
    print()
    
    # Parse command line arguments
    if len(sys.argv) > 1:
        trajectory_file = sys.argv[1]
    else:
        trajectory_file = 'trajectory.gro'
    
    # Read trajectory
    trajectory = read_trajectory(trajectory_file)
    
    if len(trajectory) < 2:
        print("ERROR: Trajectory has fewer than 2 frames. Cannot compute RMSD.")
        sys.exit(1)
    
    # Compute RMSD
    rmsd_per_atom, rmsd_avg, rmsd_std = compute_rmsd(trajectory)
    
    # Compute MSD
    msd_per_atom, msd_avg, msd_std = compute_msd(trajectory)
    
    # Create plots
    plot_rmsd(rmsd_per_atom, rmsd_avg, rmsd_std)
    plot_msd(msd_per_atom, msd_avg, msd_std)
    
    # Save results
    save_results(rmsd_per_atom, rmsd_avg, rmsd_std, msd_per_atom, msd_avg, msd_std)
    
    # Summary statistics
    print()
    print("="*80)
    print("ANALYSIS SUMMARY")
    print("="*80)
    print(f"Number of frames: {len(trajectory)}")
    print(f"Number of atoms: {len(trajectory[0])}")
    print()
    print("RMSD Statistics (at final time step):")
    print(f"  Mean RMSD: {rmsd_avg[-1]:.6f} Å")
    print(f"  Std Dev:   {rmsd_std[-1]:.6f} Å")
    print(f"  Min RMSD:  {np.min(rmsd_per_atom[-1, :]):.6f} Å")
    print(f"  Max RMSD:  {np.max(rmsd_per_atom[-1, :]):.6f} Å")
    print()
    print("MSD Statistics (at final time step):")
    print(f"  Mean MSD:  {msd_avg[-1]:.6f} Ų")
    print(f"  Std Dev:   {msd_std[-1]:.6f} Ų")
    print(f"  Min MSD:   {np.min(msd_per_atom[-1, :]):.6f} Ų")
    print(f"  Max MSD:   {np.max(msd_per_atom[-1, :]):.6f} Ų")
    print()
    print("Output files generated:")
    print("  - rmsd_analysis.png: RMSD visualization plots")
    print("  - msd_analysis.png: MSD visualization plots")
    print("  - displacement_data.npz: Raw RMSD and MSD data for further analysis")
    print()
    print("Analysis complete!")
    print("="*80)

if __name__ == "__main__":
    main()

