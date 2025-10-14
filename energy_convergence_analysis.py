#!/usr/bin/env python3
"""
Energy Convergence Analysis and Time Lapse Visualization
Analyzes LAMMPS energy data to show convergence of running averages to the total average.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.patches import Rectangle
import pandas as pd
import os

class EnergyConvergenceAnalyzer:
    def __init__(self, data_file="avg_energy.txt"):
        """Initialize the analyzer with energy data file."""
        self.data_file = data_file
        self.timesteps = []
        self.energies = []
        self.total_average = None
        self.running_averages = []
        self.convergence_threshold = None
        
    def load_data(self):
        """Load energy data from file."""
        if not os.path.exists(self.data_file):
            raise FileNotFoundError(f"Data file '{self.data_file}' not found.")
        
        print(f"Loading data from {self.data_file}...")
        
        with open(self.data_file, 'r') as file:
            for line_num, line in enumerate(file, 1):
                line = line.strip()
                
                # Skip empty lines and comment lines
                if not line or line.startswith('#'):
                    continue
                
                # Split the line and get timestep and energy
                parts = line.split()
                if len(parts) < 2:
                    continue
                
                try:
                    timestep = int(parts[0])
                    energy = float(parts[1])
                    self.timesteps.append(timestep)
                    self.energies.append(energy)
                except ValueError:
                    continue
        
        self.timesteps = np.array(self.timesteps)
        self.energies = np.array(self.energies)
        
        # Calculate total average
        self.total_average = np.mean(self.energies)
        
        # Calculate running averages
        self.running_averages = np.cumsum(self.energies) / np.arange(1, len(self.energies) + 1)
        
        print(f"Loaded {len(self.timesteps)} data points")
        print(f"Total average energy: {self.total_average:.6f}")
        
    def calculate_convergence_metrics(self, window_size=1000):
        """Calculate convergence metrics."""
        print("Calculating convergence metrics...")
        
        # Calculate rolling standard deviation
        rolling_std = pd.Series(self.running_averages).rolling(window=window_size).std()
        
        # Find convergence point (when rolling std becomes small)
        convergence_threshold = 0.01 * abs(self.total_average)  # 1% of total average
        convergence_indices = np.where(rolling_std < convergence_threshold)[0]
        
        if len(convergence_indices) > 0:
            self.convergence_threshold = convergence_indices[0]
            convergence_timestep = self.timesteps[self.convergence_threshold]
            print(f"Convergence achieved at timestep: {convergence_timestep}")
            print(f"Convergence threshold: {convergence_threshold:.6f}")
        else:
            print("No clear convergence point found with current threshold")
            
        return rolling_std
    
    def create_static_convergence_plot(self):
        """Create a static plot showing energy convergence."""
        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 10))
        
        # Plot 1: Raw energy vs time
        ax1.plot(self.timesteps, self.energies, 'b-', alpha=0.7, linewidth=0.8)
        ax1.axhline(y=self.total_average, color='r', linestyle='--', linewidth=2, 
                   label=f'Total Average: {self.total_average:.6f}')
        ax1.set_xlabel('Timestep')
        ax1.set_ylabel('Energy')
        ax1.set_title('Raw Energy vs Time')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # Plot 2: Running average vs time
        ax2.plot(self.timesteps, self.running_averages, 'g-', linewidth=2, 
                label='Running Average')
        ax2.axhline(y=self.total_average, color='r', linestyle='--', linewidth=2,
                   label=f'Total Average: {self.total_average:.6f}')
        
        if self.convergence_threshold is not None:
            ax2.axvline(x=self.timesteps[self.convergence_threshold], color='orange', 
                       linestyle=':', linewidth=2, label='Convergence Point')
        
        ax2.set_xlabel('Timestep')
        ax2.set_ylabel('Running Average Energy')
        ax2.set_title('Running Average Convergence')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        # Plot 3: Convergence error
        convergence_error = np.abs(self.running_averages - self.total_average)
        ax3.semilogy(self.timesteps, convergence_error, 'purple', linewidth=2)
        ax3.set_xlabel('Timestep')
        ax3.set_ylabel('|Running Avg - Total Avg|')
        ax3.set_title('Convergence Error (Log Scale)')
        ax3.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig('energy_convergence_static.png', dpi=300, bbox_inches='tight')
        plt.show()
        
    def create_animated_convergence(self, save_gif=False, frame_interval=50):
        """Create an animated time lapse of energy convergence."""
        print("Creating animated convergence visualization...")
        
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8))
        
        # Initialize empty lines
        line1, = ax1.plot([], [], 'b-', alpha=0.7, linewidth=0.8)
        line2, = ax2.plot([], [], 'g-', linewidth=2)
        
        # Set up the plot
        ax1.set_xlim(self.timesteps[0], self.timesteps[-1])
        ax1.set_ylim(np.min(self.energies) * 1.1, np.max(self.energies) * 1.1)
        ax1.set_xlabel('Timestep')
        ax1.set_ylabel('Energy')
        ax1.set_title('Raw Energy vs Time')
        ax1.grid(True, alpha=0.3)
        
        ax2.set_xlim(self.timesteps[0], self.timesteps[-1])
        ax2.set_ylim(np.min(self.running_averages) * 1.1, np.max(self.running_averages) * 1.1)
        ax2.set_xlabel('Timestep')
        ax2.set_ylabel('Running Average Energy')
        ax2.set_title('Running Average Convergence')
        ax2.grid(True, alpha=0.3)
        
        # Add horizontal lines for total average
        ax1.axhline(y=self.total_average, color='r', linestyle='--', linewidth=2,
                   label=f'Total Average: {self.total_average:.6f}')
        ax2.axhline(y=self.total_average, color='r', linestyle='--', linewidth=2,
                   label=f'Total Average: {self.total_average:.6f}')
        
        ax1.legend()
        ax2.legend()
        
        # Add text for current statistics
        stats_text = ax2.text(0.02, 0.98, '', transform=ax2.transAxes, 
                            verticalalignment='top', fontsize=10,
                            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
        
        def animate(frame):
            # Update data up to current frame
            current_timesteps = self.timesteps[:frame+1]
            current_energies = self.energies[:frame+1]
            current_running_avg = self.running_averages[:frame+1]
            
            # Update lines
            line1.set_data(current_timesteps, current_energies)
            line2.set_data(current_timesteps, current_running_avg)
            
            # Update statistics text
            if frame > 0:
                current_avg = current_running_avg[-1]
                error = abs(current_avg - self.total_average)
                stats_text.set_text(f'Current Average: {current_avg:.6f}\n'
                                  f'Error: {error:.6f}\n'
                                  f'Timestep: {current_timesteps[-1]}')
            
            return line1, line2, stats_text
        
        # Create animation
        anim = animation.FuncAnimation(fig, animate, frames=len(self.timesteps),
                                     interval=frame_interval, blit=False, repeat=True)
        
        if save_gif:
            print("Saving animated GIF...")
            anim.save('energy_convergence_animation.gif', writer='pillow', fps=20)
            print("Animation saved as 'energy_convergence_animation.gif'")
        
        plt.show()
        return anim
    
    def create_energy_histogram(self):
        """Create a histogram of energy states."""
        print("Creating energy histogram...")
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
        
        # Histogram of raw energies
        ax1.hist(self.energies, bins=50, alpha=0.7, color='blue', edgecolor='black')
        ax1.axvline(self.total_average, color='red', linestyle='--', linewidth=2, 
                   label=f'Mean: {self.total_average:.6f}')
        ax1.set_xlabel('Energy')
        ax1.set_ylabel('Frequency')
        ax1.set_title('Histogram of Energy States')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # Histogram of running averages (to show convergence)
        ax2.hist(self.running_averages, bins=50, alpha=0.7, color='green', edgecolor='black')
        ax2.axvline(self.total_average, color='red', linestyle='--', linewidth=2,
                   label=f'Mean: {self.total_average:.6f}')
        ax2.set_xlabel('Running Average Energy')
        ax2.set_ylabel('Frequency')
        ax2.set_title('Histogram of Running Averages')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig('energy_histogram.png', dpi=300, bbox_inches='tight')
        plt.show()
        
        # Print histogram statistics
        print(f"\nEnergy Histogram Statistics:")
        print(f"Mean energy: {self.total_average:.6f}")
        print(f"Standard deviation: {np.std(self.energies):.6f}")
        print(f"Min energy: {np.min(self.energies):.6f}")
        print(f"Max energy: {np.max(self.energies):.6f}")
        print(f"Energy range: {np.max(self.energies) - np.min(self.energies):.6f}")
        
        # Calculate percentiles
        percentiles = [5, 25, 50, 75, 95]
        energy_percentiles = np.percentile(self.energies, percentiles)
        print(f"\nEnergy Percentiles:")
        for p, val in zip(percentiles, energy_percentiles):
            print(f"{p}th percentile: {val:.6f}")
    
    def create_convergence_analysis_report(self):
        """Create a detailed convergence analysis report."""
        print("\n" + "="*60)
        print("ENERGY CONVERGENCE ANALYSIS REPORT")
        print("="*60)
        
        print(f"Total number of data points: {len(self.timesteps)}")
        print(f"Simulation duration: {self.timesteps[0]} to {self.timesteps[-1]} timesteps")
        print(f"Total average energy: {self.total_average:.8f}")
        print(f"Energy standard deviation: {np.std(self.energies):.8f}")
        print(f"Energy range: {np.min(self.energies):.8f} to {np.max(self.energies):.8f}")
        
        # Calculate convergence metrics
        final_error = abs(self.running_averages[-1] - self.total_average)
        relative_error = final_error / abs(self.total_average) * 100
        
        print(f"\nConvergence Metrics:")
        print(f"Final running average: {self.running_averages[-1]:.8f}")
        print(f"Final absolute error: {final_error:.8f}")
        print(f"Final relative error: {relative_error:.4f}%")
        
        # Find when convergence is achieved within 1% and 0.1%
        error_1_percent = abs(self.total_average) * 0.01
        error_0_1_percent = abs(self.total_average) * 0.001
        
        convergence_1_percent = np.where(np.abs(self.running_averages - self.total_average) < error_1_percent)[0]
        convergence_0_1_percent = np.where(np.abs(self.running_averages - self.total_average) < error_0_1_percent)[0]
        
        if len(convergence_1_percent) > 0:
            timestep_1_percent = self.timesteps[convergence_1_percent[0]]
            print(f"Convergence within 1% achieved at timestep: {timestep_1_percent}")
        
        if len(convergence_0_1_percent) > 0:
            timestep_0_1_percent = self.timesteps[convergence_0_1_percent[0]]
            print(f"Convergence within 0.1% achieved at timestep: {timestep_0_1_percent}")
        
        print("="*60)

def main():
    """Main function to run the energy convergence analysis."""
    analyzer = EnergyConvergenceAnalyzer()
    
    try:
        # Load data
        analyzer.load_data()
        
        # Calculate convergence metrics
        rolling_std = analyzer.calculate_convergence_metrics()
        
        # Create static plot
        analyzer.create_static_convergence_plot()
        
        # Create animated visualization (without GIF)
        #anim = analyzer.create_animated_convergence(save_gif=False)
        
        # Create energy histogram
        analyzer.create_energy_histogram()
        
        # Generate analysis report
        analyzer.create_convergence_analysis_report()
        
        print("\nAnalysis complete! Generated files:")
        print("- energy_convergence_static.png (static plot)")
        print("- energy_histogram.png (energy states histogram)")
        print("- Interactive plots displayed (no GIF generated)")
        
    except Exception as e:
        print(f"Error during analysis: {e}")
        return 1
    
    return 0

if __name__ == "__main__":
    exit(main())
