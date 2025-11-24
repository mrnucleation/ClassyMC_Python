#!/usr/bin/env python3
"""
Visualization script for RL training progress.

Monitors and plots:
- Acceptance rate over time
- Loss curves
- Reward progression
- Energy statistics
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import pickle


class TrainingMonitor:
    """Monitor and visualize RL training progress."""
    
    def __init__(self, log_file="rl_training_log.pkl"):
        self.log_file = log_file
        self.history = {
            'iterations': [],
            'acceptance_rates': [],
            'losses': [],
            'rewards': [],
            'entropy': [],
            'energy_changes': []
        }
    
    def log_iteration(self, iteration, acceptance_rate, loss, reward, 
                     entropy, energy_change):
        """Log data from a training iteration."""
        self.history['iterations'].append(iteration)
        self.history['acceptance_rates'].append(acceptance_rate)
        self.history['losses'].append(loss)
        self.history['rewards'].append(reward)
        self.history['entropy'].append(entropy)
        self.history['energy_changes'].append(energy_change)
        
        # Save to disk
        self.save()
    
    def save(self):
        """Save history to disk."""
        with open(self.log_file, 'wb') as f:
            pickle.dump(self.history, f)
    
    def load(self):
        """Load history from disk."""
        if os.path.exists(self.log_file):
            with open(self.log_file, 'rb') as f:
                self.history = pickle.load(f)
            return True
        return False
    
    def plot_progress(self, save_path="rl_training_progress.png"):
        """Create comprehensive progress plots."""
        if len(self.history['iterations']) == 0:
            print("No training data to plot")
            return
        
        fig, axes = plt.subplots(2, 3, figsize=(15, 8))
        fig.suptitle('RL Training Progress', fontsize=16, fontweight='bold')
        
        iterations = self.history['iterations']
        
        # Plot 1: Acceptance Rate
        ax = axes[0, 0]
        ax.plot(iterations, self.history['acceptance_rates'], 'b-', linewidth=2)
        ax.set_xlabel('Iteration')
        ax.set_ylabel('Acceptance Rate')
        ax.set_title('Neural Network Move Acceptance Rate')
        ax.grid(True, alpha=0.3)
        ax.set_ylim([0, 1])
        
        # Add trend line
        if len(iterations) > 2:
            z = np.polyfit(iterations, self.history['acceptance_rates'], 1)
            p = np.poly1d(z)
            ax.plot(iterations, p(iterations), "r--", alpha=0.5, label='Trend')
            ax.legend()
        
        # Plot 2: Loss
        ax = axes[0, 1]
        ax.plot(iterations, self.history['losses'], 'r-', linewidth=2)
        ax.set_xlabel('Iteration')
        ax.set_ylabel('Loss')
        ax.set_title('Training Loss')
        ax.grid(True, alpha=0.3)
        
        # Plot 3: Reward
        ax = axes[0, 2]
        ax.plot(iterations, self.history['rewards'], 'g-', linewidth=2)
        ax.set_xlabel('Iteration')
        ax.set_ylabel('Average Reward')
        ax.set_title('Average Reward per Move')
        ax.grid(True, alpha=0.3)
        
        # Plot 4: Entropy
        ax = axes[1, 0]
        ax.plot(iterations, self.history['entropy'], 'm-', linewidth=2)
        ax.set_xlabel('Iteration')
        ax.set_ylabel('Entropy')
        ax.set_title('Policy Entropy (Exploration)')
        ax.grid(True, alpha=0.3)
        
        # Plot 5: Energy Changes
        ax = axes[1, 1]
        ax.plot(iterations, self.history['energy_changes'], 'c-', linewidth=2)
        ax.axhline(y=0, color='k', linestyle='--', alpha=0.3)
        ax.set_xlabel('Iteration')
        ax.set_ylabel('Energy Change')
        ax.set_title('System Energy Change per Episode')
        ax.grid(True, alpha=0.3)
        
        # Plot 6: Summary Statistics
        ax = axes[1, 2]
        ax.axis('off')
        
        # Calculate statistics
        final_acc = self.history['acceptance_rates'][-1]
        max_acc = max(self.history['acceptance_rates'])
        mean_acc = np.mean(self.history['acceptance_rates'])
        
        final_loss = self.history['losses'][-1]
        mean_loss = np.mean(self.history['losses'])
        
        final_reward = self.history['rewards'][-1]
        max_reward = max(self.history['rewards'])
        
        # Display statistics
        stats_text = f"""
        Training Statistics
        ══════════════════════
        
        Acceptance Rate:
          Current: {final_acc:.4f}
          Best:    {max_acc:.4f}
          Mean:    {mean_acc:.4f}
        
        Loss:
          Current: {final_loss:.6f}
          Mean:    {mean_loss:.6f}
        
        Reward:
          Current: {final_reward:.4f}
          Best:    {max_reward:.4f}
        
        Total Iterations: {len(iterations)}
        """
        
        ax.text(0.1, 0.5, stats_text, fontsize=10, verticalalignment='center',
                family='monospace')
        
        plt.tight_layout()
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"\nPlot saved to: {save_path}")
        
        return fig
    
    def print_summary(self):
        """Print summary statistics."""
        if len(self.history['iterations']) == 0:
            print("No training data available")
            return
        
        print("\n" + "="*60)
        print("  RL TRAINING SUMMARY")
        print("="*60)
        
        acc_rates = self.history['acceptance_rates']
        losses = self.history['losses']
        rewards = self.history['rewards']
        
        print(f"\nAcceptance Rate:")
        print(f"  Initial:  {acc_rates[0]:.4f} ({acc_rates[0]*100:.2f}%)")
        print(f"  Final:    {acc_rates[-1]:.4f} ({acc_rates[-1]*100:.2f}%)")
        print(f"  Best:     {max(acc_rates):.4f} ({max(acc_rates)*100:.2f}%)")
        print(f"  Mean:     {np.mean(acc_rates):.4f} ({np.mean(acc_rates)*100:.2f}%)")
        print(f"  Std Dev:  {np.std(acc_rates):.4f}")
        
        # Calculate improvement
        improvement = (acc_rates[-1] - acc_rates[0]) / acc_rates[0] * 100
        print(f"  Improvement: {improvement:+.1f}%")
        
        print(f"\nLoss:")
        print(f"  Initial:  {losses[0]:.6f}")
        print(f"  Final:    {losses[-1]:.6f}")
        print(f"  Min:      {min(losses):.6f}")
        
        print(f"\nReward:")
        print(f"  Initial:  {rewards[0]:.4f}")
        print(f"  Final:    {rewards[-1]:.4f}")
        print(f"  Max:      {max(rewards):.4f}")
        
        print(f"\nTraining Info:")
        print(f"  Total iterations: {len(self.history['iterations'])}")
        print(f"  Final entropy: {self.history['entropy'][-1]:.6f}")
        print(f"  Final energy change: {self.history['energy_changes'][-1]:.6f}")
        
        print("\n" + "="*60)


def main():
    """Main function to load and visualize training data."""
    import argparse
    
    parser = argparse.ArgumentParser(description='Visualize RL training progress')
    parser.add_argument('--log-file', default='rl_training_log.pkl',
                       help='Path to training log file')
    parser.add_argument('--output', default='rl_training_progress.png',
                       help='Output plot filename')
    parser.add_argument('--show', action='store_true',
                       help='Show plot interactively')
    
    args = parser.parse_args()
    
    # Create monitor and load data
    monitor = TrainingMonitor(args.log_file)
    
    if not monitor.load():
        print(f"No training log found at: {args.log_file}")
        print("\nTo create training logs, modify train_nn_with_rl.py to use TrainingMonitor")
        return
    
    # Print summary
    monitor.print_summary()
    
    # Create plots
    fig = monitor.plot_progress(args.output)
    
    # Show interactively if requested
    if args.show:
        plt.show()
    else:
        plt.close(fig)


if __name__ == "__main__":
    main()

