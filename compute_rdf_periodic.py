#!/usr/bin/env python3
"""
Compute Radial Distribution Function (RDF) from XSF trajectory file with periodic boundary conditions.

This script reads an XSF trajectory file and computes the RDF g(r) for all atom pairs,
taking into account periodic boundary conditions. The RDF is computed as:

g(r) = <n(r)> / (4πr²ρdr)

where:
- <n(r)> is the average number of atoms in a shell of thickness dr at distance r
- ρ is the number density of atoms
- The calculation uses minimum image convention for periodic boundaries

Usage:
    python compute_rdf_periodic.py trajectory.xsf [options]

Output:
    - rdf_periodic_data.txt: RDF data file
    - rdf_periodic_plot.png: RDF visualization
"""

import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys
import os
from typing import List, Tuple, Optional
import time

class XSFParser:
    """Parser for XSF (XCrySDen Structure File) format."""
    
    def __init__(self, filename: str):
        self.filename = filename
        self.frames = []
        self.box_vectors = []
        self.n_atoms = 0
        self.atom_types = []
        
    def parse(self) -> bool:
        """Parse the XSF file and extract trajectory data."""
        try:
            with open(self.filename, 'r') as f:
                lines = f.readlines()
            
            i = 0
            frame_count = 0
            
            while i < len(lines):
                line = lines[i].strip()
                
                if line == "CRYSTAL":
                    # Parse crystal structure
                    frame_data = self._parse_crystal_frame(lines, i)
                    if frame_data is not None:
                        self.frames.append(frame_data)
                        frame_count += 1
                    i = frame_data['next_line'] if frame_data else i + 1
                else:
                    i += 1
            
            print(f"Parsed {frame_count} frames from {self.filename}")
            if frame_count > 0:
                self.n_atoms = len(self.frames[0]['atoms'])
                self.atom_types = list(set([atom[0] for atom in self.frames[0]['atoms']]))
                print(f"Number of atoms per frame: {self.n_atoms}")
                print(f"Atom types found: {self.atom_types}")
            
            return frame_count > 0
            
        except Exception as e:
            print(f"Error parsing XSF file: {e}")
            return False
    
    def _parse_crystal_frame(self, lines: List[str], start_idx: int) -> Optional[dict]:
        """Parse a single crystal frame from XSF format."""
        try:
            i = start_idx + 1
            
            # Parse PRIMVEC
            if i >= len(lines) or lines[i].strip() != "PRIMVEC":
                return None
            
            i += 1
            box_vectors = []
            for j in range(3):
                if i + j >= len(lines):
                    return None
                vec = [float(x) for x in lines[i + j].strip().split()]
                if len(vec) != 3:
                    return None
                box_vectors.append(vec)
            
            i += 3
            
            # Parse PRIMCOORD
            if i >= len(lines) or lines[i].strip() != "PRIMCOORD":
                return None
            
            i += 1
            if i >= len(lines):
                return None
            
            # Parse number of atoms and number of atom types
            header = lines[i].strip().split()
            if len(header) < 2:
                return None
            
            n_atoms = int(header[0])
            n_types = int(header[1])
            
            i += 1
            
            # Parse atomic coordinates
            atoms = []
            for j in range(n_atoms):
                if i + j >= len(lines):
                    return None
                atom_line = lines[i + j].strip().split()
                if len(atom_line) < 4:
                    return None
                
                atom_type = int(atom_line[0])
                coords = [float(atom_line[k]) for k in range(1, 4)]
                atoms.append((atom_type, coords))
            
            i += n_atoms
            
            return {
                'box_vectors': np.array(box_vectors),
                'atoms': atoms,
                'next_line': i
            }
            
        except Exception as e:
            print(f"Error parsing crystal frame: {e}")
            return None

class RDFCalculator:
    """Calculate Radial Distribution Function with periodic boundary conditions."""
    
    def __init__(self, r_max: float = 5.0, dr: float = 0.05, n_bins: Optional[int] = None):
        self.r_max = r_max
        self.dr = dr
        self.n_bins = n_bins or int(r_max / dr)
        self.r_max_sq = r_max**2
        
        # Create bin edges
        self.bin_edges = np.linspace(0, r_max, self.n_bins + 1)
        self.bin_centers = (self.bin_edges[:-1] + self.bin_edges[1:]) / 2
        self.bin_widths = np.diff(self.bin_edges)
        
        # Histogram array
        self.histogram = np.zeros(self.n_bins)
        self.n_samples = 0
        
    def compute_rdf(self, frames: List[dict], atom_types: Optional[List[int]] = None) -> Tuple[np.ndarray, np.ndarray]:
        """
        Compute RDF from trajectory frames.
        
        Args:
            frames: List of frame data from XSF parser
            atom_types: List of atom types to include (None for all types)
            
        Returns:
            Tuple of (r_values, g_r_values)
        """
        print(f"Computing RDF from {len(frames)} frames...")
        print(f"Parameters: r_max={self.r_max}, dr={self.dr}, n_bins={self.n_bins}")
        
        start_time = time.time()
        
        for frame_idx, frame in enumerate(frames):
            if frame_idx % 100 == 0:
                print(f"Processing frame {frame_idx + 1}/{len(frames)}")
            
            self._compute_frame_rdf(frame, atom_types)
            self.n_samples += 1
        
        # Normalize histogram to get g(r)
        g_r = self._normalize_histogram(frames[0], atom_types)
        
        end_time = time.time()
        print(f"RDF computation completed in {end_time - start_time:.2f} seconds")
        
        return self.bin_centers, g_r
    
    def _compute_frame_rdf(self, frame: dict, atom_types: Optional[List[int]] = None):
        """Compute RDF contribution from a single frame."""
        atoms = frame['atoms']
        box_vectors = frame['box_vectors']
        
        # Filter atoms by type if specified
        if atom_types is not None:
            atoms = [atom for atom in atoms if atom[0] in atom_types]
        
        n_atoms = len(atoms)
        
        # Compute all pairwise distances with periodic boundary conditions
        for i in range(n_atoms):
            for j in range(i + 1, n_atoms):
                # Get coordinates
                pos_i = np.array(atoms[i][1])
                pos_j = np.array(atoms[j][1])
                
                # Compute distance with minimum image convention
                dr = pos_i - pos_j
                dr = self._apply_periodic_boundary(dr, box_vectors)
                
                r_sq = np.sum(dr**2)
                
                if r_sq < self.r_max_sq:
                    r = np.sqrt(r_sq)
                    bin_idx = int(r / self.dr)
                    if bin_idx < self.n_bins:
                        self.histogram[bin_idx] += 2.0  # Count both i-j and j-i pairs
    
    def _apply_periodic_boundary(self, dr: np.ndarray, box_vectors: np.ndarray) -> np.ndarray:
        """Apply periodic boundary conditions using minimum image convention."""
        # For cubic boxes, this is straightforward
        # For general boxes, we'd need more complex transformations
        
        # Assume cubic box for simplicity (can be extended for general boxes)
        box_lengths = np.diag(box_vectors)  # Extract diagonal elements
        
        # Apply minimum image convention
        for dim in range(3):
            if box_lengths[dim] > 0:
                half_box = box_lengths[dim] / 2.0
                while dr[dim] > half_box:
                    dr[dim] -= box_lengths[dim]
                while dr[dim] < -half_box:
                    dr[dim] += box_lengths[dim]
        
        return dr
    
    def _normalize_histogram(self, frame: dict, atom_types: Optional[List[int]] = None) -> np.ndarray:
        """Normalize histogram to get g(r)."""
        atoms = frame['atoms']
        box_vectors = frame['box_vectors']
        
        # Filter atoms by type if specified
        if atom_types is not None:
            atoms = [atom for atom in atoms if atom[0] in atom_types]
        
        n_atoms = len(atoms)
        
        # Calculate box volume
        box_volume = np.abs(np.linalg.det(box_vectors))
        density = n_atoms / box_volume
        
        # Calculate normalization factors
        g_r = np.zeros(self.n_bins)
        
        for i in range(self.n_bins):
            r = self.bin_centers[i]
            dr = self.bin_widths[i]
            
            # Volume of spherical shell: 4πr²dr
            shell_volume = 4.0 * np.pi * r**2 * dr
            
            # Ideal gas number of pairs in shell
            ideal_pairs = density * shell_volume
            
            # Normalize histogram
            if ideal_pairs > 0 and self.n_samples > 0:
                g_r[i] = self.histogram[i] / (self.n_samples * ideal_pairs * n_atoms)
            else:
                g_r[i] = 0.0
        
        return g_r

def save_rdf_data(r_values: np.ndarray, g_r_values: np.ndarray, filename: str):
    """Save RDF data to file."""
    with open(filename, 'w') as f:
        f.write("# r\tg(r)\n")
        for r, g_r in zip(r_values, g_r_values):
            f.write(f"{r:.6f}\t{g_r:.6f}\n")
    print(f"RDF data saved to {filename}")

def plot_rdf(r_values: np.ndarray, g_r_values: np.ndarray, filename: str):
    """Create RDF plot."""
    plt.figure(figsize=(10, 6))
    plt.plot(r_values, g_r_values, 'b-', linewidth=2, label='g(r)')
    plt.axhline(y=1.0, color='r', linestyle='--', alpha=0.7, label='Ideal gas (g(r)=1)')
    
    plt.xlabel('Distance r', fontsize=12)
    plt.ylabel('Radial Distribution Function g(r)', fontsize=12)
    plt.title('Radial Distribution Function with Periodic Boundary Conditions', fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.legend()
    
    # Set reasonable axis limits
    plt.xlim(0, np.max(r_values))
    plt.ylim(0, max(1.5, np.max(g_r_values) * 1.1))
    
    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"RDF plot saved to {filename}")

def main():
    """Main function."""
    parser = argparse.ArgumentParser(description='Compute RDF from XSF trajectory with periodic boundary conditions')
    parser.add_argument('trajectory_file', help='XSF trajectory file')
    parser.add_argument('--r-max', type=float, default=5.0, help='Maximum distance for RDF (default: 5.0)')
    parser.add_argument('--dr', type=float, default=0.05, help='Bin width for RDF (default: 0.05)')
    parser.add_argument('--n-bins', type=int, help='Number of bins (overrides dr if specified)')
    parser.add_argument('--atom-types', nargs='+', type=int, help='Atom types to include (default: all types)')
    parser.add_argument('--output-prefix', default='rdf_periodic', help='Output file prefix (default: rdf_periodic)')
    parser.add_argument('--max-frames', type=int, help='Maximum number of frames to process')
    
    args = parser.parse_args()
    
    # Check if input file exists
    if not os.path.exists(args.trajectory_file):
        print(f"Error: Input file '{args.trajectory_file}' not found")
        sys.exit(1)
    
    print("="*60)
    print("RDF COMPUTATION WITH PERIODIC BOUNDARY CONDITIONS")
    print("="*60)
    print(f"Input file: {args.trajectory_file}")
    print(f"Parameters:")
    print(f"  r_max: {args.r_max}")
    print(f"  dr: {args.dr}")
    if args.n_bins:
        print(f"  n_bins: {args.n_bins}")
    if args.atom_types:
        print(f"  atom_types: {args.atom_types}")
    if args.max_frames:
        print(f"  max_frames: {args.max_frames}")
    print()
    
    # Parse XSF file
    parser = XSFParser(args.trajectory_file)
    if not parser.parse():
        print("Error: Failed to parse XSF file")
        sys.exit(1)
    
    # Limit frames if specified
    frames = parser.frames
    if args.max_frames and len(frames) > args.max_frames:
        frames = frames[:args.max_frames]
        print(f"Limited to first {args.max_frames} frames")
    
    # Create RDF calculator
    rdf_calc = RDFCalculator(
        r_max=args.r_max,
        dr=args.dr,
        n_bins=args.n_bins
    )
    
    # Compute RDF
    r_values, g_r_values = rdf_calc.compute_rdf(frames, args.atom_types)
    
    # Save results
    data_filename = f"{args.output_prefix}_data.txt"
    plot_filename = f"{args.output_prefix}_plot.png"
    
    save_rdf_data(r_values, g_r_values, data_filename)
    plot_rdf(r_values, g_r_values, plot_filename)
    
    # Print summary statistics
    print("\nRDF SUMMARY:")
    print(f"  Number of frames processed: {len(frames)}")
    print(f"  Number of atoms per frame: {parser.n_atoms}")
    print(f"  Atom types: {parser.atom_types}")
    print(f"  RDF range: 0 to {args.r_max}")
    print(f"  Number of bins: {len(r_values)}")
    print(f"  Maximum g(r): {np.max(g_r_values):.3f}")
    print(f"  Output files: {data_filename}, {plot_filename}")
    
    print("\nRDF computation completed successfully!")

if __name__ == "__main__":
    main()