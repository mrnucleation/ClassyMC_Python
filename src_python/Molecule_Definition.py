#Graph library
import networkx as nx
import json
import sys
import os
import numpy as np
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
#Molecule class
class Molecule_Topo:
    #--------------------------------
    def __init__(self, json_data):
        self.name = json_data['name']
        self.atoms = json_data['atoms']
        self.atomtypes = json_data['atomtypes']
        if 'bonds' in json_data and len(self.atoms) > 1:
            self.bonds = json_data['bonds']
            # Ensure bonds are tuples of atom indices
            for i in range(len(self.bonds)):
                if isinstance(self.bonds[i], list):
                    self.bonds[i] = tuple(self.bonds[i])
        self.angles = []
        self.torsion = [] 
        
           
            
        self.create_graph()
        self.auto_generate_angles_and_dihedrals()
        assert self.valid_bonds() == True, "Bonds are not valid"
        assert self.valid_angles() == True, "Angles are not valid"
        assert self.valid_torsional() == True, "Torsional are not valid"
   
    #--------------------------------
    def valid_bonds(self):
        '''
         For a bond to defined as valid it must have the following
            1. The bond must have which two atoms are part of the member, no more no less
            2. The atom IDs must match the total number of atoms specified
            3. The bond must have a type associated with it to specify a forcefield.
            
            Bond entry format:
                [(atom1, atom2)] or [(atom1, atom2, bond_type, bond_length)]
            Bond exit format:
                [(atom1, atom2, bond_type=None, bond_length=0.0)]
                
            This is because the bond_type will be determined by the forcefield.
        '''
        
        if len(self.bonds) == 0:
            return True
        
        try:
            for bond in self.bonds:
                assert len(bond) == 2, "Bond must have two atoms"
                assert bond[0][0] < len(self.atoms), "Atom 1 is not a valid atom"
                assert bond[0][1] < len(self.atoms), "Atom 2 is not a valid atom"
                assert bond[0][0] != bond[0][1], "Atom 1 and atom 2 are the same"
                assert bond[0][0] < bond[0][1], "Atom 1 is greater than atom 2"
                assert bond[0][0] >= 0
        except AssertionError:
            return False
        
        for bond in self.bonds:
            if len(bond) == 1:
                bond.append(None)
        return True
        
    #--------------------------------
    def valid_angles(self):
        '''
        For a angle to defined as valid it must have the following
            1. The angle must have which three atoms are part of the member, no more no less
            2. The atom IDs must match the total number of atoms specified
            3. The angle must have a type associated with it to specify a forcefield.
            
            Angle entry format:
                [(atom1, atom2, atom3)] or [(atom1, atom2, atom3, angle_type, angle_value)]
            Angle exit format:
                [(atom1, atom2, atom3, angle_type=None, angle_value=0.0)]
        '''
        if len(self.angles) == 0:
            return True
        try:
            for angle in self.angles:
                assert len(angle) == 3
                assert angle[0] < len(self.atoms)
                assert angle[1] < len(self.atoms)
                assert angle[2] < len(self.atoms)
        except AssertionError:
            return False
        
        for angle in self.angles:
            if len(angle) == 1:
                angle.append(None)
        return True
    
    #--------------------------------
    def valid_torsional(self):
        '''
        For a torsional to defined as valid it must have the following
            1. The torsional must have which four atoms are part of the member, no more no less
            2. The atom IDs must match the total number of atoms specified
            3. The torsional must have a type associated with it to specify a forcefield.
            
            Torsional entry format:
                [(atom1, atom2, atom3, atom4)] or [(atom1, atom2, atom3, atom4, torsional_type)]
            Torsional exit format:
                [(atom1, atom2, atom3, atom4, torsional_type=None)]
        '''
        if len(self.torsion) == 0:
            return True
        try:
            for torsional in self.torsion:
                assert len(torsional) == 4
                assert torsional[0] < len(self.atoms)
                assert torsional[1] < len(self.atoms)
                assert torsional[2] < len(self.atoms)
                assert torsional[3] < len(self.atoms)
        except AssertionError:
            return False
        
        for torsional in self.torsion:
            if len(torsional) == 1:
                torsional.append(None)
        return True
    #--------------------------------
    def create_graph(self):
        '''
         Creates a graph of the molecule. This is used in the CBMC and related algorithms to determine 
         regrowth order. 
        '''
        G = nx.Graph()
        for atom in self.atoms:
            G.add_node(atom[0])
        for bond in self.bonds:
            G.add_edge(bond[0][0], bond[0][1])
        self.graph = G
        
        #Check if the graph is connected
        assert nx.is_connected(G), "Molecule has disonnected atoms, all atoms must be part of a bond"
 
    #--------------------------------
    def find_angles(self):
        '''
        Find all angles in the molecule graph.
        An angle is defined by three atoms (i, j, k) where i-j and j-k are bonded.
        Returns list of angles in format [(atom1, atom2, atom3)] where atom2 is the central atom.
        '''
        angles = []
        G = self.graph
        
        # Find all paths of length 2 (angles)
        for node in G.nodes():
            # Get neighbors of the central atom
            neighbors = list(G.neighbors(node))
            # For each pair of neighbors, create an angle
            for i in range(len(neighbors)):
                for j in range(i + 1, len(neighbors)):
                    # Sort the outer atoms to ensure consistent ordering
                    angle_atoms = sorted([neighbors[i], neighbors[j]])
                    angles.append([angle_atoms[0], node, angle_atoms[1]])
        
        return angles
    
    #--------------------------------
    def find_proper_dihedrals(self):
        '''
        Find all proper dihedral angles in the molecule graph.
        A proper dihedral is defined by four atoms (i, j, k, l) where i-j, j-k, and k-l are bonded.
        Returns list of dihedrals in format [(atom1, atom2, atom3, atom4)].
        '''
        dihedrals = []
        G = self.graph
        
        # Find all paths of length 3 (dihedrals)
        for edge in G.edges():
            # For each bond (j-k), find paths of length 2 from each end
            j, k = edge
            
            # Get neighbors of j (excluding k)
            j_neighbors = [n for n in G.neighbors(j) if n != k]
            # Get neighbors of k (excluding j)
            k_neighbors = [n for n in G.neighbors(k) if n != j]
            
            # Create dihedrals: i-j-k-l where i is neighbor of j, l is neighbor of k
            for i in j_neighbors:
                for l in k_neighbors:
                    # Ensure i and l are different atoms
                    if i != l:
                        dihedrals.append([i, j, k, l])
        
        return dihedrals
    
    #--------------------------------
    def auto_generate_angles_and_dihedrals(self):
        '''
        Automatically generate angles and proper dihedral angles from the graph
        if they are not provided in the JSON data.
        '''
        # Only generate if not already provided
        if len(self.angles) == 0:
            self.angles = self.find_angles()
        
        if len(self.torsion) == 0:
            self.torsion = self.find_proper_dihedrals()
    #--------------------------------
    def save_tojson(self, filename):
        '''
        Saves the molecule to a json file.
        '''
        with open(filename, 'w') as f:
            json.dump(self.__dict__, f, indent=4)
            
            
    #--------------------------------
#===============================================
class Molecule_Type:
    """
     A larger scale
    """
    def __init__(self, topo_file, atomtypes, bondtype=None, angletype=None, torsiontype=None):
        self.topo = Molecule_Topo(topo_file)

        self.n_atoms = len(self.topo.atoms)
        self.n_bonds = len(self.topo.bonds)
        self.n_angles = len(self.topo.angles)
        self.n_torsions = len(self.topo.torsion)
        
        # Set default 
        self.atomtypes = atomtypes
        
        assert len(self.atomtypes) == self.n_atoms, "Atom types must match number of atoms"

        self.bondtype = bondtype if bondtype is not None else [None] * self.n_bonds
        self.angletype = angletype if angletype is not None else [None] * self.n_angles
        self.torsiontype = torsiontype if torsiontype is not None else [None] * self.n_torsions
        
        
    #--------------------------------
    def set_bond_type(self, bond_id, bond_type):
        '''
        Sets the bond type for a given bond id.
        '''
        self.topo.bond_ff_type(bond_id, bond_type)
        self.bondtype[bond_id] = bond_type
    #--------------------------------
    def set_angle_type(self, angle_id, angle_type):
        '''
        Sets the angle type for a given angle id.
        '''
        self.topo.angle_ff_type(angle_id, angle_type)
        self.angletype[angle_id] = angle_type
    #--------------------------------
    def set_torsion_type(self, torsion_id, torsion_type):
        '''
        Sets the torsion type for a given torsion id.
        '''
        self.topo.torsion_ff_type(torsion_id, torsion_type)
        self.torsiontype[torsion_id] = torsion_type
        
#===============================================
class Molecule:
    def __init__(self, 
                 moltype_object:Molecule_Type, 
                 position_array:np.ndarray):
        self.moltype = moltype_object
        self.position_array = position_array
        
        assert len(self.position_array) == len(self.moltype.topo.atoms), "Position array must be the same length as the number of atoms"
        assert len(self.position_array.shape) == 2, "Position array must be a 2D array"
        assert self.position_array.shape[0] == len(self.moltype.topo.atoms), "Position array must be the same length as the number of atoms"
        
        # Additional molecule properties
        self.mol_id = None  # Will be set when added to box
        self.box_id = None  # Will be set when added to box
        self.energy = 0.0   # Total energy of molecule
        self.intra_energy = 0.0  # Intramolecular energy
        self.inter_energy = 0.0  # Intermolecular energy
        
    #--------------------------------
    def get_moltype(self):
        return self.moltype
    
    #--------------------------------
    def get_positions(self):
        """Get the current positions of all atoms in the molecule"""
        return self.position_array
    
    #--------------------------------
    def set_positions(self, new_positions):
        """Set new positions for all atoms in the molecule"""
        assert len(new_positions) == len(self.position_array), "New positions must match atom count"
        self.position_array = new_positions
    
    #--------------------------------
    def get_atom_count(self):
        """Get the number of atoms in this molecule"""
        return len(self.position_array)
    
    #--------------------------------
    def get_bond_count(self):
        """Get the number of bonds in this molecule"""
        return self.moltype.n_bonds
    
    #--------------------------------
    def get_angle_count(self):
        """Get the number of angles in this molecule"""
        return self.moltype.n_angles
    
    #--------------------------------
    def get_torsion_count(self):
        """Get the number of torsions in this molecule"""
        return self.moltype.n_torsions
    
    #--------------------------------
    def compute_center_of_mass(self):
        """Compute the center of mass of the molecule"""
        if len(self.position_array) == 0:
            return np.zeros(3)
        
        # For now, assume equal mass for all atoms
        # In a real implementation, this would use actual atomic masses
        return np.mean(self.position_array, axis=0)
    
    #--------------------------------
    def to_dict(self):
        """Convert molecule to dictionary representation"""
        return {
            'moltype': self.moltype.topo.name,
            'positions': self.position_array.tolist(),
            'mol_id': self.mol_id,
            'box_id': self.box_id,
            'energy': self.energy
        }
    
    #--------------------------------
    def __str__(self):
        """String representation of the molecule"""
        return f"Molecule(type={self.moltype.topo.name}, atoms={len(self.position_array)}, id={self.mol_id})"
    
    #--------------------------------
    def __repr__(self):
        return self.__str__()
