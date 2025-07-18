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
    def __init__(self, json_file):
        with open(json_file, 'r') as f:
            data = json.load(f)
        self.name = data['name']
        self.atoms = data['atoms']
        self.atomtypes = data['atomtypes']
        if 'bonds' in data:
            self.bonds = data['bonds']
        else:
            self.bonds = []
        if 'angles' in data:
            self.angles = data['angles']
        else:
            self.angles = []
        if 'torsion' in data:
            self.torsion = data['torsion']
        else:
            self.torsion = []
        
        
           
            
        assert self.valid_bonds() == True, "Bonds are not valid"
        assert self.valid_angles() == True, "Angles are not valid"
        assert self.valid_torsional() == True, "Torsional are not valid"
        self.create_graph()
    
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
    def bond_ff_type(self, bond_id, bond_type):
        '''
         Given when a bond is defined, it will have a type and length associated with it.
        '''
        for bond in self.bonds:
            if bond[0] == bond_id:
                bond[1] = bond_type
                return True
        return False
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
    def save_tojson(self, filename):
        '''
        Saves the molecule to a json file.
        '''
        with open(filename, 'w') as f:
            json.dump(self.__dict__, f, indent=4)
            
            
    #--------------------------------
#===============================================
class Molecule_Type:
    def __init__(self, topo_file):
        self.topo = Molecule_Topo(topo_file)
        self.nbonds = len(self.topo.bonds)
        self.bondtype = np.zeros(self.nbonds)
        self.nangles = len(self.topo.angles)
        self.angletype = np.zeros(self.nangles)
        self.ntorsions = len(self.topo.torsion)
        self.torsiontype = np.zeros(self.ntorsions)
        self.n_atoms = len(self.topo.atoms)
        self.n_bonds = len(self.topo.bonds)
        self.n_angles = len(self.topo.angles)
        self.n_torsions = len(self.topo.torsion)
        
        
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
