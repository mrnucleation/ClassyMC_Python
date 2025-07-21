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
        #self.name = json_data['name']
        self.atoms = json_data['atoms']
        
        for atom in self.atoms:
            assert len(atom) == 2, "Each atom must have an element ID, type"
        
        
        # Initialize bonds as empty list if not present in json_data
        self.bonds = json_data.get('bonds', [])
        if len(self.atoms) > 1 and self.bonds:
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
        
        # If there is only one atom, no need to check the graph
        if len(self.atoms) == 1:
            return 
        
        if len(self.graph.nodes()) != len(self.atoms):
            raise ValueError("Number of atoms in the graph does not match the number of atoms in the molecule definition.")
        
        if len(self.graph.edges()) != len(self.bonds):
            raise ValueError("Number of bonds in the graph does not match the number of bonds in the molecule definition.")
        
        if not nx.is_connected(G):
            raise ValueError("Molecule graph is not connected. All atoms must be part of a bond.")
        
 
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
    def __init__(self, topo_info, atomtypes, bondtype=None, angletype=None, torsiontype=None):
        if isinstance(topo_info, str):
            if not topo_info.endswith('.json'):
                raise ValueError("Topology file must be a JSON file")
            # Open the topology file and read the JSON data
            with open(topo_info, 'r') as f:
                topo_data = json.load(f)
            self.topo = Molecule_Topo(topo_data)
        elif isinstance(topo_info, dict):
            self.topo = Molecule_Topo(topo_info)
        else:
            raise TypeError("Invalid topo_info type")

        self.n_atoms = len(self.topo.atoms)
        self.n_bonds = len(self.topo.bonds)
        self.n_angles = len(self.topo.angles)
        self.n_torsions = len(self.topo.torsion)
        
        # Set default 
        self.atomtypes = atomtypes
        assert isinstance(self.atomtypes, list), "atomtypes must be a list"

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
#=================================================
