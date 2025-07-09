import numpy as np

#=================================================================
class Displacement():
    #-------------------------------------------------------------
    def __init__(self, atoms: list, 
                 newpositions: np.ndarray,
                 singlemolecule: bool = True):
        """
        Initialize the Displacement object with atoms and their displacements.
        :param atoms: List of atom indices.
        :param disp: List of displacements for each atom.
        """
        self.atoms = atoms
        self.newpositions = newpositions
        self.singlemolecule = singlemolecule
    #-------------------------------------------------------------
    def update(self, boxatoms):
        boxatoms[self.atoms] = self.newpositions
        return boxatoms
    #-------------------------------------------------------------
#=================================================================