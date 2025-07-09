"""
AVBMC Monte Carlo Move
Corresponds to Move_MC_AVBMC.f90

The AVBMC move is the intra-box variant of the move created by Bin Chen et al. 
This is a biased swap move that increases the chances of creating and breaking 
bonded configurations by using a small volume around the molecules present in the 
system as a "target" location for a newly inserted molecule.
"""

import numpy as np
import sys
import math
from random import random
from Template_MCMove import MCMove
from CoordinateTypes import Addition, Deletion
from VarPrecision import dp


class AVBMC(MCMove):
    """
    Aggregation Volume Bias Monte Carlo move.
    Corresponds to the Fortran AVBMC type.
    """
    
    def __init__(self):
        super().__init__()
        
        # AVBMC specific parameters
        self.ebias = False
        self.neilistindx = 1
        self.inatmps = 1e-30
        self.inaccpt = 0.0
        self.outatmps = 1e-30
        self.outaccpt = 0.0
        self.avbmcRad = 4.0
        self.avbmcRadSq = 16.0
        self.avbmcVol = 0.0
        
        # Displacement arrays
        self.newPart = []  # Will be allocated in prologue
        self.oldPart = [Deletion()]
        
        # Insertion point arrays
        self.insPoint = None
        self.insProb = None
    
    def constructor(self):
        """
        Corresponds to AVBMC_Constructor
        Initialize the move
        """
        # Initialize box probabilities if not set
        if len(self.boxProb) == 0:
            # Default to single box
            self.boxProb = np.array([1.0], dtype=dp)
    
    def full_move(self, trial_box):
        """
        Corresponds to AVBMC_FullMove
        Perform an AVBMC move (either swap in or swap out)
        
        Args:
            trial_box: SimBox instance
            
        Returns:
            bool: Whether the move was accepted
        """
        if random() > 0.5:
            return self.swap_in(trial_box)
        else:
            return self.swap_out(trial_box)
    
    def swap_in(self, trial_box, moveid=None, targid=None, ProbIn=None):
        """
        Corresponds to AVBMC_SwapIn
        Perform a swap-in move (add molecule)
        
        Args:
            trial_box: SimBox instance
            moveid: Optional output for move ID
            targid: Optional output for target ID
            ProbIn: Optional output for probability
            
        Returns:
            bool: Whether the move was accepted
        """
        self.atmps += 1.0
        self.inatmps += 1.0
        accept = True
        
        self.load_box_info(trial_box, self.newPart)
        
        # Select molecule type to add
        from Common_MolInfo import nMolTypes
        nType = int(nMolTypes * random() + 1)
        
        if trial_box.NMol[nType - 1] + 1 > trial_box.NMolMax[nType - 1]:
            accept = False
            if ProbIn is not None:
                ProbIn = 0.0
            return accept
        
        # Get molecule data
        nMove = trial_box.NMol[nType - 1] + 1
        nMove = trial_box.MolGlobalIndx[nType - 1, nMove - 1]
        mol_data = trial_box.get_mol_data(nMove)
        molType = mol_data['molType']
        molStart = mol_data['molStart']
        molEnd = mol_data['molEnd']
        
        # Select target molecule
        nTarget = self.uniform_molecule_select(trial_box)
        targ_data = trial_box.get_mol_data(nTarget)
        targStart = targ_data['molStart']
        targEnd = targ_data['molEnd']
        
        # Get target atom coordinates
        target_atoms = trial_box.get_coordinates(slice_range=(targStart, targEnd))
        
        # Generate insertion points
        from Common_MolInfo import MolData
        nInsPoints = MolData[molType]['molConstruct'].GetNInsertPoints()
        self.allocate_prob(nInsPoints)
        
        for iIns in range(nInsPoints):
            self.create_forward(target_atoms, self.insPoint[:, iIns])
            self.insProb[iIns] = 1.0
        
        # Set up new molecule atoms
        nAtoms = MolData[molType]['nAtoms']
        for iAtom in range(nAtoms):
            atomIndx = molStart + iAtom - 1
            self.newPart[iAtom].molType = nType - 1
            self.newPart[iAtom].molIndx = nMove
            self.newPart[iAtom].atmIndx = atomIndx
        
        # Generate molecule configuration
        ProbSub, accept = MolData[molType]['molConstruct'].GenerateConfig(
            trial_box, self.newPart[:nAtoms], self.insPoint, self.insProb
        )
        if not accept:
            if ProbIn is not None:
                ProbIn = 0.0
            return accept
        
        GenProb = ProbSub
        
        # Set up neighbor lists
        for iAtom in range(nAtoms):
            # This would need proper neighbor list implementation
            self.newPart[iAtom].listIndex = iAtom
        
        # Check constraints
        if not trial_box.check_constraint(self.newPart[:nAtoms]):
            if ProbIn is not None:
                ProbIn = 0.0
            return False
        
        # Calculate energy
        e_inter, e_intra, e_diff, accept = trial_box.compute_energy_delta(
            self.newPart[:nAtoms], self.tempList, self.tempNnei, computeintra=True
        )
        if not accept:
            if ProbIn is not None:
                ProbIn = 0.0
            return accept
        
        # Check post-energy constraints
        if not trial_box.check_post_energy(self.newPart[:nAtoms], e_diff):
            if ProbIn is not None:
                ProbIn = 0.0
            return False
        
        # Calculate reverse probability
        ProbSel = 0.0
        self.select_neigh_reverse(trial_box, nTarget, ProbSel)
        
        # Calculate generation probability
        GasProb = MolData[molType]['molConstruct'].GasConfig()
        
        # Calculate acceptance probability
        Prob = float(trial_box.nMolTotal) * self.avbmcVol
        Prob = Prob * ProbSel / float(trial_box.nMolTotal + 1)
        Prob = GasProb * Prob / GenProb
        
        # Get extra terms
        from CommonSampling import sampling
        if sampling is not None:
            extra_terms = sampling.get_extra_terms(self.newPart[:nAtoms], trial_box)
        else:
            extra_terms = 0.0
        
        # Accept/reject
        if sampling is not None:
            accept = sampling.make_decision(trial_box, e_diff, self.newPart[:nAtoms], in_prob=Prob, extra_in=extra_terms)
        else:
            accept = True
        
        if moveid is not None:
            accept = True
        
        if accept:
            self.accpt += 1.0
            self.inaccpt += 1.0
            if moveid is not None:
                moveid = self.newPart[0].molIndx
            if targid is not None:
                targid = nTarget
            if ProbIn is not None:
                ProbIn = Prob
            
            trial_box.update_energy(e_diff, e_inter, e_intra)
            trial_box.update_position(self.newPart[:nAtoms], self.tempList, self.tempNnei)
        
        return accept
    
    def swap_out(self, trial_box, forcetargid=None, forceid=None, ProbOut=None):
        """
        Corresponds to AVBMC_SwapOut
        Perform a swap-out move (remove molecule)
        
        Args:
            trial_box: SimBox instance
            forcetargid: Optional forced target ID
            forceid: Optional forced molecule ID
            ProbOut: Optional output for probability
            
        Returns:
            bool: Whether the move was accepted
        """
        self.atmps += 1.0
        self.outatmps += 1.0
        accept = True
        
        self.load_box_info(trial_box, self.oldPart)
        
        # Select target molecule
        if forcetargid is not None:
            nTarget = forcetargid
        else:
            nTarget = self.uniform_molecule_select(trial_box)
        
        # Select molecule to remove
        if forceid is not None:
            nMove = forceid
            ProbSel = 1.0
        else:
            nMove, ProbSel = self.select_neigh(trial_box, nTarget)
        
        if nMove < 1:
            accept = False
            if ProbOut is not None:
                ProbOut = 0.0
            return accept
        
        mol_data = trial_box.get_mol_data(nMove)
        molType = mol_data['molType']
        
        if trial_box.NMol[molType] - 1 < trial_box.NMolMin[molType]:
            accept = False
            if ProbOut is not None:
                ProbOut = 0.0
            return accept
        
        self.oldPart[0].molType = molType
        self.oldPart[0].molIndx = nMove
        
        # Check constraints
        if not trial_box.check_constraint(self.oldPart):
            if ProbOut is not None:
                ProbOut = 0.0
            return False
        
        # Calculate energy
        e_inter, e_intra, e_diff, accept = trial_box.compute_energy_delta(
            self.oldPart, self.tempList, self.tempNnei, computeintra=True
        )
        if not accept:
            if ProbOut is not None:
                ProbOut = 0.0
            return accept
        
        # Check post-energy constraints
        if not trial_box.check_post_energy(self.oldPart, e_diff):
            if ProbOut is not None:
                ProbOut = 0.0
            return False
        
        # Get target coordinates for reverse move
        targ_data = trial_box.get_mol_data(nTarget)
        targStart = targ_data['molStart']
        targEnd = targ_data['molEnd']
        target_atoms = trial_box.get_coordinates(slice_range=(targStart, targEnd))
        
        # Generate insertion points for reverse move
        from Common_MolInfo import MolData
        nInsPoints = MolData[molType]['molConstruct'].GetNInsertPoints()
        self.allocate_prob(nInsPoints)
        
        for iIns in range(nInsPoints):
            self.create_forward(target_atoms, self.insPoint[:, iIns])
            self.insProb[iIns] = 1.0
        
        # Calculate reverse generation probability
        ProbSub, accept = MolData[molType]['molConstruct'].ReverseConfig(
            self.oldPart, trial_box, self.insPoint, self.insProb
        )
        GenProb = ProbSub
        
        # Calculate gas probability
        GasProb = MolData[molType]['molConstruct'].GasConfig()
        
        # Calculate acceptance probability
        Prob = float(trial_box.nMolTotal)
        Prob = Prob / (float(trial_box.nMolTotal - 1) * self.avbmcVol * ProbSel)
        Prob = Prob * GenProb / GasProb
        
        # Get extra terms
        from CommonSampling import sampling
        if sampling is not None:
            extra_terms = sampling.get_extra_terms(self.oldPart, trial_box)
        else:
            extra_terms = 0.0
        
        # Accept/reject
        if sampling is not None:
            accept = sampling.make_decision(trial_box, e_diff, self.oldPart, in_prob=Prob, extra_in=extra_terms)
        else:
            accept = True
        
        if forceid is not None:
            accept = True
        
        if accept:
            self.accpt += 1.0
            self.outaccpt += 1.0
            if ProbOut is not None:
                ProbOut = Prob
            
            trial_box.update_energy(e_diff, e_inter, e_intra)
            trial_box.delete_mol(self.oldPart[0].molIndx)
        
        return accept
    
    def allocate_prob(self, nInsPoints):
        """
        Corresponds to AVBMC_AllocateProb
        Allocate insertion point arrays
        """
        if self.insPoint is None or self.insPoint.shape[1] < nInsPoints:
            self.insPoint = np.zeros((3, nInsPoints), dtype=dp)
            self.insProb = np.zeros(nInsPoints, dtype=dp)
    
    def create_forward(self, target_atoms, inspoint):
        """
        Corresponds to AVBMC_CreateForward
        Create forward insertion point
        """
        # Generate random direction
        from RandomGen import Generate_UnitSphere
        dx, dy, dz = Generate_UnitSphere()
        
        # Generate random radius
        radius = self.avbmcRad * random() ** (1.0 / 3.0)
        
        # Calculate insertion point
        inspoint[0] = target_atoms[0, 0] + radius * dx
        inspoint[1] = target_atoms[1, 0] + radius * dy
        inspoint[2] = target_atoms[2, 0] + radius * dz
    
    def select_neigh(self, trial_box, targetIndx, removeMolIndx=None, ProbOut=None, forceid=None):
        """
        Corresponds to AVBMC_SelectNeigh
        Select neighbor molecule for removal
        """
        # This is a simplified implementation
        # In practice, this would need proper neighbor list handling
        
        # For now, just return a random molecule
        import random
        nMove = random.randint(1, trial_box.nMolTotal)
        ProbSel = 1.0 / trial_box.nMolTotal
        
        if removeMolIndx is not None:
            removeMolIndx = nMove
        if ProbOut is not None:
            ProbOut = ProbSel
        
        return nMove, ProbSel
    
    def select_neigh_reverse(self, trial_box, targetIndx, ProbOut):
        """
        Corresponds to AVBMC_SelectNeighReverse
        Select neighbor for reverse move
        """
        # Simplified implementation
        ProbOut = 1.0 / (trial_box.nMolTotal + 1)
    
    def prologue(self):
        """
        Corresponds to AVBMC_Prologue
        Initialize the move before simulation
        """
        # Calculate AVBMC volume
        self.avbmcVol = (4.0 / 3.0) * math.pi * self.avbmcRad ** 3
        self.avbmcRadSq = self.avbmcRad * self.avbmcRad
        
        print(f"AVBMC Radius: {self.avbmcRad}")
        print(f"AVBMC Volume: {self.avbmcVol}")
        
        # Allocate newPart array
        from Common_MolInfo import MolData, mostAtoms
        maxAtoms = 0
        for mol in MolData:
            if mol['nAtoms'] > maxAtoms:
                maxAtoms = mol['nAtoms']
        
        self.newPart = [Addition() for _ in range(maxAtoms)]
        self.create_temp_array(mostAtoms)
    
    def epilogue(self):
        """
        Corresponds to AVBMC_Epilogue
        Finalize the move after simulation
        """
        print(f"AVBMC Moves Accepted: {int(self.accpt):15d}")
        print(f"AVBMC Moves Attempted: {int(self.atmps):15d}")
        
        accpt_rate = self.get_accept_rate()
        print(f"AVBMC Acceptance Rate: {accpt_rate:15.8f}")
        
        print(f"AVBMC Out Moves Accepted: {int(self.outaccpt):15d}")
        print(f"AVBMC Out Moves Attempted: {int(self.outatmps):15d}")
        if self.outatmps > 0:
            out_rate = 100.0 * self.outaccpt / self.outatmps
            print(f"AVBMC Out Acceptance Rate: {out_rate:15.8f}")
        
        print(f"AVBMC In Moves Accepted: {int(self.inaccpt):15d}")
        print(f"AVBMC In Moves Attempted: {int(self.inatmps):15d}")
        if self.inatmps > 0:
            in_rate = 100.0 * self.inaccpt / self.inatmps
            print(f"AVBMC In Acceptance Rate: {in_rate:15.8f}")
    
    def process_io(self, line: str) -> int:
        """
        Corresponds to AVBMC_ProcessIO
        Process input commands
        """
        parts = line.strip().split()
        if len(parts) < 5:
            return -1
        
        command = parts[3].lower()
        value_str = parts[4]
        
        try:
            if command == "energybias":
                self.ebias = value_str.lower() in ['true', '.true.', 't', '1']
            elif command == "neighlist":
                self.neilistindx = int(value_str)
            elif command == "radius":
                self.avbmcRad = float(value_str)
                self.avbmcRadSq = self.avbmcRad * self.avbmcRad
                self.avbmcVol = (4.0 / 3.0) * math.pi * self.avbmcRad ** 3
            else:
                return -1
        except (ValueError, IndexError):
            return -1
        
        return 0 