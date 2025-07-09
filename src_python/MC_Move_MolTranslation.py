from .Template_MCMove import MCMove
import numpy as np
from random import random
from .VarPrecision import dp
from .Box_SimpleBox import SimpleBox
import math

class MolTranslate(MCMove):
    """
    Performs a Monte Carlo translation move on a randomly selected molecule.

    This class is a Python translation of the Fortran `MCMove_MolTranslation`
    module.
    """
    def __init__(self, BoxArray):
        super().__init__()
        # --- Default attributes from Fortran type definition ---
        self.verbose = True
        self.proportional = True
        self.tuneMax = True
        self.limit = 5.0
        self.targAccpt = 50.0
        self.max_dist = 0.05

        # --- Rejection Counters ---
        self.ovlaprej = 0
        self.constrainrej = 0
        self.detailedrej = 0

        # --- "Allocatable" arrays, initialized below ---
        self.boxatmps = np.array([], dtype=dp)
        self.boxaccpt = np.array([], dtype=dp)
        self.boxlimit = np.array([], dtype=dp)
        self.boxtargAccpt = np.array([], dtype=dp)
        self.boxmax_dist = np.array([], dtype=dp)
        self.disp = []  # Will be a list of Displacement objects

        # --- Constructor logic from MolTrans_Constructor ---
        self._initialize_arrays()

    def _initialize_arrays(self):
        """
        Initializes arrays and other variables, equivalent to the Fortran
        `MolTrans_Constructor` subroutine.
        """
        n_boxes = len(BoxArray)
        if self.boxProb.size == 0:
            self.boxProb = np.full(n_boxes, 1.0 / n_boxes, dtype=dp)

        self.boxatmps = np.full(n_boxes, 1e-50, dtype=dp)
        self.boxaccpt = np.zeros(n_boxes, dtype=dp)

        self.boxlimit = np.full(n_boxes, self.limit, dtype=dp)
        self.boxmax_dist = np.full(n_boxes, self.max_dist, dtype=dp)
        self.boxtargAccpt = np.full(n_boxes, self.targAccpt, dtype=dp)

        max_atoms = 0
        if nMolTypes > 0:
            max_atoms = max(info.nAtoms for info in MolData.values())

        self.disp = [Displacement() for _ in range(max_atoms)]
        self.CreateTempArray(max_atoms)

    def full_move(self, trial_box: SimpleBox):
        """
        Performs one molecule translation MC move.
        Returns a boolean indicating if the move was accepted.
        """
        box_idx = trial_box.boxID - 1  # Convert 1-based ID to 0-based index

        self.atmps += 1.0
        self.boxatmps[box_idx] += 1.0

        # --- Propose move ---
        raw_indx = math.floor(trial_box.nMolTotal * random() + 1.0)
        n_move = FindMolecule(trial_box, raw_indx)
        mol_start, mol_end, mol_type = trial_box.GetMolData(n_move)
        mol_start_idx = mol_start - 1

        dx = self.boxmax_dist[box_idx] * (2.0 * random() - 1.0)
        dy = self.boxmax_dist[box_idx] * (2.0 * random() - 1.0)
        dz = self.boxmax_dist[box_idx] * (2.0 * random() - 1.0)

        n_atoms = MolData[mol_type].nAtoms
        for i_atom in range(n_atoms):
            atom_idx = mol_start_idx + i_atom
            d = self.disp[i_atom]
            d.molType, d.molIndx, d.atmIndx = mol_type, n_move, atom_idx + 1
            d.x_new = trial_box.atoms[0, atom_idx] + dx
            d.y_new = trial_box.atoms[1, atom_idx] + dy
            d.z_new = trial_box.atoms[2, atom_idx] + dz
            d.newlist, d.listIndex = False, i_atom + 1

        # --- Check constraints and calculate energy ---
        active_disp = self.disp[:n_atoms]
        if not trial_box.CheckConstraint(active_disp):
            self.constrainrej += 1
            return False

        e_inter, e_intra, e_diff, accept_energy = trial_box.ComputeEnergyDelta(
            active_disp, self.tempList, self.tempNnei, computeintra=False
        )
        if not accept_energy:
            self.ovlaprej += 1
            return False

        if not trial_box.CheckPostEnergy(active_disp, e_diff):
            self.constrainrej += 1
            return False

        # --- Accept/Reject ---
        accept = sampling.MakeDecision(trial_box, e_diff, active_disp, inProb=1.0)

        if accept:
            self.accpt += 1.0
            self.boxaccpt[box_idx] += 1.0
            trial_box.UpdateEnergy(e_diff, e_inter)
            trial_box.UpdatePosition(active_disp, self.tempList, self.tempNnei)
        else:
            self.detailedrej += 1
        
        return accept

    def maintenance(self):
        """
        Tunes the maximum displacement to achieve a target acceptance rate.
        """
        if not self.tuneMax:
            return

        for i_box in range(len(self.boxatmps)):
            if self.boxatmps[i_box] < 0.5:
                continue
            
            acc_rate = 100.0 * self.boxaccpt[i_box] / self.boxatmps[i_box]

            if acc_rate > self.boxtargAccpt[i_box]:
                new_dist = self.boxmax_dist[i_box] * 1.01
                self.boxmax_dist[i_box] = min(new_dist, self.boxlimit[i_box])
            else:
                self.boxmax_dist[i_box] *= 0.99

    def prologue(self):
        """Prints information at the start of a simulation block."""
        print(f"(Molecule Translate) Maximum Displacement: {self.max_dist:15.8f}", file=nout)

    def epilogue(self):
        """Prints summary statistics at the end of a simulation block."""
        print(file=nout)
        print(f"Molecule Translation Moves Accepted:  {round(self.accpt):15d}", file=nout)
        print(f"Molecule Translation Moves Attempted: {round(self.atmps):15d}", file=nout)
        
        accpt_rate = self.GetAcceptRate()
        print(f"Molecule Translation Acceptance Rate: {accpt_rate:15.8f}", file=nout)
        
        if self.tuneMax:
            dist_str = " ".join([f"{d:15.8f}" for d in self.boxmax_dist])
            print(f"Final Maximum Displacement: {dist_str}", file=nout)

        if self.verbose:
            print(f"Molecule Translation, Rejections due to overlap:         {self.ovlaprej:15d}", file=nout)
            print(f"Molecule Translation, Rejections due to constraint:      {self.constrainrej:15d}", file=nout)
            print(f"Molecule Translation, Rejections due to detailed balance:{self.detailedrej:15d}", file=nout)

    def update(self):
        """Updates box probabilities based on the number of molecules."""
        if self.proportional:
            n_mols = np.array([bw.box.nMolTotal for bw in BoxArray], dtype=dp)
            norm = np.sum(n_mols)
            if norm > 0:
                self.boxProb = n_mols / norm

    def process_io(self, line: str):
        """
        Parses a line from an input file to set parameters.
        Assumes format like: "... ... ... command value"
        Returns 0 for success, -1 for failure.
        """
        parts = line.strip().lower().split()
        if len(parts) < 5: return -1

        command, value_str = parts[3], parts[4]
        try:
            if command == "tunemax":
                self.tuneMax = value_str.lower() in ['true', '.true.', 't', '1']
            elif command == "dynamiclimit":
                self.limit = float(value_str)
            elif command == "dynamictarget":
                self.targAccpt = float(value_str)
            elif command == "maxdisplace":
                self.max_dist = float(value_str)
            elif command == "proportional":
                self.proportional = value_str.lower() in ['true', '.true.', 't', '1']
            elif command == "updatefreq":
                self.maintFreq = int(value_str)
            else:
                return -1  # Command not recognized
        except (ValueError, IndexError):
            return -1  # Error parsing value
        return 0

# =============================================================================
# Example Usage
# =============================================================================
if __name__ == '__main__':
    # 1. Setup the mock simulation environment
    box1_atoms = np.random.rand(3, 200) * 10
    box2_atoms = np.random.rand(3, 300) * 10
    box1 = SimpleBox(boxID=1, nMolTotal=100, atoms=box1_atoms)
    box2 = SimpleBox(boxID=2, nMolTotal=150, atoms=box2_atoms)
    BoxArray = [BoxWrapper(box1), BoxWrapper(box2)]
    MolData = {1: MolInfo(nAtoms=2)}
    nMolTypes = 1

    # 2. Create an instance of the move
    mol_translate_move = MolTranslate()

    # 3. Process some IO commands (mimicking input file parsing)
    mol_translate_move.process_io("MC_MOVE MOLTRANS TRANSLATE tunemax true")
    mol_translate_move.process_io("MC_MOVE MOLTRANS TRANSLATE maxdisplace 0.1")
    mol_translate_move.process_io("MC_MOVE MOLTRANS TRANSLATE updatefreq 50")
    print(f"Initial max_dist set by IO: {mol_translate_move.max_dist}")
    print(f"Initial boxmax_dist: {mol_translate_move.boxmax_dist}")
    print("-" * 20)

    # 4. Run a mock simulation loop
    mol_translate_move.prologue()
    for i in range(1, 1001):
        # In a real simulation, a box would be chosen based on boxProb
        chosen_box = np.random.choice(BoxArray, p=mol_translate_move.boxProb).box
        mol_translate_move.full_move(chosen_box)
        
        if i % mol_translate_move.maintFreq == 0:
            mol_translate_move.maintenance()
            mol_translate_move.update()

    print("-" * 20)
    mol_translate_move.epilogue()
