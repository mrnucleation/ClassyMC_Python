# SimMonteCarlo.py
import time

# Note: The following imports assume corresponding Python modules exist.


from src_python

from random import choice
from typing import List, Optional



# =============================================================================
class SimMonteCarlo:
    def __init__(self, nCycles, nMoves, screenfreq, configfreq, energyCheck):
        self.nCycles = nCycles
        self.nMoves = nMoves
        self.screenfreq = screenfreq
        self.configfreq = configfreq
        self.energyCheck = energyCheck
        
        # Initialize required attributes that were missing
        self.BoxList = []           # List of simulation boxes
        self.Moves = []             # List of Monte Carlo moves
        self.Sampling = None        # Sampling rule object
        self.AnalysisArray = None   # Analysis functions
        self.TrajArray = None       # Trajectory output functions
        self.EnergyCalculator = []  # Energy calculation functions
        self.MolData = []           # Molecular data
        self.TimeStart = time.time()  # Simulation start time
    
    # ---------------------------------------------------------------------------
    def run_monte_carlo(self):
        iCycle = 0
        iMove = 0
        print("Starting Pre-Simulation Checks....")
        self.prologue()
        self.safety_check()
        nBoxes = len(self.BoxList)
        boxProb = [0] * nBoxes
        boxNum = 1
        accept = True

        # Initial analysis calls
        self.analyze(iCycle, iMove, accept, True)
        self.analyze(iCycle, iMove, accept, False)

        if self.nCycles < 1:
            print("============================================")
            print("Number of Cycles is less than 1!")
            print("Run Command was issued, but nothing is done")
            print("============================================")
            return

        print("============================================")
        print("       Simulation Start!")
        print("============================================")

        time_start = time.time()
        self.screen_out(iCycle, iMove)
        screenfreq = 100
        configfreq = 1000
        energyCheck = 1000

        for iCycle in range(self.nCycles):
            for iMove in range(self.nMoves):
                accept = True
                if len(self.Moves) > 0:
                    moveNum = choice(range(len(self.Moves)))
                    curmove = self.Moves[moveNum]

                    # Perform selected move
                    self.Moves[moveNum].GetBoxProb(boxProb)
                    box = choice(self.BoxList, p=boxProb)
                    accept = curmove.FullMove(box, self.Sampling, accept)

                    if accept:
                        self.update(accept)

                    # Per-move analysis and neighbor-list checks
                    self.analyze(iCycle, iMove, accept, True)
                    if accept:
                        for iBox in range(nBoxes):
                            if boxNum < 0 or boxNum == iBox + 1:
                                self.BoxList[iBox].CheckLists()

            # Periodic outputs and checks
            if iCycle % screenfreq == 0:
                self.screen_out(iCycle, iMove)
            if energyCheck > 0 and iCycle % energyCheck == 0:
                for boxIdx in range(nBoxes):
                    self.BoxList[boxIdx].EnergySafetyCheck()

            # Per-cycle analysis and maintenance
            self.analyze(iCycle, iMove, accept, False)
            self.maintenance(iCycle, iMove)
            self.trajectory(iCycle, iMove)

        # End of main loop
        self.screen_out(iCycle, iMove)
        print("=======================================")
        print("     Simulation End")
        print("=======================================")
        print("     Beginning Epilogue.....")

        self.epilogue()

        # Finalize analysis
        if self.AnalysisArray is not None:
            for analysis in self.AnalysisArray:
                analysis.func.Finalize()
    # ---------------------------------------------------------------------------
    
    # ---------------------------------------------------------------------------
    def analyze(self, iCycle, iMove, accept, moveloop):
        if self.AnalysisArray is not None:
            for analysis in self.AnalysisArray:
                if analysis.func.perMove == moveloop:
                    if (iCycle % analysis.func.UpdateFreq == 0) or analysis.func.perMove:
                        analysis.func.Compute(accept)
    # ---------------------------------------------------------------------------
    
    # ---------------------------------------------------------------------------
    def screen_out(self, iCycle, iMove):
        print(f"    ----------------- Cycle Number: {iCycle} ----------------------")
        current_time = time.time()
        print(f"   Simulation Time: {current_time - self.TimeStart:.3f} sec")
        if iCycle > 0:
            time_per_cycle = (current_time - self.TimeStart) / iCycle
            if time_per_cycle > 1e-3:
                print(f"   Time per Cycle: {time_per_cycle:.3f} sec")
            else:
                print(f"   Time per Cycle: {time_per_cycle:.3e} sec")

        if self.Sampling is not None:
            self.Sampling.ScreenOut()
        if self.AnalysisArray is not None:
            for analysis in self.AnalysisArray:
                analysis.func.ScreenOut()
        if self.TrajArray is not None:
            for traj in self.TrajArray:
                traj.ScreenOut()

        for ecalc in self.EnergyCalculator:
            ecalc.ScreenOut()
        for box in self.BoxList:
            box.ScreenOut()
        for move in self.Moves:
            move.ScreenOut()
        print()
    # ---------------------------------------------------------------------------
    
    # ---------------------------------------------------------------------------
    def maintenance(self, iCycle, iMove):
        if self.Sampling is not None and hasattr(self.Sampling, 'maintFreq'):
            if iCycle % self.Sampling.maintFreq == 0:
                self.Sampling.Maintenance()

        if self.AnalysisArray is not None:
            for analysis in self.AnalysisArray:
                if iCycle % analysis.func.maintFreq == 0:
                    analysis.func.Maintenance()

        if self.TrajArray is not None:
            for traj in self.TrajArray:
                if iCycle % traj.traj.maintFreq == 0:
                    traj.traj.Maintenance()

        for box in self.BoxList:
            box.box.Maintenance()
        for move in self.Moves:
            if iCycle % move.move.maintFreq == 0:
                move.move.Maintenance()
    # ---------------------------------------------------------------------------
    
    # ---------------------------------------------------------------------------
    def trajectory(self, iCycle, iMove):
        if self.TrajArray is not None:
            for traj in self.TrajArray:
                if iCycle % traj.traj.outfreq == 0:
                    traj.WriteFrame(iCycle)
    # ---------------------------------------------------------------------------
    
    # ---------------------------------------------------------------------------
    def update(self, accept):
        assert isinstance(accept, bool), "Accept must be a boolean value"
        if not accept:
            return
        if self.Sampling is not None:
            self.Sampling.Update()
        if self.AnalysisArray is not None:
            for analysis in self.AnalysisArray:
                analysis.Update()
        if self.TrajArray is not None:
            for traj in self.TrajArray:
                traj.Update()
        for ecalc in self.EnergyCalculator:
            ecalc.Update()
        for box in self.BoxList:
            box.Update()
        for move in self.Moves:
            move.Update()
    # ---------------------------------------------------------------------------
    
    # ---------------------------------------------------------------------------
    def safety_check(self):
        if not self.MolData:
            raise RuntimeError("CRITICAL ERROR! Molecular Topology Information has not been defined!")
        for mol in self.MolData:
            if not getattr(mol, 'molConstruct', None):
                print("WARNING! Molecule reconstructor is not defined in the forcefield file!")
            else:
                mol.molConstruct.SafetyCheck()

        if not self.Sampling:
            raise RuntimeError("CRITICAL ERROR! Sampling Rule has not been defined!")
        self.Sampling.SafetyCheck()

        if self.AnalysisArray is not None:
            for analysis in self.AnalysisArray:
                analysis.func.SafetyCheck()
        if self.TrajArray is not None:
            for traj in self.TrajArray:
                traj.traj.SafetyCheck()

        if not self.BoxList:
            raise RuntimeError("CRITICAL ERROR! No Simulation Boxes have been defined!")
        for box in self.BoxList:
            box.box.SafetyCheck()

        if self.Moves:
            for move in self.Moves:
                move.move.SafetyCheck()
        else:
            print("WARNING! No Monte Carlo moves have been defined!")
    # ---------------------------------------------------------------------------
    
    # ---------------------------------------------------------------------------
    def prologue(self):
        for mol in self.MolData:
            if getattr(mol, 'molConstruct', None):
                mol.molConstruct.Prologue()
            if hasattr(mol, 'miscdata'):
                for misc in mol.MiscData:
                    misc.miscFF.Prologue()

        if self.AnalysisArray is not None:
            for analysis in self.AnalysisArray:
                analysis.func.Prologue()
        if self.Sampling is not None:
            self.Sampling.Prologue()

        if self.TrajArray is not None:
            for traj in self.TrajArray:
                traj.traj.Prologue()

        for ecalc in self.EnergyCalculator:
            ecalc.method.Prologue()
        for box in self.BoxList:
            box.box.Prologue()
        for move in self.Moves:
            move.move.Prologue()
    # ---------------------------------------------------------------------------
    
    # ---------------------------------------------------------------------------
    def epilogue(self):
        for box in self.BoxList:
            box.Epilogue()
        print("-----------------------")

        if self.Sampling is not None:
            self.Sampling.Epilogue()
        for mol in self.MolData:
            if getattr(mol, 'molConstruct', None):
                mol.Epilogue()
            if hasattr(mol, 'miscdata'):
                for misc in mol.MiscData:
                    misc.Epilogue()

        if self.AnalysisArray is not None:
            for analysis in self.AnalysisArray:
                analysis.Epilogue()
        if self.TrajArray is not None:
            for traj in self.TrajArray:
                traj.Epilogue()

        for move in self.Moves:
            move.Epilogue()
    # ---------------------------------------------------------------------------
# =============================================================================
