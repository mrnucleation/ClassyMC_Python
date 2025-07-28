# SimMonteCarlo.py
import time

# Note: The following imports assume corresponding Python modules exist.

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
            #if energyCheck > 0 and iCycle % energyCheck == 0:
            #    for boxIdx in range(nBoxes):
            #        self.BoxList[boxIdx].Energysafety_check()

            # Per-cycle analysis and maintenance
            self.analyze(iCycle, iMove, accept, False)
            self.maintenance(iCycle, iMove)
            self.trajectory(iCycle, iMove)

        # End of main loop
        self.screen_out(iCycle, iMove)
        print("=======================================")
        print("     Simulation End")
        print("=======================================")
        print("     Beginning epilogue.....")

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
                    if (iCycle % analysis.func.updateFreq == 0) or analysis.func.perMove:
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
            self.Sampling.screenout()
        if self.AnalysisArray is not None:
            for analysis in self.AnalysisArray:
                analysis.func.screenout()
        if self.TrajArray is not None:
            for traj in self.TrajArray:
                traj.screenout()

        for ecalc in self.EnergyCalculator:
            ecalc.screenout()
        for box in self.BoxList:
            box.screenout()
        for move in self.Moves:
            move.screenout()
        print()
    # ---------------------------------------------------------------------------
    
    # ---------------------------------------------------------------------------
    def maintenance(self, iCycle, iMove):
        if self.Sampling is not None and hasattr(self.Sampling, 'maintFreq'):
            if iCycle % self.Sampling.maintFreq == 0:
                self.Sampling.maintenance()

        if self.AnalysisArray is not None:
            for analysis in self.AnalysisArray:
                if iCycle % analysis.func.maintFreq == 0:
                    analysis.func.maintenance()

        if self.TrajArray is not None:
            for traj in self.TrajArray:
                if iCycle % traj.traj.maintFreq == 0:
                    traj.traj.maintenance()

        for box in self.BoxList:
            box.maintenance()
        for move in self.Moves:
            if iCycle % move.move.maintFreq == 0:
                move.move.maintenance()
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
            self.Sampling.update()
        if self.AnalysisArray is not None:
            for analysis in self.AnalysisArray:
                analysis.update()
        if self.TrajArray is not None:
            for traj in self.TrajArray:
                traj.update()
        for ecalc in self.EnergyCalculator:
            ecalc.update()
        for box in self.BoxList:
            box.update()
        for move in self.Moves:
            move.update()
    # ---------------------------------------------------------------------------
    
    # ---------------------------------------------------------------------------
    def safety_check(self):
        #if not self.MolData:
        #    raise RuntimeError("CRITICAL ERROR! Molecular Topology Information has not been defined!")
        #for mol in self.MolData:
        #    if not getattr(mol, 'molConstruct', None):
        #        print("WARNING! Molecule reconstructor is not defined in the forcefield file!")
        #    else:
        #        mol.molConstruct.safety_check()

        if not self.Sampling:
            raise RuntimeError("CRITICAL ERROR! Sampling Rule has not been defined!")
        self.Sampling.safety_check()

        if self.AnalysisArray is not None:
            for analysis in self.AnalysisArray:
                analysis.safety_check()
        if self.TrajArray is not None:
            for traj in self.TrajArray:
                traj.safety_check()

        if not self.BoxList:
            raise RuntimeError("CRITICAL ERROR! No Simulation Boxes have been defined!")
        for box in self.BoxList:
            box.safety_check()

        if self.Moves:
            for move in self.Moves:
                move.safety_check()
        else:
            print("WARNING! No Monte Carlo moves have been defined!")
    # ---------------------------------------------------------------------------
    
    # ---------------------------------------------------------------------------
    def prologue(self):
        #for mol in self.MolData:
        #    if getattr(mol, 'molConstruct', None):
        #        mol.molConstruct.prologue()
        #    if hasattr(mol, 'miscdata'):
        #        for misc in mol.MiscData:
        #            misc.miscFF.prologue()

        if self.AnalysisArray is not None:
            for analysis in self.AnalysisArray:
                analysis.prologue()
        if self.Sampling is not None:
            self.Sampling.prologue()

        if self.TrajArray is not None:
            for traj in self.TrajArray:
                traj.prologue()

        for ecalc in self.EnergyCalculator:
            ecalc.prologue()
        for box in self.BoxList:
            box.prologue()
        for move in self.Moves:
            move.prologue()
    # ---------------------------------------------------------------------------
    
    # ---------------------------------------------------------------------------
    def epilogue(self):
        for box in self.BoxList:
            box.epilogue()
        print("-----------------------")

        if self.Sampling is not None:
            self.Sampling.epilogue()
        for mol in self.MolData:
            if getattr(mol, 'molConstruct', None):
                mol.epilogue()
            if hasattr(mol, 'miscdata'):
                for misc in mol.MiscData:
                    misc.epilogue()

        if self.AnalysisArray is not None:
            for analysis in self.AnalysisArray:
                analysis.epilogue()
        if self.TrajArray is not None:
            for traj in self.TrajArray:
                traj.epilogue()

        for move in self.Moves:
            move.epilogue()
    # ---------------------------------------------------------------------------
# =============================================================================
