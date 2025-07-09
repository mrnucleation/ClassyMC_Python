
import numpy as np
#========================================================================
class EasyPairCut:
    #----------------------------------------------------------------------------------
    def __init__(self):
        self.usetailcorrection = False
        self.rMin = None
        self.rMinTable = None
        self.rCut = 5.0
        self.rCutSq = 25.0
        
    #----------------------------------------------------------------------------------
    def constructor_easy_pair_cut(self):
        """Constructor equivalent"""
        self.rMinTable = 0.5  # Will need proper initialization based on atom types
        self.rCut = 5.0
        self.rCutSq = 25.0
        
    #----------------------------------------------------------------------------------
    def pair_function(self, rsq, atmtype1, atmtype2):
        """Base pair function - to be overridden in child classes"""
        return 0.0
        
    #----------------------------------------------------------------------------------
    def tail_correction(self, curbox, disp=None):
        """Calculate tail correction energy"""
        return 0.0
    #----------------------------------------------------------------------------------
    def detailed_calc(self, curbox):
        """Detailed energy calculation for all pairs"""
        atoms = curbox.get_coordinates()
        
        E_Total = 0.0
        accept = True
        
        for iAtom in range(curbox.nMaxAtoms - 1):
            atmType1 = curbox.AtomType[iAtom]
            if not curbox.is_active(iAtom):
                continue
                
            for jAtom in range(iAtom + 1, curbox.nMaxAtoms):
                if not curbox.is_active(jAtom):
                    continue
                if curbox.MolIndx[jAtom] == curbox.MolIndx[iAtom]:
                    continue
                    
                rx = atoms[0][iAtom] - atoms[0][jAtom]
                ry = atoms[1][iAtom] - atoms[1][jAtom]
                rz = atoms[2][iAtom] - atoms[2][jAtom]
                
                curbox.boundary(rx, ry, rz)
                rsq = rx**2 + ry**2 + rz**2
                
                if rsq < self.rCutSq:
                    atmType2 = curbox.AtomType[jAtom]
                    rmin_ij = self.rMinTable[atmType1][atmType2]
                    
                    if rsq < rmin_ij:
                        print(f"ERROR! Overlapping atoms found: {rsq**0.5}")
                        print(f"Atoms: {iAtom}, {jAtom}")
                        print(f"Positions: {atoms[0][iAtom]}, {atoms[1][iAtom]}, {atoms[2][iAtom]}")
                        print(f"Positions: {atoms[0][jAtom]}, {atoms[1][jAtom]}, {atoms[2][jAtom]}")
                        accept = False
                        
                    E_Pair = self.pair_function(rsq, atmType1, atmType2)
                    E_Total += E_Pair
                    curbox.ETable[iAtom] += E_Pair
                    curbox.ETable[jAtom] += E_Pair
                    
        print(f"Total Pair Energy: {E_Total}")
        
        if self.usetailcorrection:
            E_Corr = self.tail_correction(curbox)
            E_Total += E_Corr
            print(f"Total Tail Corrections: {E_Corr}")
            
        return E_Total, accept
    #----------------------------------------------------------------------------------
    def diff_calc(self, curbox, disp, tempList, tempNNei):
        """Calculate energy difference for perturbations"""
        accept = True
        curbox.dETable = [0.0] * curbox.nMaxAtoms
        E_Diff = 0.0
        
        # Handle different displacement types
        if isinstance(disp[0], Displacement):
            if disp[0].newlist:
                E_Diff, accept = self.shift_calc_single(curbox, disp, tempList, tempNNei)
            else:
                E_Diff, accept = self.shift_calc_single(curbox, disp)
                
        if self.usetailcorrection:
            E_Corr = self.tail_correction(curbox, disp)
            E_Diff += E_Corr
            
        return E_Diff, accept
    #----------------------------------------------------------------------------------
    def shift_calc_single(self, curbox, disp, tempList=None, tempNNei=None):
        """Calculate energy change for atom displacement"""
        atoms = curbox.get_coordinates()
        
        if tempList is not None and tempNNei is not None:
            neighlist = tempList
            nNeigh = tempNNei
        else:
            neighlist, nNeigh = curbox.get_neighlist()
            
        dispLen = len(disp)
        E_Diff = 0.0
        accept = True
        
        for iDisp in range(dispLen):
            if tempList is not None:
                iAtom = iDisp
            else:
                iAtom = disp[iDisp].atmIndx
                
            atmType1 = curbox.AtomType[iAtom]
            
            for jNei in range(nNeigh[iAtom]):
                jAtom = neighlist[jNei][iAtom]
                
                # New position calculation
                rx = disp[iDisp].x_new - atoms[0][jAtom]
                rx = curbox.boundary(rx)
                rsq = np.sum(rx**2)
                
                atmType2 = curbox.AtomType[jAtom]
                if rsq < self.rCutSq:
                    rmin_ij = self.rMinTable[atmType2][atmType1]
                    if rsq < rmin_ij:
                        accept = False
                        return E_Diff, accept
                        
                    E_Pair = self.pair_function(rsq, atmType1, atmType2)
                    E_Diff += E_Pair
                    
                # Old position calculation
                rx = atoms[iAtom, :] - atoms[jAtom, :]
                rx = curbox.boundary(rx)
                rsq = np.sum(rx**2)
                
                if rsq < self.rCutSq:
                    E_Pair = self.pair_function(rsq, atmType1, atmType2)
                    E_Diff -= E_Pair
        return E_Diff, accept
    #----------------------------------------------------------------------------------
    def new_calc(self, curbox, disp, tempList, tempNNei):
        """Calculate energy for new atoms being added"""
        atoms = curbox.get_coordinates()
        
        dispLen = len(disp)
        E_Diff = 0.0
        accept = True
        
        for iDisp in range(dispLen):
            iAtom = disp[iDisp].atmIndx
            atmType1 = curbox.AtomType[iAtom]
            
            listIndx = disp[iDisp].listIndex
            maxNei = tempNNei[listIndx]
            
            for jNei in range(maxNei):
                jAtom = tempList[jNei][listIndx]
                rx = disp[iDisp].x_new - atoms[0][jAtom]
                ry = disp[iDisp].y_new - atoms[1][jAtom]
                rz = disp[iDisp].z_new - atoms[2][jAtom]
                curbox.boundary(rx, ry, rz)
                rsq = rx*rx + ry*ry + rz*rz
                
                if rsq < self.rCutSq:
                    atmType2 = curbox.AtomType[jAtom]
                    rmin_ij = self.rMinTable[atmType2][atmType1]
                    
                    if rsq < rmin_ij:
                        accept = False
                        return E_Diff, accept
                        
                    E_Pair = self.pair_function(rsq, atmType1, atmType2)
                    E_Diff += E_Pair
                    curbox.dETable[iAtom] += E_Pair
                    curbox.dETable[jAtom] += E_Pair
                    
        return E_Diff, accept
        
    #----------------------------------------------------------------------------------
    def old_calc(self, curbox, disp):
        """Calculate energy for atoms being deleted"""
        neighlist, nNeigh = curbox.get_neighlist()
        atoms = curbox.get_coordinates()
        
        E_Diff = 0.0
        molStart, molEnd = curbox.get_mol_data(disp[0].molIndx)
        
        for iAtom in range(molStart, molEnd + 1):
            atmType1 = curbox.AtomType[iAtom]
            
            for jNei in range(nNeigh[iAtom]):
                jAtom = neighlist[jNei][iAtom]
                
                if iAtom == jAtom:
                    raise ValueError("Neighborlist error!")
                    
                atmType2 = curbox.AtomType[jAtom]
                rx = atoms[0][iAtom] - atoms[0][jAtom]
                ry = atoms[1][iAtom] - atoms[1][jAtom]
                rz = atoms[2][iAtom] - atoms[2][jAtom]
                curbox.boundary(rx, ry, rz)
                rsq = rx*rx + ry*ry + rz*rz
                
                if rsq < self.rCutSq:
                    E_Pair = self.pair_function(rsq, atmType1, atmType2)
                    E_Diff -= E_Pair
                    curbox.dETable[iAtom] -= E_Pair
                    curbox.dETable[jAtom] -= E_Pair
                    
        return E_Diff
        
    #----------------------------------------------------------------------------------
    def ortho_vol_calc(self, curbox, disp):
        """Calculate energy change for volume changes"""
        # Implementation would follow the Fortran logic
        # This is a complex method that would need the full context
        E_Diff = 0.0
        accept = True
        return E_Diff, accept
        
    #----------------------------------------------------------------------------------
    def atom_exchange(self, curbox, disp):
        """Calculate energy change for atom type exchange"""
        neighlist, nNeigh = curbox.get_neighlist()
        atoms = curbox.get_coordinates()
        
        E_Diff = 0.0
        accept = True
        
        iAtomNew = disp[0].newAtmIndx
        iAtomOld = disp[0].oldAtmIndx
        newType1 = curbox.AtomType[iAtomNew]
        oldType1 = curbox.AtomType[iAtomOld]
        
        for jNei in range(nNeigh[iAtomOld]):
            jAtom = neighlist[jNei][iAtomOld]
            atmType2 = curbox.AtomType[jAtom]
            
            rx = curbox.atoms[0][iAtomOld] - curbox.atoms[0][jAtom]
            ry = curbox.atoms[1][iAtomOld] - curbox.atoms[1][jAtom]
            rz = curbox.atoms[2][iAtomOld] - curbox.atoms[2][jAtom]
            curbox.boundary(rx, ry, rz)
            rsq = rx*rx + ry*ry + rz*rz
            
            if rsq < self.rCutSq:
                rmin_ij = self.rMinTable[atmType2][newType1]
                if rsq < rmin_ij:
                    accept = False
                    return E_Diff, accept
                    
                E_Pair_new = self.pair_function(rsq, newType1, atmType2)
                E_Pair_Old = self.pair_function(rsq, oldType1, atmType2)
                
                E_Diff += E_Pair_new - E_Pair_Old
                curbox.dETable[iAtomOld] += E_Pair_new - E_Pair_Old
                curbox.dETable[jAtom] += E_Pair_new - E_Pair_Old
                
        return E_Diff, accept
        
    #----------------------------------------------------------------------------------
    def get_cutoff(self):
        """Return the cutoff distance"""
        return self.rCut
    #----------------------------------------------------------------------------------
    def process_io(self, line):
        """Process input/output line"""
        pass
    #----------------------------------------------------------------------------------
#========================================================================