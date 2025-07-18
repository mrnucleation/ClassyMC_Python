!===================================================================================
! This module contains a simple neighborlist
!===================================================================================
module RSqListDef
use VarPrecision
use CoordinateTypes
use Template_SimBox, only: SimBox
use SimpleSimBox, only: SimpleBox
use Template_NeighList, only: NeighListDef

  type, public, extends(NeighListDef) :: RSqList
!      logical :: Sorted = .false.
!      logical :: Strict = .false.
!      integer, allocatable :: list(:,:)
!      integer, allocatable :: nNeigh(:)
!      integer :: maxNei
!      real(dp) :: rCut, rCutSq
!      logical :: restrictType = .false.
!      integer, allocatable :: allowed(:)
!      integer :: safetyCheck = .false.

      class(SimpleBox), pointer :: parent => null()
    contains
      procedure, pass :: Constructor => RSqList_Constructor 
      procedure, pass :: BuildList => RSqList_BuildList 
      procedure, pass :: GetNewList => RSqList_GetNewList
      procedure, pass :: AddMol => RSqList_AddMol
      procedure, pass :: SwapAtomType => RSqList_SwapAtomType
      procedure, pass :: GetNeighCount => RSqList_GetNeighCount
      procedure, pass :: ProcessIO => RSqList_ProcessIO
!      procedure, pass :: TransferList
      procedure, pass :: DeleteMol => RSqList_DeleteMol
      procedure, pass :: Prologue => RSqList_Prologue
      procedure, pass :: Update => RSqList_Update
  end type

!===================================================================================
  contains
!===================================================================================
  subroutine RSqList_Constructor(self, parentID, rCut)
    use BoxData, only: BoxArray
    use Common_NeighData, only: neighSkin
    use Common_MolInfo, only: nAtomTypes, nMolTypes, MolData
    use ParallelVar, only: nout
    implicit none
    class(RSqList), intent(inout) :: self
    integer, intent(in) :: parentID
    real(dp), intent(in), optional :: rCut
    real(dp), parameter :: atomRadius = 0.85E0_dp  !Used to estimate an approximate volume of 
    integer :: iType
    integer :: AllocateStatus
    integer :: maxAtoms

    maxAtoms = 0
    do iType = 1, nMolTypes
      if(MolData(iType)%nAtoms > maxAtoms) then
        maxAtoms = MolData(iType)%nAtoms 
      endif
    enddo

    self%parent => BoxArray(parentID)%box
    if(.not. allocated(self%parent%atoms) ) then
      error stop
    endif

!     If no rCut value is given by the subroutine call attempt to pull
!     the rSq value from the parent box's energy function. The assumption being
!     that the neighborlist is used for the energy calculation routines.
    if( present(rCut) ) then
      self % rCut = rCut
      self % rCutSq = rCut * rCut
      self % maxNei = ceiling(rCut**3/atomRadius**3)
    else
      if(self%rCut > 0E0_dp) then
        self % rCutSq = (self%rCut)**2
        self % maxNei = ceiling(self%rCut**3/atomRadius**3)

      else
        self % rCut = self % parent % EFunc % Method % GetCutOff() + neighSkin
        self % rCutSq = (self%rCut)**2
        self % maxNei = ceiling(self%rCut**3/atomRadius**3)
      endif
    endif
 
    write(nout,*) "Neighbor List CutOff:", self%rCut
    if(self%maxNei > self%parent%nMaxAtoms) then
      self%maxNei = self%parent%nMaxAtoms
    endif
    write(nout,*) "Neighbor List Maximum Neighbors:", self%maxNei

    allocate( self%list(1:self%maxNei, 1:self%parent%nMaxAtoms), stat=AllocateStatus )
    allocate( self%nNeigh(1:self%parent%nMaxAtoms), stat=AllocateStatus )

    allocate( self%templist(1:self%maxNei+1, 1:maxAtoms), stat=AllocateStatus )
    allocate( self%tempNNeigh(1:maxAtoms), stat=AllocateStatus )

    if(.not. allocated(self%allowed) ) then
      allocate(self%allowed(1:nAtomTypes), stat=AllocateStatus )
      self%allowed = .true.
    endif

    self%list = 0
    self%nNeigh = 0 
    IF (AllocateStatus /= 0) error STOP "*** NeighRSQList: Not enough memory ***"

    self%restrictType = .false.
  end subroutine
!===================================================================================
  subroutine RSqList_Prologue(self)
    implicit none
    class(RSqList), intent(inout) :: self


!    call self%DumpList(2)
  end subroutine
!===================================================================================
  subroutine RSqList_Update(self)
    implicit none
    class(RSqList), intent(inout) :: self


!    call self%DumpList(2)
  end subroutine
!===================================================================================
  subroutine RSqList_BuildList(self, listindx)
    implicit none
    class(RSqList), intent(inout) :: self
    integer, intent(in) :: listindx

    if(listindx > 1) then
      return
    endif


    call Builder_RSq(self%parent)
  end subroutine
!===================================================================================
  function RSqList_GetNeighCount(self, nAtom, rCount) result(nCount)
    implicit none
    class(RSqList), intent(inout) :: self
    integer, intent(in) :: nAtom
    real(dp), intent(in) ,optional :: rCount
    integer :: iNei
    real(dp) :: rCut, rCutSq
    real(dp) :: rx, ry, rz, rsq
    integer :: nCount, jAtom

    if(present(rCount)) then
      rCut = rCount
      rCutSq = rCount*rCount
    else
      rCut = self%rCut
      rCutSq = self%rCutSq
    endif

    nCount = 0
    do iNei = 1, self % nNeigh(nAtom)
      jAtom = self%list(iNei, nAtom) 
      rx = self%parent%atoms(1, nAtom) - self%parent%atoms(1, jAtom)
      ry = self%parent%atoms(2, nAtom) - self%parent%atoms(2, jAtom)
      rz = self%parent%atoms(3, nAtom) - self%parent%atoms(3, jAtom)
      call self%parent%Boundary(rx,ry,rz)
      rsq = rx*rx + ry*ry + rz*rz
      if(rsq < rCutSq) then
        nCount = nCount + 1
      endif
    enddo

  end function
!===================================================================================
  subroutine RSqList_AddMol(self, disp, tempList, tempNNei)
    implicit none
    class(RSqList), intent(inout) :: self
    class(Perturbation), intent(in) :: disp(:)
    integer, intent(in) :: tempList(:,:), tempNNei(:)

    select type(disp)

      class is(Addition)
        call UpdateList_AddMol_RSq(self%parent, disp, tempList, tempNNei)
    end select

  end subroutine
!===================================================================================
  subroutine RSqList_DeleteMol(self, molIndx, topIndx)
    use Common_MolInfo, only: nMolTypes, MolData
    use SearchSort, only: BinarySearch, SimpleSearch
    implicit none
    class(RSqList), intent(inout) :: self
    integer, intent(in) :: molIndx, topIndx
    integer :: iAtom, iNei, jNei, nType, j
    integer :: nStart, topStart
    integer :: nEnd, topEnd
    integer :: atmIndx, topAtom
    integer :: curNei, curIndx, nNei

    nStart = self % parent % MolStartIndx(molIndx)
    topStart = self % parent % MolStartIndx(topIndx)

    nEnd = self % parent % MolEndIndx(molIndx)
    topEnd = self % parent % MolEndIndx(topIndx)
    nType = self % parent % MolType(nStart)

!    write(2,*) "----------------------------"
!    write(2,*) "Delete"
!    write(2,*) "Removed Mol:", molIndx
!    do iAtom = 1, self%parent%nMaxatoms
!      write(2,"(I3,A,1000(I3))") iAtom,"|", (self%list(j, iAtom) ,j=1,self%nNeigh(iAtom))
!    enddo
!    write(2,*) "Sorted?:", self%sorted

    do iAtom = 1, MolData(nType)%nAtoms
!      atmIndx = nStart + iAtom - 1
!      topAtom = topStart + iAtom - 1


      atmIndx = nEnd - iAtom + 1
      topAtom = topEnd - iAtom + 1
!      write(2,*) "iAtom", atmIndx, "topAtom", topAtom
      !Remove the deleted from the list of it's neighbors
      do iNei = 1, self % nNeigh(atmIndx)
        curNei = self % list(iNei, atmIndx)
        nNei = self%nNeigh(curNei)
        if(nNei == 0) then
          cycle
        endif

        if(self%sorted) then
          curIndx = BinarySearch( atmIndx, self%list(1:nNei, curNei) )
        else
          curIndx = SimpleSearch( atmIndx, self%list(1:nNei, curNei) )
        endif

!        if(nNei <= 1) then
!          self%nNeigh(curNei) = 0
!          self%list(:, curNei) = 0
!          curIndx = 0
!          cycle
!        else
!
!        endif
!        write(2,*) "curNei", curNei, "Indexes", curIndx, atmIndx
        if(curIndx /= 0) then
          if(nNei > 2) then
            self%list(1:nNei, curNei ) = [self%list(1:curIndx-1, curNei), &
                                          self%list(curIndx+1:nNei, curNei) ]
          else
            if(curIndx == 1) then
              self%list(1, curNei) = self%list(2,curNei)
            endif
          endif
          self%nNeigh(curNei) = self%nNeigh(curNei) - 1 

        endif
      enddo
    enddo
!    do iAtom = 1, self%parent%nMaxatoms
!      write(2,"(I3,A,1000(I3))") iAtom,"|", (self%list(j, iAtom) ,j=1,self%nNeigh(iAtom))
!    enddo

    do iAtom = 1, MolData(nType)%nAtoms
      atmIndx = nEnd - iAtom + 1
      topAtom = topEnd - iAtom + 1     
      !Move the top atom into the spot formerly taken up by the deleted atom.
      do iNei = 1, self % nNeigh(topAtom)
        self%list(iNei, atmIndx) = self%list(iNei, topAtom)
      enddo
      self%nNeigh(atmIndx) = self%nNeigh(topAtom)

      !Re-index the neighbor's of the top atom to it's new array location
      do iNei = 1, self % nNeigh(topAtom)
        curNei = self % list(iNei, topAtom)
        nNei = self%nNeigh(curNei)
        if(nNei == 0 ) then
          cycle
        endif
        if(self%sorted) then
          curIndx = BinarySearch( topAtom, self%list(1:nNei, curNei) )
          self%sorted = .false.
        else
          curIndx = SimpleSearch( topAtom, self%list(1:nNei, curNei) )
        endif
!        write(2,*) "Second", "curNei", curNei, "Indexes", curIndx, topAtom
        if(curIndx /= 0) then
          self % list(curIndx, curNei) = atmIndx
        endif
      enddo
      self%nNeigh(topAtom) = 0
!      do iNei = 1, self%parent%nMaxatoms
!        write(2,"(I3,A,1000(I3))") iNei,"|", (self%list(j, iNei) ,j=1,self%nNeigh(iNei))
!      enddo
    enddo

!    do iAtom = 1, self%parent%nMaxatoms
!      write(2,"(I3,A,1000(I3))") iAtom,"|", (self%list(j, iAtom) ,j=1,self%nNeigh(iAtom))
!    enddo
!    write(2,*) "---------------------------------"


    self % sorted = .false.

  end subroutine
!===================================================================================
  subroutine RSqList_SwapAtomType(self, disp, topIndx)
    !----------------
    ! Routine used to change the type of an atom.  Primarily used for atom swap moves.
    ! Not currently designed to work with molecules. 
    ! Variables
    !    input
    !        disp => Displacement class variable which contains information related to
    !                what changed in the system.  For example this might contain the old/new
    !                positions of an atom that was shifted. This routine expects an
    !                AtomExchange class which contains the array location of 
    !                the atom being changed and it's new array location.
    !    
    !    function variables
    !       iAtomNew => Index of the atom's new position in the array
    !       iAtomOld => Index of the atom's old position in the array
    !       iNei => Loop integer for looping over the number of neighbors in the neighborlist list
    !       cellIndx => Cell ID of the current atom. 
    !       curIndx => Atom Index returned by the search algorithm
    !
    !---------
    use Common_MolInfo, only: nMolTypes, MolData
    use SearchSort, only: SimpleSearch
    implicit none
    class(RSqList), intent(inout) :: self
    class(AtomExchange), intent(in) :: disp(:)
    integer, intent(in) :: topIndx
    integer :: iAtomNew, iAtomOld, iNei, jAtom, iAtom, nNei
    integer :: topAtom
    integer :: cellIndx, curIndx
    integer :: molIndxNew, molIndxOld



    iAtomNew = disp(1)%newAtmIndx
    iAtomOld = disp(1)%oldAtmIndx
    topAtom = self%parent%MolStartIndx(topIndx)

    molIndxNew = self%parent%MolIndx(iAtomNew)
    molIndxOld = self%parent%MolIndx(iAtomOld)

    !Search through the neighborlists and replace the old atom index for
    !a new atom index
    nNei = self%nNeigh(iAtomOld)
    do iNei = 1, nNei
      jAtom = self % list(iNei, iAtomOld)
      curIndx = SimpleSearch( iAtomOld, self%list(1:nNei, jAtom) )
      self%list(curIndx, jAtom) = iAtomNew
    enddo


    nNei = self%nNeigh(topAtom)
    do iNei = 1, nNei
      jAtom = self % list(iNei, topAtom)
      if(jAtom == iAtomNew) then
        cycle
      endif
      curIndx = SimpleSearch( topAtom, self%list(1:nNei, jAtom) )
      if(curIndx /= 0) then
        self%list(curIndx, jAtom) = iAtomOld
      endif
    enddo
    curIndx = SimpleSearch( topAtom, self%list(1:nNei, iAtomOld) )
    if(curIndx /= 0) then
      self%list(curIndx, iAtomOld) = iAtomOld
    endif

    !Copy the neighbor of the atom being swapped list from the old location
    !to the new location
    do iNei = 1, self%nNeigh(iAtomOld)
      self % list(iNei, iAtomNew) =  self % list(iNei, iAtomOld)
    enddo
    do iNei = 1, self%nNeigh(topAtom)
      self % list(iNei, iAtomOld) =  self % list(iNei, topAtom)
    enddo


    self%nNeigh(iAtomNew) = self%nNeigh(iAtomOld) 
    self%nNeigh(iAtomOld) = self%nNeigh(topAtom) 
    self%nNeigh(topAtom) = 0
   
    self%sorted = .false.



  end subroutine

!===================================================================================
  subroutine RSqList_GetNewList(self, iDisp, tempList, tempNNei, disp, nCount, rCount)
    use Common_MolInfo, only: nMolTypes
    implicit none
    class(RSqList), intent(inout) :: self
    integer, intent(in) :: iDisp
    class(Perturbation), intent(inout) :: disp
    integer, intent(inout) :: tempList(:,:), tempNNei(:)
    integer, optional :: nCount
    real(dp), optional :: rCount
    integer :: jType, jAtom, j, iAtom
    integer :: jUp, jLow, molIndx, jMol
    real(dp) :: xn, yn, zn
    real(dp) :: rx, ry, rz, rsq

!    if(.not. (associated(tempList) .and. associated(tempNNei)))then
!      error stop "Unassociated Temporary List passed into CellRSq_GetNewList"
!    endif

    if(present(nCount)) then
      nCount = 0
    endif
    select type(disp)
      class is (Addition)
!        disp % newlist = .true.
        disp % listIndex = iDisp
        iAtom = disp%atmIndx
        molIndx = self%parent%MolIndx(iAtom)
        xn = disp%x_new
        yn = disp%y_new
        zn = disp%z_new
      class is (Displacement)
        disp % newlist = .true.
        disp % listIndex = iDisp
        iAtom = disp%atmIndx
        molIndx = self%parent%MolIndx(iAtom)
        xn = disp%x_new
        yn = disp%y_new
        zn = disp%z_new
    end select

    templist(:, iDisp) = 0
    tempNNei(iDisp) = 0

!    molStart = 1
!    do jType = 1, nMolTypes
      do jAtom = 1, self%parent%nMaxAtoms
        if( self%parent%MolSubIndx(jAtom) == molIndx ) then
          cycle
        endif
        if( self%parent%MolSubIndx(jAtom) > self%parent%NMol(self%parent%MolType(jAtom)) ) then
          cycle
        endif
        rx = xn - self%parent%atoms(1, jAtom)
        ry = yn - self%parent%atoms(2, jAtom)
        rz = zn - self%parent%atoms(3, jAtom)
        call self%parent%Boundary(rx,ry,rz)
        rsq = rx*rx + ry*ry + rz*rz
        if(rsq < self%rCutSq) then
          tempNNei(iDisp) = tempNNei(iDisp) + 1
          templist(tempNNei(iDisp), iDisp) = jAtom
        endif
        if(present(rCount)) then
          if(rsq < rCount*rCount) then
            nCount = nCount + 1
          endif
        endif
      enddo
!      molStart = molStart + self%parent%NMolMax(jType)
!    enddo
!    write(2,"(A, 1000(I3))") "New",   templist(1:tempNNei(iDisp), iDisp)
  end subroutine
!====================================================================
  subroutine RSqList_ProcessIO(self, line, lineStat)
    use Input_Format, only: GetAllCommands, GetXCommand,maxLineLen
    use Common_MolInfo, only: nAtomTypes
    implicit none
    class(RSqList), intent(inout) :: self
    integer, intent(out) :: lineStat
    character(len=maxLineLen), intent(in) :: line   

    integer :: i, intVal, nPar
    real(dp) :: realVal

    character(len=30) :: command 
    character(len=30), allocatable :: parlist(:)


    lineStat = 0
    call GetXCommand(line, command, 6, lineStat)
    select case( trim(adjustl(command)) )
      case("rcut")
        call GetXCommand(line, command, 7, lineStat)
        read(command,*) realVal
        self%rCut = realVal
        self%rCutSq = realVal * realVal

      case("restricttype")
        call GetAllCommands(line, parlist, nPar, lineStat)
        self%restrictType = .true.
        if(.not. allocated(self%allowed) ) then
          allocate(self%allowed(1:nAtomTypes) )
        endif
        self%allowed = .false.
        do i = 7, size(parList)
          read(parList(i), *) intVal
          self%allowed(intVal) = .true.
        enddo

      case default
        lineStat = -1
    end select

  end subroutine
!===================================================================================
! End Type Bound
!===================================================================================
  subroutine Builder_RSq(trialBox)
    use Common_MolInfo, only: nMolTypes, MolData
    use Common_NeighData, only: neighSkin
    use ParallelVar, only: nout
    implicit none
    class(SimpleBox), intent(inout) :: trialBox
    integer :: iList
    integer :: iType, jType, iAtom, jAtom, j
    integer :: iUp, iLow, jUp, jLow, molStart, jMolStart, jMolEnd, atmType
    integer :: nNeigh
!    integer, allocatable :: oldlist(:)
    real(dp) :: rx, ry, rz, rsq

!    if(trialBox%NeighList(1)%safetyCheck) then
!      allocate(oldlist(1:trialBix%NeighList(1)%maxNei)
!    endif

!    write(*,*) "Building"
    do iList = 1, size(trialBox%NeighList)
      trialBox%NeighList(iList)%nNeigh = 0
      trialBox%NeighList(iList)%list = 0
    enddo

    do iAtom = 1, trialBox%nMaxAtoms-1
      if( trialBox%MolSubIndx(iAtom) > trialBox%NMol(trialBox%MolType(iAtom)) ) then
        cycle
      endif
      do jAtom = iAtom+1, trialBox%nMaxAtoms
        if( trialBox%MolSubIndx(jAtom) > trialBox%NMol(trialBox%MolType(jAtom)) ) then
          cycle
        endif

        if( trialBox%MolIndx(iAtom) == trialBox%MolIndx(jAtom)) then
          cycle
        endif

        rx = trialBox%atoms(1, iAtom) - trialBox%atoms(1, jAtom)
        ry = trialBox%atoms(2, iAtom) - trialBox%atoms(2, jAtom)
        rz = trialBox%atoms(3, iAtom) - trialBox%atoms(3, jAtom)
        call trialBox%Boundary(rx, ry, rz)
        rsq = rx*rx + ry*ry + rz*rz
        do iList = 1, size(trialBox%NeighList)
          if( trialBox % NeighList(iList) % restrictType ) then
            atmType = trialBox % atomType(iAtom)
            if( trialBox%NeighList(iList)%allowed(atmType)  ) then
              cycle
            endif

            atmType = trialBox % atomType(jAtom)
            if( trialBox%NeighList(iList)%allowed(atmType)  ) then
              cycle
            endif
          endif
          if( rsq <= trialBox%NeighList(iList)%rCutSq ) then 

            trialBox%NeighList(iList)%nNeigh(iAtom) = trialBox%NeighList(iList)%nNeigh(iAtom) + 1
            if(trialBox%NeighList(iList)%nNeigh(iAtom) > trialBox%NeighList(iList)%maxNei) then
              write(nout, *) "Neighborlist overflow!"
            endif
            trialBox%NeighList(iList)%list( trialBox%NeighList(iList)%nNeigh(iAtom), iAtom ) = jAtom

            trialBox%NeighList(iList)%nNeigh(jAtom) = trialBox%NeighList(iList)%nNeigh(jAtom) + 1
            if(trialBox%NeighList(iList)%nNeigh(jAtom) > trialBox%NeighList(iList)%maxNei) then
              write(nout, *) "Neighborlist overflow!"
            endif
            trialBox%NeighList(iList)%list( trialBox%NeighList(iList)%nNeigh(jAtom), jAtom ) = iAtom
          endif
        enddo
      enddo  
    enddo

!    write(2,*) "----------------------------"
!    do iAtom = 1, trialBox%nMaxAtoms
!      if( trialBox%MolSubIndx(iAtom) > trialBox%NMol(trialBox%MolType(iAtom)) ) then
!        cycle
!      endif
!      nNeigh = trialBox%NeighList(1)%nNeigh(iAtom)
!      write(2,"(I6,A,999(I4))") iAtom, "|", trialBox%NeighList(1)%list(1:nNeigh, iAtom)
!    enddo


    do iList = 1, size(trialBox%NeighList)
      trialBox%NeighList(iList)%sorted = .true.
    enddo
  end subroutine
!===================================================================================
  subroutine UpdateList_AddMol_RSq(trialBox, disp, tempList, tempNNei)
    use Common_MolInfo, only: nMolTypes, MolData
    implicit none
    class(SimpleBox), intent(inout) :: trialBox
    class(Addition), intent(in) :: disp(:)
    integer, intent(in) :: tempList(:,:), tempNNei(:)
    integer :: iList, iDisp, iAtom, iNei, nNei, neiIndx, j
    real(dp) :: rx, ry, rz, rsq

    do iList = 1, size(trialBox%NeighList)
      if(iList == 1) then
!        write(2,*) "----------------------------"
!        write(2,*) "Add"
!        do iAtom = 1, trialBox%nMaxatoms
!          write(2,"(I3,A,1000(I3))") iAtom,"|", (trialBox % NeighList(iList)%list(j, iAtom) ,j=1,trialBox % NeighList(iList)%nNeigh(iAtom))
!1        enddo
!        write(2,*)
        do iDisp = 1, size(disp)
!          write(2,"(A,A,1000(I3))") "NewList","|", (tempList(j, iDisp) ,j=1,tempNNei(iDisp))
          iAtom = disp(iDisp)%atmIndx
          trialBox % NeighList(iList) % nNeigh(iAtom) = tempNNei(iDisp)
          do iNei = 1, tempNNei(iDisp)
            neiIndx = tempList(iNei, iDisp)
            trialBox % NeighList(iList) % list(iNei, iAtom) =  neiIndx
            trialBox % NeighList(iList) % list( trialBox%NeighList(iList)%nNeigh(neiIndx)+1, neiIndx ) = iAtom
            trialBox%NeighList(iList)%nNeigh(neiIndx)= trialBox%NeighList(iList)%nNeigh(neiIndx) + 1
          enddo
        enddo
!        do iAtom = 1, trialBox%nMaxatoms
!          write(2,"(I3,A,1000(I3))") iAtom,"|", (trialBox % NeighList(iList)%list(j, iAtom) ,j=1,trialBox % NeighList(iList)%nNeigh(iAtom))
!        enddo
!        write(2,*)
      endif
!      write(2,*) "N", trialBox%NeighList(iList)%nNeigh(:)     
    enddo



  end subroutine
!===================================================================================
end module
!===================================================================================
