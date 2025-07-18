!========================================================
module MCMove_MolTranslation
use CoordinateTypes, only: Displacement
use MoveClassDef
use SimpleSimBox, only: SimpleBox
use VarPrecision

  type, public, extends(MCMove) :: MolTranslate
!    real(dp) :: atmps = 1E-30_dp
!    real(dp) :: accpt = 0E0_dp

    real(dp), allocatable :: boxatmps(:)
    real(dp) , allocatable:: boxaccpt(:)
    logical :: verbose = .true.
    logical :: proportional = .true.
    logical :: tuneMax = .true.
    real(dp) :: limit = 5.00E0_dp
    real(dp) :: targAccpt = 50E0_dp
    real(dp) :: max_dist = 0.05E0_dp
    real(dp), allocatable :: boxlimit(:)
    real(dp), allocatable :: boxtargAccpt(:)
    real(dp), allocatable :: boxmax_dist(:)

    type(Displacement), allocatable :: disp(:)

    !Rejection Counters
    integer :: ovlaprej = 0 
    integer :: constrainrej = 0 
    integer :: detailedrej = 0 

!    integer, allocatable :: tempNnei(:)
!    integer, allocatable :: tempList(:, :)

    contains
      procedure, pass :: Constructor => MolTrans_Constructor
!      procedure, pass :: GeneratePosition => MolTrans_GeneratePosition
      procedure, pass :: FullMove => MolTrans_FullMove
      procedure, pass :: Maintenance => MolTrans_Maintenance
      procedure, pass :: Prologue => MolTrans_Prologue
      procedure, pass :: Update => MolTrans_Update
      procedure, pass :: Epilogue => MolTrans_Epilogue
      procedure, pass :: ProcessIO => MolTrans_ProcessIO
  end type
!========================================================
 contains
!========================================================
  subroutine MolTrans_Constructor(self)
    use Common_MolInfo, only: MolData, nMolTypes
    use BoxData, only: BoxArray
    implicit none
    class(MolTranslate), intent(inout) :: self
    integer :: iType, maxAtoms, nBoxes

    nBoxes = size(boxArray)
    if(.not. allocated(self%boxProb)) then
      allocate( self%boxProb(1:nBoxes) )
      self%boxProb = 1E0_dp/real(nBoxes,dp)
    endif

    allocate( self%boxatmps(1:nBoxes) )
    self%boxatmps = 1e-50_dp
    allocate( self%boxaccpt(1:nBoxes) )
    self%boxaccpt = 0E0_dp

    allocate( self%boxLimit(1:nBoxes) )
    self%boxLimit = self%limit
    allocate( self%boxmax_dist(1:nBoxes) )
    self%boxmax_dist = self%max_dist
    allocate( self%boxtargAccpt(1:nBoxes) )
    self%boxtargAccpt = self%targAccpt

    maxAtoms = 0
    do iType = 1, nMolTypes
      if(MolData(iType)%nAtoms > maxAtoms) then
        maxAtoms = MolData(iType)%nAtoms 
      endif
    enddo


    allocate( self%disp(1:maxAtoms) )

    call self%CreateTempArray(maxAtoms)
  end subroutine
!========================================================
!  subroutine MolTrans_GeneratePosition(self, disp)
!    use RandomGen, only: grnd
!    implicit none
!    class(MolTranslate), intent(in) :: self
!    type(Displacement), intent(inout) :: disp
!    real(dp) :: dx, dy, dz
!      dx = self % max_dist * (2E0_dp*grnd() - 1E0_dp)
!      dy = self % max_dist * (2E0_dp*grnd() - 1E0_dp)
!      dz = self % max_dist * (2E0_dp*grnd() - 1E0_dp)
!  end subroutine
!===============================================
  subroutine MolTrans_FullMove(self, trialBox, accept) 
    use Box_Utility, only: FindAtom, FindMolecule
    use CommonSampling, only: sampling
    use Common_NeighData, only: neighSkin
    use Common_MolInfo, only: MolData, nMolTypes
    use ForcefieldData, only: EnergyCalculator
    use RandomGen, only: grnd
    implicit none
    class(MolTranslate), intent(inout) :: self
    class(SimpleBox), intent(inout) :: trialBox
    logical, intent(out) :: accept
    integer :: boxID, iAtom, nAtoms, atomIndx
    integer :: nMove, rawIndx, iConstrain
    integer :: CalcIndex, molStart, molEnd, molType
    real(dp) :: dx, dy, dz
    real(dp) :: E_Diff, E_Inter, E_Intra, biasE
    real(dp), parameter :: Prob = 1E0_dp

    boxID = trialBox % boxID
    self % atmps = self % atmps + 1E0_dp
    self % boxatmps(boxID) = self % boxatmps(boxID) + 1E0_dp
    accept = .true.
    call self%LoadBoxInfo(trialBox, self%disp)
    !Propose move
    rawIndx = floor( trialBox%nMolTotal * grnd() + 1E0_dp )
    call FindMolecule(trialbox, rawIndx, nMove)
    call trialBox % GetMolData(nMove, molStart=molStart, molEnd=molEnd, &
                               molType=molType)

    dx = self % boxmax_dist(boxID) * (2E0_dp * grnd() - 1E0_dp)
    dy = self % boxmax_dist(boxID) * (2E0_dp * grnd() - 1E0_dp)
    dz = self % boxmax_dist(boxID) * (2E0_dp * grnd() - 1E0_dp)
 
    nAtoms = MolData(molType)%nAtoms
    do iAtom = 1, nAtoms
      atomIndx = molStart + iAtom - 1

      self%disp(iAtom)%molType = molType
      self%disp(iAtom)%molIndx = nMove
      self%disp(iAtom)%atmIndx = atomIndx

      self%disp(iAtom)%x_new = trialBox%atoms(1, atomIndx) + dx
      self%disp(iAtom)%y_new = trialBox%atoms(2, atomIndx) + dy
      self%disp(iAtom)%z_new = trialBox%atoms(3, atomIndx) + dz

      self%disp(iAtom)%newlist = .false.
      self%disp(iAtom)%listIndex = iAtom
    enddo
!    write(*,*) dx, dy, dz

    !If the particle moved a large distance get a temporary neighborlist
!    if(any([dx,dy,dz] > neighSkin)) then
!      call trialBox % NeighList(1) % GetNewList(1, self%tempList, self%tempNNei, self%disp(1))
!      self%disp(1)%newlist = .true.
!    else

!    endif

    !Check Constraint
    accept = trialBox % CheckConstraint( self%disp(1:nAtoms) )
    if(.not. accept) then
      self%constrainrej = self%constrainrej + 1
      return
    endif

    !Energy Calculation
!    call trialbox% EFunc % Method % DiffECalc(trialBox, self%disp(1:nAtoms), self%tempList, self%tempNNei, E_Diff, accept)
    call trialBox%ComputeEnergyDelta(self%disp(1:nAtoms),&
                                     self%templist,&
                                     self%tempNNei, &
                                     E_Inter, &
                                     E_Intra, &
                                     E_Diff, &
                                     accept, &
                                     computeintra=.false.)
    if(.not. accept) then
      self%ovlaprej = self%ovlaprej + 1
      return
    endif

    !Check Post Energy Constraint
    accept = trialBox % CheckPostEnergy( self%disp(1:nAtoms), E_Diff )
    if(.not. accept) then
      self%constrainrej = self%constrainrej + 1
      return
    endif




    !Accept/Reject
    accept = sampling % MakeDecision(trialBox, E_Diff, self%disp(1:nAtoms), inProb=Prob)
    if(accept) then
      self % accpt = self % accpt + 1E0_dp
      self % boxaccpt(boxID) = self % boxaccpt(boxID) + 1E0_dp
      call trialBox % UpdateEnergy(E_Diff, E_Inter)
      call trialBox % UpdatePosition(self%disp(1:nAtoms), self%tempList, self%tempNNei)
    else
      self%detailedrej = self%detailedrej + 1
!      write(*,*) E_Diff, trialBox%beta, Prob
    endif

  end subroutine
!=========================================================================
  subroutine MolTrans_Maintenance(self)
    implicit none
    class(MolTranslate), intent(inout) :: self
    integer :: iBox
    real(dp) :: accRate
!    real(dp), parameter :: lowerlimit = 0.1E0_dp
      
    if(self%tuneMax) then
      do iBox = 1, size(self%boxatmps)
        if(self%boxatmps(iBox) < 0.5E0_dp) then
          cycle
        endif
        accRate = 1e2_dp*self%boxaccpt(iBox)/self%boxatmps(iBox)

        if(accRate > self%boxtargAccpt(iBox)) then
          if(self%boxmax_dist(iBox)*1.01E0_dp < self%boxlimit(iBox)) then
            self%boxmax_dist(iBox) = self%boxmax_dist(iBox) * 1.01E0_dp
          else 
            self%boxmax_dist(iBox) = self%boxlimit(iBox)
          endif
        else
          self%boxmax_dist(iBox) = self%boxmax_dist(iBox) * 0.99E0_dp
        endif
      enddo
    endif

  end subroutine
!=========================================================================
  subroutine MolTrans_Prologue(self)
    use ParallelVar, only: nout
    implicit none
    class(MolTranslate), intent(inout) :: self

    if(.not. allocated(self%disp)) then
      call self % Constructor
    endif
      

    write(nout,"(1x,A,F15.8)") "(Molecule Translate) Maximum Displacement: ", self%max_dist

  end subroutine
!=========================================================================
  subroutine MolTrans_Epilogue(self)
    use ParallelVar, only: nout
    implicit none
    class(MolTranslate), intent(inout) :: self
    real(dp) :: accptRate
      
    write(nout,*) 
    write(nout,"(1x,A,I15)") "Molecule Translation Moves Accepted: ", nint(self%accpt)
    write(nout,"(1x,A,I15)") "Molecule Translation Moves Attempted: ", nint(self%atmps)
    accptRate = self%GetAcceptRate()
    write(nout,"(1x,A,F15.8)") "Molecule Translation Acceptance Rate: ", accptRate
    if(self%tunemax) then
      write(nout,"(1x,A,100F15.8)") "Final Maximum Displacement: ", self%boxmax_dist(1:)
    endif

    if(self%verbose) then
      write(nout, "(1x,A,I15)") "Molecule Translation, Rejections due to overlap:", self%ovlaprej
      write(nout, "(1x,A,I15)") "Molecule Translation, Rejections due to constraint:", self%constrainrej
      write(nout, "(1x,A,I15)") "Molecule Translation, Rejections due to detailed balance:", self%detailedrej
    endif

  end subroutine
!=========================================================================
  subroutine MolTrans_Update(self)
    use BoxData, only: BoxArray
    use ParallelVar, only: nout
    implicit none
    class(MolTranslate), intent(inout) :: self
    integer :: iBox
    real(dp) :: norm

      if(self%proportional) then
        norm = 0E0_dp
        do iBox = 1, size(BoxArray)
          norm = norm + real(BoxArray(ibox)%box%nMolTotal,dp)
        enddo
        do iBox = 1, size(BoxArray)
          self%boxprob(iBox) = real(BoxArray(ibox)%box%nMolTotal,dp)/norm
        enddo
      endif
  end subroutine
!=========================================================================
  subroutine MolTrans_ProcessIO(self, line, lineStat)
    use Input_Format, only: GetXCommand, maxLineLen
    implicit none
    class(MolTranslate), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line
    integer, intent(out) :: lineStat
    character(len=30) :: command
    logical :: logicVal
    integer :: intVal
    real(dp) :: realVal

    call GetXCommand(line, command, 4, lineStat)
    select case( trim(adjustl(command)) )
      case("tunemax")
        call GetXCommand(line, command, 5, lineStat)
        read(command, *) logicVal
        self%tunemax = logicVal

      case("dynamiclimit")
        call GetXCommand(line, command, 5, lineStat)
        read(command, *) realVal
        self%limit = realVal

      case("dynamictarget")
        call GetXCommand(line, command, 5, lineStat)
        read(command, *) realVal
        self%targAccpt = realVal

      case("maxdisplace")
        call GetXCommand(line, command, 5, lineStat)
        read(command, *) realVal
        self%max_dist = realVal

      case("proportional")
        call GetXCommand(line, command, 5, lineStat)
        read(command, *) logicVal
        self%proportional = logicVal


      case("updatefreq")
        call GetXCommand(line, command, 5, lineStat)
        read(command, *) intVal
        self%maintFreq = intVal

      case default
        lineStat = -1
        return

    end select
    lineStat = 0

  end subroutine
!========================================================
end module
!========================================================
