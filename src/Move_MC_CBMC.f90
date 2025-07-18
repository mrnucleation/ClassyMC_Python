!========================================================
! Move which cuts and regrows a molecule in conjunction with
! the CBMC MolCon Regrowth Functions.  
!========================================================
module MCMove_CBMC
use CoordinateTypes, only: Displacement
use MoveClassDef
use SimpleSimBox, only: SimpleBox
use VarPrecision

  type, public, extends(MCMove) :: CBMC
!    real(dp) :: atmps = 1E-30_dp
!    real(dp) :: accpt = 0E0_dp

    real(dp), private, allocatable :: boxatmps(:)
    real(dp), private , allocatable:: boxaccpt(:)
    logical, private :: verbose = .true.
    logical, private :: proportional = .true.

    logical, private, allocatable :: validtype(:)
    integer, private, allocatable :: patharrays(:, :)

    type(Displacement), allocatable :: disp(:)
    real(dp), private , allocatable:: newpos(:, :)

    !Rejection Counters
    integer, private :: ovlaprej = 0 
    integer, private :: constrainrej = 0 
    integer, private :: detailedrej = 0 

!    integer, allocatable :: tempNnei(:)
!    integer, allocatable :: tempList(:, :)

    contains
      procedure, pass :: Constructor => CBMC_Constructor
!      procedure, pass :: GeneratePosition => CBMC_GeneratePosition
      procedure, pass :: FullMove => CBMC_FullMove
!      procedure, pass :: Maintenance => CBMC_Maintenance
      procedure, pass :: Prologue => CBMC_Prologue
      procedure, pass :: Update => CBMC_Update
      procedure, pass :: Epilogue => CBMC_Epilogue
      procedure, pass :: ProcessIO => CBMC_ProcessIO
  end type
!========================================================
 contains
!========================================================
  subroutine CBMC_Constructor(self)
    use BoxData, only: BoxArray
    use Common_MolInfo, only: MolData, nMolTypes, mostAtoms
    use MolCon_LinearCBMC, only: LinearCBMC
    use Template_MolConstructor, only: MolConstructor
    use ParallelVar, only: nout
    implicit none
    class(CBMC), intent(inout) :: self
    integer :: iType, maxAtoms, nAtoms, nBoxes

    nBoxes = size(boxArray)
    if(.not. allocated(self%boxProb)) then
      allocate( self%boxProb(1:nBoxes) )
      self%boxProb = 1E0_dp/real(nBoxes,dp)
    endif

    allocate( self%validtype(1:nMolTypes) )

    allocate( self%boxatmps(1:nBoxes) )
    allocate( self%boxaccpt(1:nBoxes) )
    self%boxatmps = 1e-50_dp
    self%boxaccpt = 0E0_dp

    maxAtoms = 0
    do iType = 1, nMolTypes
      if(MolData(iType)%nAtoms > maxAtoms) then
        maxAtoms = MolData(iType)%nAtoms 
      endif
    enddo


    allocate( self%disp(1:maxAtoms) )
    allocate( self%patharrays(1:maxAtoms, 1:nMolTypes) )
    self%patharrays = 0
    self%validtype = .true.

    do iType = 1, nMolTypes
      select type(molcon => MolData(iType)%molConstruct)
        class is(LinearCBMC)
          nAtoms = MolData(iType)%nAtoms 
          call molcon % GetPath( self%patharrays(1:nAtoms, iType) )

        class default
          self%validtype(iType) = .false.
      end select
    enddo

    write(nout,*) "CBMC Valid Molecule Types: ", self%validtype(1:nMolTypes)
    if( all(.not. self%validtype) ) then
      stop "CBMC Moves have been specified on molecules which do not use a CBMC regrowth style"
    endif


    call self%CreateTempArray(mostAtoms)
  end subroutine
!===============================================
  subroutine CBMC_FullMove(self, trialBox, accept) 
    use Box_Utility, only: FindAtom, FindMolecule
    use CommonSampling, only: sampling
    use Common_NeighData, only: neighSkin
    use Common_MolInfo, only: MolData, nMolTypes
    use ForcefieldData, only: EnergyCalculator
    use RandomGen, only: grnd
    implicit none
    class(CBMC), intent(inout) :: self
    class(SimpleBox), intent(inout) :: trialBox
    logical, intent(out) :: accept
    logical :: reverse
    integer :: boxID, iAtom, nAtoms, atomIndx
    integer :: nMove, rawIndx, iConstrain
    integer :: CalcIndex, molStart, molEnd, molType
    integer :: nRegrow, cutpoint
    integer :: lowIndx, highIndx, iDisp
    real(dp) :: dx, dy, dz
    real(dp) :: E_Diff, E_Inter, E_Intra, biasE
    real(dp) :: Prob, ProbFor, ProbRev
    
!    integer :: slice(1:2)
    real(dp), pointer :: atoms(:,:) => null()

    boxID = trialBox % boxID
    self % atmps = self % atmps + 1E0_dp
    self % boxatmps(boxID) = self % boxatmps(boxID) + 1E0_dp
    accept = .true.
    call self%LoadBoxInfo(trialBox, self%disp)

    !Propose move
    nMove = self%UniformMoleculeSelect(trialBox, restrict=self%validtype(1:nMolTypes))
    call trialBox % GetMolData(nMove, molStart=molStart, molEnd=molEnd, &
                               molType=molType)
    nAtoms = MolData(molType)%nAtoms

    cutpoint = floor( nAtoms * grnd() + 1E0_dp )
!    cutpoint = 3
!    reverse = .true.
    
    if(cutpoint == 1) then
      reverse = .true.
    elseif(cutpoint == nAtoms) then
      reverse = .false.
    else
      if( grnd() < 0.5E0_dp) then
        reverse = .false.
      else
        reverse = .true.
      endif
    endif

 
    if(reverse) then
      lowIndx = 1
      highIndx = cutpoint
    else
      lowIndx = cutpoint
      highIndx = nAtoms
    endif

    iDisp = 0
    nRegrow = 0
    do iAtom = lowIndx, highIndx
      iDisp = iDisp + 1
      atomIndx = molStart + self%patharrays(iAtom, molType) - 1
      nRegrow = nRegrow + 1

      self%disp(iDisp)%molType = molType
      self%disp(iDisp)%molIndx = nMove
      self%disp(iDisp)%atmIndx = atomIndx

      self%disp(iDisp)%x_new = 0E0_dp
      self%disp(iDisp)%y_new = 0E0_dp
      self%disp(iDisp)%z_new = 0E0_dp

      self%disp(iDisp)%newlist = .false.
      self%disp(iDisp)%listIndex = iDisp
    enddo

    call MolData(molType) % molConstruct % GenerateConfig(trialBox, self%disp(1:nRegrow), ProbFor, accept)
    if(.not. accept) then
      return
    endif
    Prob = 1E0_dp/ProbFor


    !If the particle moved a large distance get a temporary neighborlist

!    if(any([dx,dy,dz] > neighSkin)) then
!      call trialBox % NeighList(1) % GetTempListArray(tempList, tempNNei)
!      do iAtom = 1, nAtoms
!        call trialBox % NeighList(1) % GetNewList(iAtom, self%tempList, self%tempNNei, &
!                                                  self%disp(iAtom))
!        self%disp(iAtom)%listIndex = iAtom
!      enddo 
!    endif

    !Check Constraint
    accept = trialBox % CheckConstraint( self%disp(1:nRegrow) )
    if(.not. accept) then
      self%constrainrej = self%constrainrej + 1
      return
    endif

    !Energy Calculation
    call trialBox%ComputeEnergyDelta(self%disp(1:nRegrow),&
                                     self%templist,&
                                     self%tempNNei, &
                                     E_Inter, &
                                     E_Intra, &
                                     E_Diff, &
                                     accept, &
                                     computeintra=.true.)
    if(.not. accept) then
      self%ovlaprej = self%ovlaprej + 1
      return
    endif

    !Check Post Energy Constraint
    accept = trialBox % CheckPostEnergy( self%disp(1:nRegrow), E_Diff )
    if(.not. accept) then
      self%constrainrej = self%constrainrej + 1
      return
    endif



    call MolData(molType) % molConstruct % ReverseConfig(self%disp(1:nRegrow), trialBox, ProbRev, accept)

    Prob = Prob * ProbRev


    !Accept/Reject
!    if(abs(log(Prob) -trialBox%beta * E_Intra) > 1E-10) then
!      write(0,*) "?", log(Prob) -trialBox%beta * E_Intra
!      write(0,*) "reverse", -log(ProbRev)/trialbox%beta, trialBox%ETotal
!      write(0,*) "reverse", -log(ProbRev)/trialbox%beta, trialBox%ETotal
!      write(0,*) nRegrow, E_Diff, exp(-trialBox%beta * E_Diff), 1E0_dp/Prob 
!      write(0,*) nRegrow, E_Diff, exp(-trialBox%beta * E_Diff), 1E0_dp/Prob 
!      error stop
!    endif
    accept = sampling % MakeDecision(trialBox, E_Diff, self%disp(1:nRegrow), inProb=Prob)

    if(accept) then
!      write(*,*) "Accept", nRegrow

      self % accpt = self % accpt + 1E0_dp
      self % boxaccpt(boxID) = self % boxaccpt(boxID) + 1E0_dp
      call trialBox % UpdateEnergy(E_Diff, E_Inter, E_Intra)
      call trialBox % UpdatePosition(self%disp(1:nRegrow), self%tempList, self%tempNNei)
    else
!      slice(1) = molStart
!      slice(2) = molEnd
!      call trialbox%GetCoordinates(atoms, slice=slice)
!      write(*,*) "Reject", nRegrow, ProbFor, ProbRev
!      write(*,*) reverse, cutpoint
!      write(*,*) self%patharrays(lowIndx:highIndx, molType)
!        write(2,*) nAtoms
!         write(2,*) -log(ProbFor)/trialBox%beta, -log(ProbRev)/trialBox%beta, self%patharrays(lowIndx:highIndx, molType)
!         do iDisp = 1, nAtoms
!            write(2,*) "C", atoms(1:3, iDisp)
!         enddo
!         do iDisp = 1, nRegrow
!            write(2,*) "P", self%disp(iDisp)%x_new, self%disp(iDisp)%y_new, self%disp(iDisp)%z_new
!         enddo
!         write(2,*)
!      write(*,*) "forward", trialBox%ETotal-log(ProbFor)/trialbox%beta, trialBox%ETotal+E_diff
!      write(*,*) "reverse", -log(ProbRev)/trialbox%beta, trialBox%ETotal
!      write(*,*) nRegrow, E_Diff, exp(-trialBox%beta * E_Diff), Prob 
!      write(*,*) exp(-trialBox%beta * E_Diff)/Prob 
      self%detailedrej = self%detailedrej + 1
    endif

  end subroutine
!=========================================================================
  subroutine CBMC_Prologue(self)
    use ParallelVar, only: nout
    implicit none
    class(CBMC), intent(inout) :: self

    if(.not. allocated(self%disp)) then
      call self % Constructor
    endif
      

  end subroutine
!=========================================================================
  subroutine CBMC_Epilogue(self)
    use ParallelVar, only: nout
    implicit none
    class(CBMC), intent(inout) :: self
    real(dp) :: accptRate
      
    write(nout,*) 
    write(nout,"(1x,A,I15)") "CBMC Moves Accepted: ", nint(self%accpt)
    write(nout,"(1x,A,I15)") "CBMC Moves Attempted: ", nint(self%atmps)
    accptRate = self%GetAcceptRate()
    write(nout,"(1x,A,F15.8)") "CBMC Acceptance Rate: ", accptRate


    if(self%verbose) then
      write(nout, "(1x,A,I15)") "CBMC, Rejections due to overlap:", self%ovlaprej
      write(nout, "(1x,A,I15)") "CBMC, Rejections due to constraint:", self%constrainrej
      write(nout, "(1x,A,I15)") "CBMC, Rejections due to detailed balance:", self%detailedrej
    endif

  end subroutine
!=========================================================================
  subroutine CBMC_Update(self)
    use BoxData, only: BoxArray
    use ParallelVar, only: nout
    implicit none
    class(CBMC), intent(inout) :: self
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
  subroutine CBMC_ProcessIO(self, line, lineStat)
    use Input_Format, only: GetXCommand, maxLineLen
    implicit none
    class(CBMC), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line
    integer, intent(out) :: lineStat
    character(len=30) :: command
    logical :: logicVal
    integer :: intVal
    real(dp) :: realVal

    call GetXCommand(line, command, 4, lineStat)
    select case( trim(adjustl(command)) )
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
