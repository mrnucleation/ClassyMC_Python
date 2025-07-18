!===========================================================================
! Main Monte Carlo simulation subroutine.  This makes all the calls
! to the various sub components of this simulation.
!===========================================================================
#define __StdErr__ 0
!===========================================================================
module SimMonteCarlo
  use ParallelVar, only: myid, ierror, nout
  use VarPrecision
!===========================================================================
contains
!===========================================================================
  subroutine RunMonteCarlo
#ifdef MPIPARALLEL
    use MPI
#endif

    use AnalysisData, only: AnalysisArray
    use BoxData, only: BoxArray
    use CommonSampling, only: Sampling
    use Common_MolInfo, only: nAtomTypes, nMolTypes, MolData
!    use Debug, only: Debug_DumpNeiList
!    use ForcefieldData, only: EnergyCalculator
    use MCMoveData, only: Moves, MoveProb
    use MoveClassDef, only: MCMove
    use MultiBoxMoveDef, only: MCMultiBoxMove
    use Output_DumpCoords, only: Output_DumpData
    use RandomGen, only: sgrnd, grnd, ListRNG
    use SimControl, only: nMoves, nCycles, screenfreq, configfreq, energyCheck
    use SimControl, only: TimeStart
    use Units, only: outEngUnit

    implicit none
 
    logical :: accept
    integer :: i, j, nBoxes, iBox
    integer :: iAtom, moveNum, boxNum
    integer(kind=8) :: iCycle, iMove
    character(len=50) :: fileName
!    class(MCMove), pointer :: curMove
    real(dp), allocatable :: boxProb(:)

    iCycle = 0
    iMove = 0



    write(nout,*) "Starting Pre-Simulation Checks...."
    call Prologue
    call SafetyCheck
    nBoxes = size(BoxArray)
    allocate( boxProb(1:nBoxes) )
    boxNum = 1

    accept = .true.

    call Analyze(iCycle, iMove, accept, .true.)
    call Analyze(iCycle, iMove, accept, .false.)

    if(nCycles < 1) then
      write(nout, *) "============================================"
      write(nout,*) "Number of Cycles is less than 1!"
      write(nout,*) "Run Command was issued, but nothing is done"
      write(nout, *) "============================================"
      return
    endif
    !-------Main Monte Carlo Simulation Loop-------
    write(nout, *) "============================================"
    write(nout, *) "       Simulation Start!"
    write(nout, *) "============================================"

    flush(nout)


    call CPU_TIME(TimeStart)
    call ScreenOut(iCycle, iMove)

    do iCycle = 1, nCycles

      !-----Start Move Loop
      do iMove = 1, nMoves
        accept = .true.
        moveNum = ListRNG(MoveProb) !Randomly select a move to perform
        select type( curmove => Moves(moveNum) % Move )
          class is (MCMultiBoxMove) ! Mutli Box Move
!            call Moves(moveNum) % Move % MultiBox (accept)
            call curmove % MultiBox (accept)
            boxNum = -1 !Set to -1 since multiple boxes are changes
          class is (MCMove) ! Single Box Move
            if(nBoxes > 1) then
              call Moves(moveNum) % Move % GetBoxProb(boxProb)
              boxNum = ListRNG(boxProb)
            else
              boxNum = 1
            endif
            call Moves(moveNum) % Move % FullMove(BoxArray(boxNum)%box, accept)
        end select
        call Sampling%UpdateStatistics(accept)
        if(accept) then
          call Update(accept)
        endif

        call Analyze(iCycle, iMove, accept, .true.) !Per Move Analysis
        if(accept) then
          do iBox = 1, size(BoxArray)
            if( (boxNum < 0) .or. (boxNum == iBox) ) then
              call BoxArray(iBox)%box%CheckLists
            endif
          enddo
        endif
      enddo 
      !------End Move Loop
      if(mod(iCycle, int(screenfreq,8)) == 0) then
        call ScreenOut(iCycle, iMove)
        flush(nout)
      endif

      if( mod(iCycle, int(configfreq, 8) ) == 0) then
        call Output_DumpData
      endif

      if(energyCheck > 0) then
        if( mod(iCycle, int(energyCheck, 8) ) == 0) then
          do boxNum = 1, size(BoxArray)
            call BoxArray(boxNum)%box%EnergySafetyCheck
          enddo
        endif
      endif



      call Analyze(iCycle, iMove, accept, .false.) !Per Cycle Analysis
      call Maintenance(iCycle, iMove)
      call Trajectory(iCycle, iMove)
    enddo
    !-------End of Main Monte Carlo Simulation Loop-------
 
    call ScreenOut(iCycle, iMove)
    write(nout,*) "======================================="
    write(nout,*) "     Simulation End"
    write(nout,*) "======================================="

    write(nout,*) "     Beginning Epilogue....."


    call Epilogue

#ifdef MPIPARALLEL
    call MPI_BARRIER(MPI_COMM_WORLD, ierror)       
#endif

    if( allocated(AnalysisArray) ) then
      do i = 1, size(AnalysisArray)
        call AnalysisArray(i) % func % Finalize
      enddo
    endif

    call Output_DumpData
      
  end subroutine
!===========================================================================
  subroutine Analyze(iCycle, iMove, accept, moveloop)
    use AnalysisData, only: AnalysisArray
    implicit none
    integer(kind=8), intent(in) :: iCycle, iMove
    logical, intent(in) :: accept
    logical, intent(in) :: moveloop
    integer :: iAn

    if( allocated(AnalysisArray) ) then
      do iAn = 1, size(AnalysisArray)
        if( AnalysisArray(iAn)%func%perMove .eqv. moveloop) then
          if((mod(iCycle, int(AnalysisArray(iAn)%func%UpdateFreq,8)) == 0) .or.  AnalysisArray(iAn)%func%perMove) then
            call AnalysisArray(iAn) % func % Compute(accept)
          endif
        endif
      enddo
    endif
 
  end subroutine
!===========================================================================
  subroutine ScreenOut(iCycle, iMove)
    use AnalysisData, only: AnalysisArray
    use BoxData, only: BoxArray
    use CommonSampling, only: Sampling
    use ForcefieldData, only: EnergyCalculator
    use Input_Format, only: ReplaceText
    use MCMoveData, only: Moves, MoveProb
    use SimControl, only: TimeStart
    use TrajData, only: TrajArray
    implicit none
    integer(kind=8), intent(in) :: iCycle, iMove
    integer :: i
    character(len=80) :: tempStr, tempStr2
    real(dp) :: CurrentTime, timepercycle

    write(tempStr, "(A)") "    ----------------- Cycle Number: %s ----------------------"
    write(tempStr2, "(I60)") iCycle
    tempStr = ReplaceText(tempStr, "%s", trim(adjustl(tempStr2)))
    write(nout, "(A)") tempStr

    call CPU_TIME(CurrentTime)

    write(tempStr, "(A)") "   Simulation Time: %s sec "
    write(tempStr2, "(F60.3)") CurrentTime - TimeStart
    tempStr = ReplaceText(tempStr, "%s", trim(adjustl(tempStr2)))
    write(nout, "(A)") tempStr

    if(iCycle > 0) then
        write(tempStr, "(A)") "   Time per Cycle: %s sec "
        timepercycle =  (CurrentTime - TimeStart)/real(iCycle, dp)
        if(timepercycle > 1e-3_dp) then
          write(tempStr2, "(F60.3)") timepercycle
        else
          write(tempStr2, "(E60.3)") timepercycle
        endif
        tempStr = ReplaceText(tempStr, "%s", trim(adjustl(tempStr2)))
        write(nout, "(A)") tempStr
    endif

    call Sampling % ScreenOut

    if( allocated(AnalysisArray) ) then
      do i = 1, size(AnalysisArray)
        call AnalysisArray(i) % func % ScreenOut
      enddo
    endif

    if( allocated(TrajArray) ) then
      do i = 1, size(TrajArray)
        call TrajArray(i) % traj % ScreenOut
      enddo
    endif

    do i = 1, size(EnergyCalculator)
      call EnergyCalculator(i)%method%ScreenOut
    enddo


    do i = 1, size(BoxArray)
      call BoxArray(i) % box % ScreenOut
    enddo

    do i = 1, size(Moves)
      call Moves(i) % move % ScreenOut
    enddo

    write(nout, *) 

  end subroutine
!===========================================================================
  subroutine Maintenance(iCycle, iMove)
    use AnalysisData, only: AnalysisArray
    use BoxData, only: BoxArray
    use MCMoveData, only: Moves, MoveProb
    use TrajData, only: TrajArray
    use CommonSampling, only: Sampling
    implicit none
    integer(kind=8), intent(in) :: iCycle, iMove
    integer :: i

    if(mod(iCycle, int(Sampling%maintFreq,8)) == 0 ) then
      call Sampling % Maintenance
    endif

    if( allocated(AnalysisArray) ) then
      do i = 1, size(AnalysisArray)
        if(mod(iCycle, int(AnalysisArray(i)%func%maintFreq,8)) == 0) then
          call AnalysisArray(i) % func % Maintenance
        endif
      enddo
    endif

    if( allocated(TrajArray) ) then
      do i = 1, size(TrajArray)
        if(mod(iCycle, int(TrajArray(i)%traj%maintFreq,8)) == 0) then
          call TrajArray(i) % traj % Maintenance
        endif
      enddo
    endif

    do i = 1, size(BoxArray)
!      if(mod(iCycle, BoxArray(i)%box%maintFreq) == 0) then
        call BoxArray(i) % box % Maintenance
!      endif
    enddo


    do i = 1, size(Moves)
!      if(mod(iCycle, int(Moves(i)%move%maintFreq,8)) == 0) then
      if(mod(iCycle, Moves(i)%move%maintFreq) == 0) then
        call Moves(i) % move % Maintenance
      endif
    enddo

  end subroutine
!===========================================================================
  subroutine Trajectory(iCycle, iMove)
    use TrajData, only: TrajArray
    implicit none
    integer(kind=8), intent(in) :: iCycle, iMove
    integer :: iTraj

    if( allocated(TrajArray) ) then
      do iTraj = 1, size(TrajArray)
        if(mod(iCycle, int(TrajArray(iTraj)%traj%outfreq,8)) == 0) then
          call TrajArray(iTraj) % traj % WriteFrame(iCycle)
        endif
      enddo
    endif

  end subroutine
!===========================================================================
  subroutine Update(accept)
    use AnalysisData, only: AnalysisArray
    use BoxData, only: BoxArray
    use ForcefieldData, only: EnergyCalculator
    use MCMoveData, only: Moves, MoveProb
    use TrajData, only: TrajArray
    use CommonSampling, only: Sampling
    implicit none
    logical, intent(in) :: accept
    integer :: i

    if( .not. accept) then
      return
    endif

    call Sampling % Update

    if( allocated(AnalysisArray) ) then
      do i = 1, size(AnalysisArray)
        call AnalysisArray(i) % func % Update
      enddo
    endif

    if( allocated(TrajArray) ) then
      do i = 1, size(TrajArray)
        call TrajArray(i) % traj % Update
      enddo
    endif

    do i = 1, size(EnergyCalculator)
      call EnergyCalculator(i) % method % Update
    enddo

    do i = 1, size(BoxArray)
      call BoxArray(i) % box % Update
    enddo

    do i = 1, size(Moves)
      call Moves(i) % move % Update
    enddo

  end subroutine
!===========================================================================
  subroutine SafetyCheck
    use AnalysisData, only: AnalysisArray
    use BoxData, only: BoxArray
    use Common_MolInfo, only: MolData
    use MCMoveData, only: Moves, MoveProb
    use TrajData, only: TrajArray
    use CommonSampling, only: Sampling
    use ParallelVar, only: nout
    implicit none
    integer :: i

    if(.not. allocated(MolData)) then
      write(__StdErr__,*) "*******************************************************************************"
      write(__StdErr__,*) "  CRITICAL ERROR! Molecular Topology Information has not been defined!"
      write(__StdErr__,*) "*******************************************************************************"
      error stop
    endif
    do i = 1, size(MolData)
      if(.not. allocated(MolData(i)%molConstruct)) then
        write(__StdErr__,*) "*******************************************************************************"
        write(__StdErr__,*) "  WARNING! Molecule reconstructor is not defined in the forcefield file!"
        write(__StdErr__,*) "  Swap moves and any move which generates a new configuration from scratch will not work!"
        write(__StdErr__,*) "*******************************************************************************"
        write(nout,*) "*******************************************************************************"
        write(nout,*) "  WARNING! Molecule reconstructor is not defined in the forcefield file!"
        write(nout,*) "  Swap moves and any move which generates a new configuration from scratch will not work!"
        write(nout,*) "*******************************************************************************"

      else
        call MolData(i) % molConstruct % SafetyCheck
      endif
    enddo

    if(.not. allocated(Sampling) ) then
      write(__StdErr__, *) "*******************************************************************************"
      write(__StdErr__, *) "  CRITICAL ERROR! Sampling Rule has not been defined!"
      write(__StdErr__, *) "*******************************************************************************"
      error stop
    endif

    call Sampling % SafetyCheck

    if( allocated(AnalysisArray) ) then
      do i = 1, size(AnalysisArray)
        call AnalysisArray(i) % func % SafetyCheck
      enddo
    endif

    if( allocated(TrajArray) ) then
      do i = 1, size(TrajArray)
        call TrajArray(i) % traj % SafetyCheck
      enddo
    endif

    if( allocated(BoxArray) ) then
      do i = 1, size(BoxArray)
        call BoxArray(i) % box % SafetyCheck
      enddo
    else
      write(__StdErr__,*) "*******************************************************************************"
      write(__StdErr__,*) "  CRITICAL ERROR! No Simulation Boxes have not been defined!"
      write(__StdErr__,*) "*******************************************************************************"
      error stop
    endif

    if( allocated(BoxArray) ) then
      do i = 1, size(Moves)
        call Moves(i) % move % SafetyCheck
      enddo
    else
      write(__StdErr__,*) "*******************************************************************************"
      write(__StdErr__,*) "  WARNING! No Monte Carlo Moves have been defined!"
      write(__StdErr__,*) "  Nothing will move! Are you ok with this?"
      write(__StdErr__,*) "*******************************************************************************"
    endif

  end subroutine
!===========================================================================
  subroutine Prologue
    use AnalysisData, only: AnalysisArray
    use BoxData, only: BoxArray
    use CommonSampling, only: Sampling
    use ForcefieldData, only: EnergyCalculator
    use Common_MolInfo, only: MolData
    use MCMoveData, only: Moves, MoveProb
    use TrajData, only: TrajArray
    implicit none
    integer :: i, iMisc

    do i = 1, size(MolData)
      if(allocated(MolData(i) % molConstruct) ) then
        call MolData(i) % molConstruct % Prologue
      endif
      if(allocated(MolData(i) % miscdata) ) then
        do iMisc = 1, MolData(i) % nMisc
          call MolData(i) % MiscData(iMisc) % miscFF % Prologue
        enddo
      endif
    enddo


    if( allocated(AnalysisArray) ) then
      do i = 1, size(AnalysisArray)
        call AnalysisArray(i) % func % Prologue
      enddo
    endif

    call Sampling % Prologue
!    write(nout, *) "Traj"
    if( allocated(TrajArray) ) then
      do i = 1, size(TrajArray)
        call TrajArray(i) % traj % Prologue
      enddo
    endif

!    write(nout,*) "ECalc"
    do i = 1, size(EnergyCalculator)
      call EnergyCalculator(i)%method%Prologue
    enddo


!    write(nout,*) "Box"
    do i = 1, size(BoxArray)
      call BoxArray(i) % box % Prologue
    enddo

!    write(nout,*) "Moves"
    flush(nout)
    do i = 1, size(Moves)
      call Moves(i) % move % Prologue
    enddo

  end subroutine
!===========================================================================
  subroutine Epilogue
    use AnalysisData, only: AnalysisArray
    use BoxData, only: BoxArray
    use Common_MolInfo, only: MolData
    use MCMoveData, only: Moves, MoveProb
    use TrajData, only: TrajArray
    use CommonSampling, only: Sampling
    use ParallelVar, only: nout
    implicit none
    integer :: i, iMisc

    do i = 1, size(BoxArray)
      call BoxArray(i) % box % Epilogue
    enddo
    write(nout, *) "-----------------------"

    call Sampling % Epilogue
    
    do i = 1, size(MolData)
      if(allocated(MolData(i) % molConstruct)) then
        call MolData(i) % molConstruct % Epilogue
      endif
      if(allocated(MolData(i) % miscdata) ) then
        do iMisc = 1, MolData(i) % nMisc
          call MolData(i) % MiscData(iMisc) % miscFF % Epilogue
        enddo
      endif
    enddo

    if( allocated(AnalysisArray) ) then
      do i = 1, size(AnalysisArray)
        call AnalysisArray(i) % func % Epilogue
      enddo
    endif

    if( allocated(TrajArray) ) then
      do i = 1, size(TrajArray)
        call TrajArray(i) % traj % Epilogue
      enddo
    endif


    do i = 1, size(Moves)
      call Moves(i) % move % Epilogue
    enddo

  end subroutine

!===========================================================================
end module
!===========================================================================
