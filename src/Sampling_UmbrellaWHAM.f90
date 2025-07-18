!====================================================================
module UmbrellaWHAMRule
  use VarPrecision
  use Template_AcceptRule, only: acceptrule
  use CoordinateTypes, only: Displacement, Perturbation, Addition, Deletion, VolChange
 
  type, public, extends(AcceptRule) :: UmbrellaWHAM

    logical :: lastaccept = .false.
    integer :: nBiasVar = 0
    integer, allocatable :: AnalysisIndex(:)

    integer :: umbrellaLimit = 0
    integer, private :: potfile = 0
    integer, private :: histfile = 0
    integer, private :: whamfile = 0
    integer, allocatable :: varType(:)
    integer, allocatable :: indexCoeff(:)
    integer, allocatable :: binMin(:), binMax(:)
    integer, allocatable :: binIndx(:), nBins(:)
    integer, allocatable :: UArray(:)
    integer :: oldIndx, newIndx
    real(dp), allocatable :: valMin(:), valMax(:)
    real(dp), allocatable :: UBias(:)
    real(dp), allocatable :: UHist(:)
    real(dp), allocatable :: UBinSize(:)
    real(dp), allocatable :: varValues(:)
    character(len=50) :: fileName = "umbrella.dat"

    real(dp), allocatable :: refVals(:)
    integer :: refBin = 1
    integer :: maxSelfConsist = 3000
    real(dp) :: tolLimit = 1E-5_dp
!    integer :: equilInterval
!    integer :: whamEstInterval

    integer :: nWhamItter, nCurWhamItter

    real(dp), allocatable :: WHAM_Numerator(:)
    real(dp), allocatable :: WHAM_Denominator(:,:)
    real(dp), allocatable :: HistStorage(:)
    real(dp), allocatable :: FreeEnergyEst(:)
    real(dp), allocatable :: BiasStorage(:,:)
    real(dp), allocatable :: NewBias(:)
    real(dp), allocatable :: ProbArray(:)      
    real(dp), allocatable :: TempHist(:)
 

    contains
       procedure, pass :: Constructor => UmbrellaWHAM_Constructor
       procedure, pass :: MakeDecision => UmbrellaWHAM_MakeDecision
       procedure, pass :: MakeDecision2Box => UmbrellaWHAM_MakeDecision2Box
       procedure, pass :: UpdateStatistics => UmbrellaWHAM_UpdateStatistics

       procedure, pass :: GetBiasIndex => UmbrellaWHAM_GetBiasIndex
       procedure, pass :: GetNewBiasIndex => UmbrellaWHAM_GetNewBiasIndex
       procedure, pass :: ReadInitialBias => UmbrellaWHAM_ReadInitialBias
       procedure, pass :: GetUIndexArray => UmbrellaWHAM_GetUIndexArray
       procedure, pass :: FindVarValues => UmbrellaWHAM_FindVarValues
       procedure, pass :: OutputUmbrellaHist => UmbrellaWHAM_OutputUmbrellaHist
       procedure, pass :: AdjustHist => UmbrellaWHAM_AdjustHist
       procedure, pass :: WHAMSimOutput => UmbrellaWHAM_WHAMSimOutput
       procedure, pass :: Maintenance => UmbrellaWHAM_Maintenance

       procedure, pass :: ProcessIO => UmbrellaWHAM_ProcessIO
       procedure, pass :: Epilogue => UmbrellaWHAM_Epilogue
       procedure, pass :: Prologue => UmbrellaWHAM_Prologue
!       procedure, pass :: Update => UmbrellaWHAM_Update

  end type
!====================================================================
  contains
!====================================================================
  subroutine UmbrellaWHAM_Constructor(self)
    implicit none
    class(UmbrellaWHAM), intent(inout) :: self
    integer :: AllocationStat

    if(self%nBiasVar < 1) then
      write(0,*) "No Bias Variables have been specified for the UmbrellaWHAM sampling method"
      write(0,*) "Please defined at least one."
      stop
    endif

    allocate( self%refVals(1:self%nBiasVar), stat=AllocationStat ) 
    allocate( self%AnalysisIndex(1:self%nBiasVar), stat=AllocationStat ) 
    allocate( self%varType(1:self%nBiasVar), stat=AllocationStat ) 
    allocate( self%indexCoeff(1:self%nBiasVar), stat=AllocationStat ) 
    allocate( self%binMax(1:self%nBiasVar), stat=AllocationStat) 
    allocate( self%binMin(1:self%nBiasVar), stat=AllocationStat) 
    allocate( self%nBins(1:self%nBiasVar), stat=AllocationStat) 
    allocate( self%binIndx(1:self%nBiasVar), stat=AllocationStat) 

    allocate( self%UBinSize(1:self%nBiasVar), STAT =AllocationStat )
    allocate( self%UArray(1:self%nBiasVar), STAT =AllocationStat )
    allocate( self%varValues(1:self%nBiasVar), STAT =AllocationStat )
    allocate( self%valMin(1:self%nBiasVar), STAT =AllocationStat )
    allocate( self%valMax(1:self%nBiasVar), STAT =AllocationStat )
    self%varType = 1
    self%UBinSize = 0
    self%UArray = 0
    self%refVals = 0E0_dp
    self%valMin = 0E0_dp
    self%valMax = 0E0_dp
    self%varValues = 0E0_dp
  end subroutine
!====================================================================
  subroutine UmbrellaWHAM_Prologue(self)
    use AnalysisData, only: AnalysisArray
    use SimControl, only: nCycles
    use ParallelVar, only: nout, myid
    implicit none
    class(UmbrellaWHAM), intent(inout) :: self
    integer :: i,j, indx, stat, AllocateStatus

    if(self%nBiasVar < 1) then
      write(0,*) "No Bias Variables have been specified!"
      write(0,*) "Number of variables is set to 0!"
      stop
    endif

    do i = 1, self%nBiasVar
      indx = self%AnalysisIndex(i)
      if((indx > size(AnalysisArray)) .or. (indx < 1) ) then
        write(nout, *) "ERROR! The Umbrella Sampling routine has been directed to an "
        write(nout, *) "invalid Analysis fucntion"
        write(nout, *) "Chosen Function:", indx
        write(nout, *) "Number of Analysis Functions:", size(AnalysisArray)
        error stop "Error dectected in Umbrella Sampling"
      endif
    enddo

    do i = 1, self%nBiasVar
      if(self%valMin(i) > self%valMax(i) ) then
        write(nout,*) "ERROR! The given bounds for one of the umbrella variables does not make sense!"
        write(nout,*) "Smallest bin is larger than the largest bin" 
        write(nout,*) "Minimum Value:", self%valMin(i)
        write(nout,*) "Maximum Value:", self%valMax(i)
        error stop
      endif
    enddo

    do i = 1, self%nBiasVar
      indx = self%AnalysisIndex(i)
      AnalysisArray(indx)%func%usedInMove = .true.
      AnalysisArray(indx)%func%permove = .true.
    enddo

     ! Since the number of biasing variables is only known at run time, the bias matrix
     ! must be stored in a 1D array.  The index coefficient variables are used to emulate a N dimension matrix
     ! using the linear mapping equation of this form:
     ! U =  a1*x1 + a2*x2 + a3*x3 + .....
     ! Which maps a N dimension matrix onto a 1D array. 

    self%indexCoeff(1) = 1
    do i = 2, self%nBiasVar 
      self%indexCoeff(i) = 1
      do j = 1, i-1
        self%indexCoeff(i) = self%indexCoeff(i) + self%indexCoeff(j) * self%nBins(j) 
      enddo
    enddo      
    self%umbrellaLimit = 1
    do i = 1, self%nBiasVar 
      self%umbrellaLimit = self%umbrellaLimit + self%indexCoeff(i) * self%nBins(i)
    enddo
    self%umbrellaLimit = max(self%umbrellaLimit, 1)

    write(nout,*) "Sampling Style: Histogram based Umbrella Sampling \w Auto WHAM method"
    write(nout,*) "Number of Umbrella Bins:", self%umbrellaLimit
       
    allocate(self%UBias(1:self%umbrellaLimit+1), STAT = AllocateStatus)
    allocate(self%UHist(1:self%umbrellaLimit+1), STAT = AllocateStatus)
    allocate(self%TempHist(1:self%umbrellaLimit+1), STAT = AllocateStatus)      
    allocate(self%NewBias(1:self%umbrellaLimit+1), STAT = AllocateStatus)
    self%NewBias = 0E0_dp
    self%TempHist = 0E0_dp


!    write(nout,*) self%binMax
    self%UBias = 0E0_dp
    self%UHist = 0E0_dp

    call self%GetUIndexArray(self%refVals, i, stat)
    if(stat /= 0) then
      error stop "Index Error Encountered in Umbrella Sampling"
    endif
    self%refBin = i

    call self%ReadInitialBias
    self%nWhamItter = ceiling(dble(nCycles)/dble(self%maintFreq))
    ! Allocation of the WHAM variables
    if(myid .eq. 0) then
      allocate(self%WHAM_Numerator(1:self%umbrellaLimit+1), STAT = AllocateStatus)
      allocate(self%WHAM_Denominator(1:self%umbrellaLimit+1,1:self%nWhamItter+1), STAT = AllocateStatus)
      allocate(self%HistStorage(1:self%umbrellaLimit+1), STAT = AllocateStatus)
      allocate(self%BiasStorage(1:self%umbrellaLimit+1,1:self%nWhamItter+1), STAT = AllocateStatus)
      allocate(self%FreeEnergyEst(1:self%umbrellaLimit+1), STAT = AllocateStatus)
      allocate(self%ProbArray(1:self%umbrellaLimit+1), STAT = AllocateStatus)

      self%FreeEnergyEst = 0E0_dp
      self%WHAM_Numerator = 0E0_dp
      self%WHAM_Denominator = 0E0_dp
      self%HistStorage = 0E0_dp
      self%BiasStorage = 0E0_dp
      self%ProbArray = 0E0_dp
      self%nCurWhamItter = 1
!        tolLimit = 1d-2
      open(newunit = self%histfile, file="Umbrella_Histogram.dat")
      open(newunit = self%potfile, file="Umbrella_Potential.dat")
      open(newunit = self%whamfile, file="Umbrella_WHAMDG.dat")
    endif


    write(nout,*) self%refVals, i 
    write(nout,*) "Bin Size:", self%UBinSize
!    self%oldIndx = self % GetBiasIndex()
    self%oldIndx = 1
  end subroutine
!====================================================================
  function UmbrellaWHAM_MakeDecision(self, trialBox, E_Diff, disp, inProb, logProb, extraIn) result(accept)
    use AnalysisData, only: AnalysisArray
    use Template_SimBox, only: SimBox
    use RandomGen, only: grnd
    implicit none
    class(UmbrellaWHAM), intent(inout) :: self
    class(SimBox), intent(in) :: trialBox
    class(Perturbation), intent(in) :: disp(:)
    real(dp), intent(in), optional :: inProb, logProb
    real(dp), intent(in), optional:: extraIn
    real(dp), intent(in) :: E_Diff

    logical :: accept, indxchanged
    logical :: usebias
    integer :: iBias, indx
    real(dp) :: biasE, biasOld, biasNew
    real(dp) :: extraTerms, probTerm

    accept = .true.
    if(present(inProb)) then
      if(inProb <= 0E0_dp) then
        accept = .false.
        return
      endif
      probTerm = log(inProb)
    elseif(present(logProb)) then
      probTerm = logProb
    else
      write(0,*) "Coding Error! Probability has not been passed into Sampling "
      error stop
    endif

    indxchanged = .false.
    do iBias = 1, self%nBiasVar
      indx = self%AnalysisIndex(iBias)
      call AnalysisArray(indx) % func % CalcNewState(disp, accept=usebias)
      if(usebias) indxchanged = .true.
    enddo

    self%oldIndx = self % GetBiasIndex()
    if(indxchanged) then
      call self%GetNewBiasIndex(self%newIndx, accept)
    else
      self%newIndx = self%oldIndx
    endif

    if(.not. accept) then
!      write(*,*) "???"
      return
    endif

    biasOld = self%UBias(self%oldIndx)
    biasNew = self%UBias(self%newIndx)

    if(present(extraIn)) then
      extraTerms = extraIn
    else
      extraTerms = 0E0_dp
    endif


    accept = .false.
    biasE = -trialBox%beta * E_Diff + probTerm + (biasNew-biasOld) + extraTerms
!    if( present(extraIn) ) then
!      write(*,*) "Bias1B", E_diff, -trialBox%beta * E_Diff
!      write(*,*) "    ", biasOld, biasNew,  probTerm, extraTerms
!      write(*,*) "    ",  biasE
!      write(*,*)
!    endif
    if(biasE >= 0.0E0_dp) then
      accept = .true.
    elseif(biasE > log(grnd()) ) then
      accept = .true.
    endif

  end function
!====================================================================
  function UmbrellaWHAM_MakeDecision2Box(self, trialBox1,  trialBox2, E_Diff1, E_Diff2, &
                                       disp1, disp2, inProb, &
                                       logProb, extraIn ) result(accept)
    use Template_SimBox, only: SimBox
    use AnalysisData, only: AnalysisArray
    use RandomGen, only: grnd
    use ParallelVar, only: nout
    implicit none
    class(UmbrellaWHAM), intent(inout) :: self
    class(SimBox), intent(in) :: trialBox1, trialBox2
    class(Perturbation), intent(in) :: disp1(:), disp2(:)
    real(dp), intent(in) :: E_Diff1, E_Diff2
    real(dp), intent(in), optional :: inProb, logProb, extraIn
    logical :: accept
    logical :: usebias
    integer :: iDisp, iBias
    integer :: indx
    integer :: biasNew, biasOld
    real(dp) :: biasE, chemPot, extraTerms, probTerm, rannum



    accept = .false.
    if(present(inProb)) then
      if(inProb <= 0E0_dp) then
        return
      endif
    endif

    if(present(extraIn)) then
      extraTerms = extraIn
    else
      extraTerms = 0E0_dp
    endif

    if(present(inProb)) then
      probTerm = log(inProb)
    elseif(present(logProb)) then
      probTerm = logProb
    else
      write(0,*) "Coding Error! Probability has not been passed into Sampling "
      error stop
    endif

    do iBias = 1, self%nBiasVar
      indx = self%AnalysisIndex(iBias)
      call AnalysisArray(indx) % func % CalcNewState(disp1, usebias)
      call AnalysisArray(indx) % func % CalcNewState(disp2, usebias)
    enddo

    self%oldIndx = self % GetBiasIndex()
    call self%GetNewBiasIndex(self%newIndx, accept)
!    write(*,*) self%oldIndx, self%newIndx

    if(.not. accept) then
      return
    endif

    biasOld = self%UBias(self%oldIndx)
    biasNew = self%UBias(self%newIndx)


!    write(*,*) "Bias2b", E_diff1, E_Diff2
!    write(*,*) "Bias", -trialBox1%beta*E_diff1, -trialBox2%beta*E_Diff2
!    write(*,*) "    ", biasOld, biasNew
!    write(*,*) "    ", probTerm, extraTerms
    biasE = -trialBox1%beta*E_Diff1 - trialBox2%beta*E_Diff2 + probTerm + extraTerms + (biasNew-biasOld)
    rannum = grnd()
!    write(*,*) "    ",  biasE, log(rannum)
    accept = .false.
    if(biasE >= 0.0E0_dp) then
      accept = .true.
!      write(*,*) "???"
    elseif( biasE >= log(rannum) ) then
      accept = .true.
    endif
!    write(*,*) "Accept: ", accept
!    write(*,*)

  end function

!==========================================================================
  function UmbrellaWHAM_GetBiasIndex(self)  result(biasIndx)
    use AnalysisData, only: AnalysisArray, analyCommon
    implicit none
    class(UmbrellaWHAM), intent(inout) :: self
    integer :: biasIndx

    integer :: analyIndx
    integer :: iBias, intVal
    real(dp) :: biasVal
    
   ! Get the variables that are used for the biasing and figure out which histogram bin they
   ! fall into. 
    do iBias = 1, self%nBiasVar
      analyIndx = self%AnalysisIndex(iBias)
      biasVal = AnalysisArray(analyIndx)%func%GetResult()
!      write(2,*) biasVal
!      self%binIndx(iBias) = nint(biasVal/self%UBinSize(iBias))
      select type( biasVar => analyCommon(analyIndx)%val )
        type is(integer)
            intVal = nint(biasVal)
            self%binIndx(iBias) = intVal - nint(self%valMin(iBias))
!            write(2,*) intVal, self%binIndx(iBias), nint(self%valMin(iBias))

        type is(real(dp))
!            biasVal = biasVar
            self%binIndx(iBias) = floor((biasVal-self%valMin(iBias))/self%UBinSize(iBias))
       end select


!      self%binIndx(iBias) = floor((biasVal-self%valMin(iBias))/self%UBinSize(iBias))
    enddo


   ! Using the bin values from each biasing variable, determine which
    biasIndx = 1
    do iBias = 1, self%nBiasVar
!      biasIndx = biasIndx + self%indexCoeff(iBias) * ( self%binIndx(iBias) - self%binMin(iBias) )
      biasIndx = biasIndx + self%indexCoeff(iBias) * self%binIndx(iBias)
    enddo
!    write(2,*) biasIndx
!    write(2,*) 

  if( (biasIndx < 1) .or. (biasIndx > self%umbrellaLimit) ) then
    write(0,*) "ERROR! Umbrella Bias Index out of bounds!"
    write(0,*) "This may be due to the initial system configuration being outside the Umbrella sampling range"
    write(0,*) "or a calculation error in the analysis module being used."
    error stop "ERROR! Umbrella Bias Index out of bounds!"
  endif


  end function
!==========================================================================
  subroutine UmbrellaWHAM_GetNewBiasIndex(self, biasIndx, accept)
    use AnalysisData, only: analyCommon
    implicit none
    class(UmbrellaWHAM), intent(inout) :: self
    logical, intent(out) :: accept
    integer, intent(out) :: biasIndx

    integer :: analyIndx
    integer :: iBias
    real(dp) :: biasVal
    
   ! Get the variables that are used for the biasing and figure out which histogram bin they
   ! fall into. 
    accept = .true.
    biasVal = 0E0_dp
    do iBias = 1, self%nBiasVar
      analyIndx = self%AnalysisIndex(iBias)
      select type( biasVar => analyCommon(analyIndx)%val )
        type is(integer)
!             write(*,*) "new vals", biasVar, self%valMax(iBias),  self%valMin(iBias)
            biasVal = real(biasVar, dp)
        type is(real(dp))
            biasVal = real(biasVar, dp)
!            write(*,*) biasVar
        class default
             write(*,*) "WHAT TYPE IS THIS?! I DON'T KNOW, HELP!"
      end select

!      write(*,*) self%valMax(iBias), biasVal
      if(biasVal > self%valMax(iBias) ) then
        accept = .false.
        return
      endif
      if(biasVal < self%valMin(iBias) ) then
        accept = .false.
        return
      endif

      select type( biasVar => analyCommon(analyIndx)%val )
        type is(integer)
            self%binIndx(iBias) = biasVar-nint(self%valMin(iBias))
!            write(*,*) biasVar, self%binIndx(iBias), self%valMin(iBias),self%UBinSize(iBias)
        type is(real(dp))
            self%binIndx(iBias) = floor((biasVar-self%valMin(iBias))/self%UBinSize(iBias))
!            write(*,*) biasVar, self%binIndx(iBias), self%valMin(iBias),self%UBinSize(iBias)
      end select
!      self%binIndx(iBias) = floor((biasVal-self%valMin(iBias))/self%UBinSize(iBias))
    enddo


   ! Using the bin values from each biasing variable, determine which
    biasIndx = 1
    do iBias = 1, self%nBiasVar
!       write(*,*) biasIndx, self%indexCoeff(iBias), self%binIndx(iBias), self%binMin(iBias) 
!      biasIndx = biasIndx + self%indexCoeff(iBias) * ( self%binIndx(iBias) - self%binMin(iBias) )
      biasIndx = biasIndx + self%indexCoeff(iBias) * self%binIndx(iBias)
    enddo

    if(biasIndx > self%umbrellaLimit) then
      accept = .false.
      return
    endif

  end subroutine
!==========================================================================================
    subroutine UmbrellaWHAM_ReadInitialBias(self)
    implicit none
    class(UmbrellaWHAM), intent(inout) :: self
    integer :: AllocateStatus
    integer :: j, inStat, biasIndx, iUmbrella
    real(dp), allocatable :: varValue(:)
    real(dp) :: curBias, refVal

    open(unit=36, file=trim(adjustl(self%filename)) )

    allocate(varValue(1:self%nBiasVar), STAT = AllocateStatus )
    IF (AllocateStatus /= 0) error STOP "*** UmbrellaWHAM: Allocation Error ***"

    self%UBias = 0E0_dp
    biasIndx = 0
    do 
      read(36, *, IOSTAT=inStat) (varValue(j), j=1,self%nBiasVar), curBias

      if(inStat < 0) then
        exit
      endif

      call self%GetUIndexArray(varValue, biasIndx, inStat) 
      if(inStat == 1) then
        cycle
      endif

      self%UBias(biasIndx) = curBias
    enddo

    if(inStat < 0) then
      return
    endif

    

    refVal = self%UBias(self%refBin)
    do iUmbrella = 1, self%umbrellaLimit
      self%UBias(biasIndx) = self%UBias(biasIndx) - refVal
    enddo


    deallocate(varValue)

    close(36)

    end subroutine
!====================================================================
! This function is designed to take a set of input values for each biasing
! variable (varArray) and uses that information to generate a 1D-bias index U.
! Math: Given values (x1,x2,..xm) find the sub bins values (n1, n2,..nm) and find a U such that
! U - U0 = a1*n1 + a1*n2 + ... am*nm
  subroutine UmbrellaWHAM_GetUIndexArray(self, varArray, biasIndx, stat) 
    implicit none
    class(UmbrellaWHAM), intent(inout) :: self
    real(dp), intent(in) :: varArray(:)
    integer, intent(out) :: biasIndx, stat
    integer :: iBias
      

    stat = 0
    biasIndx = 1
    do iBias = 1, self%nBiasVar
      if(varArray(iBias) > self%valMax(iBias)) then
!      if(self%binIndx(iBias) > self%binMax(iBias)) then
        stat = 1
        return
      endif
      if(varArray(iBias) < self%valMin(iBias)) then
!      if(self%binIndx(iBias) < self%binMin(iBias)) then
        stat = -1
        return
      endif
!      binIndx(iBias) = floor( (varArray(iBias)-self%valMin(iBias)) / UBinSize(iBias) + 1E-8 )
      self%binIndx(iBias) = floor( (varArray(iBias) - self%valMin(iBias) ) / self%UBinSize(iBias) )
    enddo

    biasIndx = 1
    do iBias = 1, self%nBiasVar
      biasIndx = biasIndx + self%indexCoeff(iBias) * self%binIndx(iBias) 
    enddo


   end subroutine
!==========================================================================================
   subroutine UmbrellaWHAM_OutputUmbrellaHist(self)
     implicit none
     class(UmbrellaWHAM), intent(inout) :: self
     integer :: iUmbrella, iBias
     character(len = 100) :: outputString

     write(outputString, *) "(", ("F12.8, 2x", iBias =1,self%nBiasVar), "2x, F18.1)"
     open(unit=60, file="UmbrellaHist.txt")
      
     do iUmbrella = 1, self%umbrellaLimit
       call self%FindVarValues(iUmbrella, self%UArray)  
       do iBias = 1, self%nBiasVar        
         self%varValues(iBias) = real(self% UArray(iBias), dp) * self%UBinSize(iBias)          
       enddo

       if(self%UHist(iUmbrella) .ne. 0E0_dp) then
         write(60,outputString) (self%varValues(iBias), iBias=1,self%nBiasVar), self%UHist(iUmbrella)
       endif
        
     enddo 
    
     flush(60)
     close(60)
      
  end subroutine
!====================================================================
! This function is designed to take an input bin value for the biasing array
! (UIndx) and back converts it to give the bin indicies of the N-dimensional biasing
! array. 
! Math Problem: Given a value U find the values for (n1,n2,..nm) such that.
! U - U0 = a1*n1 + a1*n2 + ... am*nm
!
  subroutine UmbrellaWHAM_FindVarValues(self, UIndx, UArray) 
    implicit none 
    class(UmbrellaWHAM), intent(inout) :: self
    integer, intent(in) :: UIndx 
    integer, intent(inout) :: UArray(:) 
    integer :: i, iBias 
    integer :: remainder, curVal 
                            
    remainder = UIndx - 1 
    do i = 1, self%nBiasVar 
      iBias = self%nBiasVar - i + 1 
      curVal = int( real(remainder, dp)/real(self%indexCoeff(iBias),dp) ) 
!      self%UArray(iBias) = curVal + self%binMin(iBias) 
      self%UArray(iBias) = curVal 
      remainder = remainder - curVal * self%indexCoeff(iBias) 
    enddo 
                                                                             
  end subroutine
!====================================================================
  subroutine UmbrellaWHAM_ProcessIO(self, line, linestat) 
    use Input_Format, only: GetXCommand, maxLineLen
    implicit none
    class(UmbrellaWHAM), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line   
    integer, intent(out) :: lineStat

    integer :: i, intVal, intVal2
    real(dp) :: realVal
    character(len=30) :: command

    lineStat  = 0
    call GetXCommand(line, command, 3, lineStat)
    if(lineStat < 0) then
      return
    endif
    select case(trim(adjustl(command)))
      case("biasvariables")
        call GetXCommand(line, command, 4, lineStat)
        read(command, *) intVal
        self%nBiasVar = intVal
        if(.not. allocated(self%AnalysisIndex) ) then
          call self%Constructor
        endif

      case("analysis")
        if(allocated(self%AnalysisIndex) ) then
          call GetXCommand(line, command, 4, lineStat)
          read(command, *) intVal

          call GetXCommand(line, command, 5, lineStat)
          read(command, *) intVal2
 
          self%AnalysisIndex(intVal) = intVal2
        endif

      case("bounds")
        if(allocated(self%AnalysisIndex) ) then
          call GetXCommand(line, command, 4, lineStat)
          read(command, *) intVal

          call GetXCommand(line, command, 5, lineStat)
          read(command, *) intVal2
          self%nBins(intVal) = intVal2

          call GetXCommand(line, command, 6, lineStat)
          read(command, *) realVal
          self%valMin(intVal) = realVal

          call GetXCommand(line, command, 7, lineStat)
          read(command, *) realVal
          self%valMax(intVal) = realVal
          self%UBinSize(intVal) = (self%valMax(intVal)-self%valMin(intVal))/real(self%nBins(intVal), dp) 
!          self%binMin(intVal) = nint(realVal/self%UBinSize(intVal))
!          self%binMax(intVal) = nint(realVal/self%UBinSize(intVal))
        endif

      case("reference")
        do i = 1, self%nBiasVar
          call GetXCommand(line, command, 4+i-1, lineStat)
          read(command, *) realVal
          self%refVals(i) = realVal
        enddo

      case("whamfreq")
          call GetXCommand(line, command, 4, lineStat)
          read(command, *) realVal
          self%maintFreq = floor(realVal)

      case default
        lineStat = -1
    end select



   end subroutine
!====================================================================
  subroutine UmbrellaWHAM_Maintenance(self)
    use AnalysisData, only: AnalysisArray
    use ParallelVar, only: nout
    implicit none
    class(UmbrellaWHAM), intent(inout) :: self

    call self%Adjusthist
  end subroutine
!====================================================================
  subroutine UmbrellaWHAM_Epilogue(self)
    use AnalysisData, only: AnalysisArray
    use ParallelVar, only: nout, myid
    implicit none
    class(UmbrellaWHAM), intent(inout) :: self

    write(nout,*) "Writing Umbrella Sampling Histogram..."
  end subroutine
!====================================================================
  subroutine UmbrellaWHAM_UpdateStatistics(self, accept)
    use Template_SimBox, only: SimBox
    implicit none
    class(UmbrellaWHAM), intent(inout) :: self
    logical, intent(in) :: accept

    if(accept) then
      self%UHist(self%newIndx) = self%UHist(self%newIndx) + 1E0_dp
      self%oldIndx = self%newIndx
    else
!      self%oldIndx = self % GetBiasIndex()
      self%UHist(self%oldIndx) = self%UHist(self%oldIndx) + 1E0_dp
    endif


  end subroutine
!==================================================================================
  subroutine UmbrellaWHAM_WHAMSimOutput(self)
    use ParallelVar
    implicit none
    class(UmbrellaWHAM), intent(inout) :: self

    integer :: i,j
    integer :: UArray(1:self%nBiasVar)
    real(dp) :: varVal
    character(len = 100) :: outputString

    write(outputString, *) "(", ("2x, F10.4,", j =1,self%nBiasVar), "2x, F22.1)"

!        This block exports the histogram 
    rewind(self%histfile)
    do i = 1, self%umbrellaLimit
      if(self%HistStorage(i) /= 0E0_dp ) then
        call self%FindVarValues(i, self%UArray)
        write(self%histfile, outputString) (self%UArray(j)*self%UBinSize(j)+self%valMin(j), j=1,self%nBiasVar), self%HistStorage(i)
      endif
    enddo
    flush(self%histfile)

    write(outputString, *) "(", ("2x, F10.4,", j =1,self%nBiasVar), "2x, F18.8)"

!      This block exports the current umbrella bias
    rewind(self%potfile)
    do i = 1, self%umbrellaLimit
      call self%FindVarValues(i, self%UArray)
      write(self%potfile, outputString) (self%UArray(j)*self%UBinSize(j)+self%valMin(j),j=1,self%nBiasVar), self%UBias(i)
    enddo
    flush(self%potfile)


!        This block exports the calculated free energy to a file
    rewind(self%whamfile)
    do i = 1, self%umbrellaLimit
      if(self%ProbArray(i) /= 0E0_dp ) then
        call self%FindVarValues(i, self%UArray)
        write(self%whamfile, outputString) (self%UArray(j)*self%UBinSize(j)+self%valMin(j) ,j=1,self%nBiasVar), self%FreeEnergyEst(i)
      endif
    enddo
    flush(self%whamfile)



           
  end subroutine
!=========================================================================
!     This subroutine periodically adjusts the Umbrella Sampling Bias
!     by collecting histogram data from across each thread. 
  subroutine UmbrellaWHAM_AdjustHist(self)
    use ParallelVar, only: myid, nout, ierror
#ifdef MPIPARALLEL
    use MPI
#endif
    implicit none
    class(UmbrellaWHAM), intent(inout) :: self
    integer :: arraySize, i, j, cnt, maxbin, maxbin2
    real(dp) :: norm, maxBias
    real(qp) :: F_Estimate(1:self%nWhamItter), F_Old(1:self%nWhamItter), fSum, denomSum
    real(dp) :: tol, refBias

#ifdef MPIPARALLEL
    write(nout,*) "Halting for WHAM"
    call MPI_BARRIER(MPI_COMM_WORLD, ierror) 
#endif
    
!      This block condences the histogram data from all the different processors
!      into one collective array on the root (myid = 0) processor.        
#ifdef MPIPARALLEL
    arraySize = size(self%UHist)     
    if(myid .eq. 0) then
      self%TempHist = 0E0_dp
    endif
    call MPI_REDUCE(self%UHist, self%TempHist, arraySize, &
              MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierror)
#else
    self%TempHist = self%UHist
#endif

    

    if(myid .eq. 0) then
!        This block calculates the terms needed to 
      norm = sum(self%TempHist)
      do i = 1, self%umbrellaLimit
        self%BiasStorage(i, self%nCurWhamItter) = self%UBias(i) 
        if(self%TempHist(i) /= 0E0_dp) then
          self%WHAM_Numerator(i) = self%WHAM_Numerator(i) + self%TempHist(i)*self%TempHist(i)/norm
          self%WHAM_Denominator(i, self%nCurWhamItter) = self%TempHist(i)*exp(self%UBias(i))
          self%HistStorage(i) = self%HistStorage(i) + self%TempHist(i)
        endif
      enddo 

!        This block solves for the free energy terms required by WHAM.  This is done
!        self-consistently.
      self%ProbArray = 0E0_dp
      F_Estimate = 0E0_dp
      tol = self%tolLimit + 1E0_dp
      cnt = 0
      do while(tol > self%tolLimit)
        cnt = cnt + 1
!          Infinite Loop protection
        if(cnt > self%maxSelfConsist) then 
          exit
        endif

        self%ProbArray = 0E0_dp
        do j = 1, self%nCurWhamItter
          F_Old(j) = F_Estimate(j)
        enddo
!            If bin #i has been sampled at any point in the simulation, estimate the unbiased probability
!            based on the current guess value for F
        do i = 1, self%umbrellaLimit
          if(self%WHAM_Numerator(i) /= 0E0_dp) then
            denomSum = 0E0_dp
            do j = 1, self%nCurWhamItter
              if(self%WHAM_Denominator(i,j) > 0E0_dp) then
                denomSum = denomSum + self%WHAM_Denominator(i,j)*exp(-F_Estimate(j))
              endif
            enddo
            if(denomSum /= 0E0_dp) then
              self%ProbArray(i) = self%WHAM_Numerator(i)/denomSum
            endif
          else
            self%ProbArray(i) = 0E0_dp
          endif
        enddo 

        norm = sum(self%ProbArray)
        do i = 1, self%umbrellaLimit
          self%ProbArray(i) = self%ProbArray(i)/norm
        enddo 
!          Once all the unbiased probabilities have been estimated, use these unbiased estimates
!          to calculate a new estimate for F
        do j = 1, self%nCurWhamItter
          fSum = 0E0_dp
          do i = 1, self%umbrellaLimit
            if(self%ProbArray(i) /= 0E0_dp) then
              fSum = fSum + self%ProbArray(i)*exp(self%BiasStorage(i,j))
            endif
          enddo
          F_Estimate(j) = log(fSum)
          F_Estimate(j) = (F_Estimate(j) + F_Old(j))*0.5E0_dp
       enddo 
!         Calculate the average change in F from the previous estimate and determine 
!         if there has been a significant change to the F values.
       tol = 0E0_dp
       do j = 1, self%nCurWhamItter
         tol = tol + abs(F_Estimate(j) - F_Old(j))
       enddo
     enddo

!        Using the new estimates for the unbiased probability, calculate the free energy of the histogram
!        and modify the umbrella sampling bias such that it's equal to the new estimate.
      self%NewBias = 0E0_dp
      maxbin = maxloc(self%HistStorage,1)
      maxbin2 = maxloc(self%TempHist,1)
      write(nout,*) "Largest Bin", maxbin, maxbin2
      write(nout,*) "Largest Value", self%ProbArray(maxbin)
      if(self%ProbArray(maxbin) > 1E-300_dp) then
        do i = 1, self%umbrellaLimit
          if(self%ProbArray(i) > 0E0_dp) then
            self%FreeEnergyEst(i) = -log(self%ProbArray(i)/self%ProbArray(maxbin))
          endif
        enddo
      else
        self%FreeEnergyEst = 0E0_dp
      endif


      do i = 1, self%umbrellaLimit
        if(self%TempHist(i) < 1E0_dp) then
          self%NewBias(i) = self%UBias(i) - self%UBias(maxbin2) + log(self%TempHist(maxbin2))
        else
          self%NewBias(i) = self%UBias(i) - self%UBias(maxbin2) - log(self%TempHist(i)/self%TempHist(maxbin2))
        endif
      enddo
!        Rescale the pontential such that the reference free energy is set to 0
      refBias = self%NewBias(self%refBin)
      do i = 1, self%umbrellaLimit
        self%NewBias(i) = self%NewBias(i) - refBias
      enddo
      refBias = self%FreeEnergyEst(self%refBin)
      do i = 1, self%umbrellaLimit
        self%FreeEnergyEst(i) = self%FreeEnergyEst(i) - refBias
      enddo

    endif      !End of processor 0 only block


#ifdef MPIPARALLEL
    call MPI_BARRIER(MPI_COMM_WORLD, ierror) 
!      Distribute the new free energy estimate to all threads so that they can continue the simulation
!      with the new free energy. 
    arraySize = size(self%NewBias)      
    call MPI_BCast(self%NewBias, arraySize, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror) 
#endif

    do i = 1, self%umbrellaLimit
      self%UBias(i) = self%NewBias(i)
      self%UHist(i) = 0E0_dp
    enddo 
    if(myid .eq. 0) then
      call self%WHAMSimOutput
    endif

    self%nCurWhamItter = self%nCurWhamItter + 1
     
  end subroutine
!====================================================================
end module
!====================================================================
