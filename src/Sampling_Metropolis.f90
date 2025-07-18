!====================================================================
! Basic Metropolis Rule for molecular systems. 
!====================================================================
module MetropolisRule
  use VarPrecision
  use CoordinateTypes, only: Perturbation, Addition, Deletion, VolChange
  use Template_AcceptRule, only: acceptrule
 
  type, public, extends(acceptrule) :: metropolis
    contains
       procedure, pass :: MakeDecision => Metropolis_MakeDecision
       procedure, pass :: MakeDecision2Box => Metropolis_MakeDecision2Box
!       procedure, pass :: Maintenance => Metropolis_Maintenance
!       procedure, pass :: ProcessIO => Metropolis_ProcessIO
  end type
!====================================================================
  contains
!====================================================================
  function Metropolis_MakeDecision(self, trialBox, E_Diff, disp, inProb, logProb, extraIn) result(accept)
    use Template_SimBox, only: SimBox
    use ErrorChecking, only: IsNan, IsInf
    use RandomGen, only: grnd
    use ParallelVar, only: nout
    implicit none
    class(metropolis), intent(inout) :: self
    class(simBox), intent(in) :: trialBox
    class(Perturbation), intent(in) :: disp(:)
    real(dp), intent(in), optional:: inProb, logProb
    real(dp), intent(in), optional:: extraIn
    real(dp), intent(in) :: E_Diff
    logical :: accept
    integer :: iDisp
    real(dp) :: biasE = 1E0_dp
    real(dp) :: chemPot, extraTerms, probTerm


    if(IsNan(E_Diff)) then
      write(0,*) "ERROR! Invalid energy has been passed to the sampling routine!"
      write(0,*) E_Diff
      error stop
    endif

    accept = .false.
    if(present(inProb)) then
      if(inProb <= 0E0_dp) then
        return
!      write(nout,*) "Probability:", inProb
!      error stop "CRITICAL ERROR! Probability passed to the Metropolis Sampling Function is zero or negative!"
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
!    write(*,*) "Bias", E_diff
!    write(*,*) "    ", probTerm, extraTerms

    biasE = -trialBox%beta * E_Diff + probTerm + extraTerms
!    write(*,*) "    ",  biasE
!    write(*,*)
    if(biasE >= 0.0E0_dp) then
      accept = .true.
    elseif( biasE > log(grnd()) ) then
      accept = .true.
    endif

  end function
!====================================================================
  function Metropolis_MakeDecision2Box(self, trialBox1,  trialBox2, E_Diff1, E_Diff2, &
                                       disp1, disp2, inProb, &
                                       logProb, extraIn ) result(accept)
    use Template_SimBox, only: SimBox
    use RandomGen, only: grnd
    use ParallelVar, only: nout
    implicit none
    class(Metropolis), intent(inout) :: self
    class(SimBox), intent(in) :: trialBox1, trialBox2
    class(Perturbation), intent(in) :: disp1(:), disp2(:)
    real(dp), intent(in) :: E_Diff1, E_Diff2
    real(dp), intent(in), optional :: inProb, logProb, extraIn
    logical :: accept
    integer :: iDisp
    real(dp) :: biasE, chemPot, extraTerms, probTerm



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

!    write(*,*) "Bias", E_diff1, E_Diff2
!    write(*,*) "    ", probTerm, extraTerms

    biasE = -trialBox1%beta*E_Diff1 - trialBox2%beta*E_Diff2 + probTerm + extraTerms
!    write(*,*) "    ",  biasE
!    write(*,*)
    if(biasE > 0.0E0_dp) then
      accept = .true.
    elseif( biasE > log(grnd()) ) then
      accept = .true.
    endif

  end function

!====================================================================
end module
!====================================================================
