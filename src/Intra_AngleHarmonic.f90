!=============================================================================+
module IntraAngle_Harmonic
  use Template_IntraAngle, only: Angle_FF
  use VarPrecision
  use Template_SimBox, only: SimBox
  use CoordinateTypes
  use Units, only: inAngUnit, inEngUnit

  type, public, extends(Angle_FF) :: HarmonicAngle
!    real(dp) :: theta0
    real(dp), private :: k0
    contains
      procedure, pass :: EFunc => HarmonicAngle_EFunc
      procedure, pass :: Constructor => HarmonicAngle_Constructor
!      procedure, pass :: DetailedECalc => HarmonicAngle_DetailedECalc
      procedure, pass :: GenerateDist => HarmonicAngle_GenerateDist
      procedure, pass :: GenerateTrial => HarmonicAngle_GenerateTrial
!      procedure, pass :: GenerateReverseDist => HarmonicAngle_GenerateReverseDist
      procedure, pass :: ProcessIO => HarmonicAngle_ProcessIO
  end type

  contains
!=============================================================================+
  function HarmonicAngle_EFunc(self, angle) result(E_Angle)
    implicit none
    class(HarmonicAngle), intent(in) :: self
    real(dp), intent(in) :: angle
    real(dp) :: E_Angle

    E_Angle = 0.5E0_dp * self%k0 * (angle - self%theta0)**2
  end function
!=============================================================================+
  subroutine HarmonicAngle_Constructor(self)
    use Common_MolInfo, only: nMolTypes
    implicit none
    class(HarmonicAngle), intent(inout) :: self


  end subroutine
!==========================================================================
  subroutine HarmonicAngle_GenerateTrial(self, beta, val, bounds )
    use RandomGen, only: grnd, Gaussian
    use ClassyConstants, only: two_pi
    implicit none
    class(HarmonicAngle), intent(inout) :: self
    real(dp), intent(in) :: beta
    real(dp), intent(out) :: val
    real(dp), intent(in), optional :: bounds(1:2)
    real(dp) :: EGen
    real(dp) :: lb, ub

    if(.not. present(bounds)) then
      lb = 0E0_dp
      ub = two_pi
    else
      lb = bounds(1)
      ub = bounds(2)
    endif

    do
      val = Gaussian()
      val = self%theta0 + val/sqrt(self%k0*beta)
      if( (val > ub) .or. (val < lb) ) cycle
      exit
    enddo
  end subroutine
!==========================================================================
  subroutine HarmonicAngle_GenerateDist(self, beta, val, probgen, E_T)
    use ClassyConstants, only: two_pi
    use RandomGen, only: grnd
    implicit none
    class(HarmonicAngle), intent(inout) :: self
    real(dp), intent(in) :: beta
    real(dp), intent(out) :: val
    real(dp), intent(out) :: probgen
    real(dp), intent(out), optional :: E_T
    real(dp) :: EGen

    do
!      val = Gaussian()
!      val = self%theta0 + val/sqrt(self%k0*beta)
!      if( (val > two_pi) .or. (val < 0.0E0_dp) ) then
!        cycle
!      endif
      call self%GenerateTrial(beta, val)
      probgen = sin(val)
      if(probgen > grnd()) then
        exit
      endif
    enddo
    
    EGen = self%EFunc(val)
    probgen = exp(-beta*EGen)

    if(present(E_T)) then     
      E_T = EGen
    endif

  end subroutine
!=============================================================================+
  subroutine HarmonicAngle_ProcessIO(self, line)
    use Input_Format, only: maxLineLen, GetXCommand 
    implicit none
    class(HarmonicAngle), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line
    character(len=30) :: command
    integer :: lineStat

    call GetXCommand(line, command, 2, lineStat)
    read(command, *) self%theta0
    self%theta0 = self%theta0*inAngUnit

    call GetXCommand(line, command, 3, lineStat)
    read(command, *) self%k0
    self%k0 = self%k0*inEngUnit

!    write(*,*) self%theta0, self%k0

    if(lineStat /= 0) then
      write(*,*) "Missing input required for the haromic Angle style"
    endif
  end subroutine
!=============================================================================+
end module
!=============================================================================+
