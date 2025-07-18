!================================================================================
module Input_RegrowType
use Input_Format, only: LowerCaseLine, GetXCommand
use VarPrecision
!================================================================================
contains
!================================================================================
  subroutine Script_RegrowType(line, molNum, lineStat)
    use Common_MolInfo, only: MolData
    use ForcefieldData, only: nForceFields
    use MolCon_LinearCBMC, only: LinearCBMC
    use MolCon_SimpleRegrowth, only: SimpleRegrowth
    use MolCon_RidgidRegrowth, only: RidgidRegrowth
    use ParallelVar, only: nout
    implicit none
    character(len=*), intent(in) :: line
    integer, intent(in) :: MolNum
    integer, intent(out) :: lineStat

    character(len=30) :: dummy, command, Regrow_Type
!    character(len=30) :: fileName    
    logical :: logicValue
    integer :: j
    real(dp) :: realValue

    lineStat  = 0
    call GetXCommand(line, command, 2, lineStat)
    read(command, *) Regrow_Type

    !Safety check to ensure that the index number is within proper bounds
    select case(trim(adjustl(Regrow_Type)))
      case("ridgid")
        allocate(RidgidRegrowth :: MolData(molNum)%molConstruct )
        MolData(molNum)%ridgid = .true.
        write(nout,*) "Molecule uses ridgid regrowth"

      case("simple")
        allocate(SimpleRegrowth :: MolData(molNum)%molConstruct )
        MolData(molNum)%ridgid = .false.
        write(nout,*) "Molecule uses simple regrowth"

      case("linearcbmc")
        allocate(LinearCBMC :: MolData(molNum)%molConstruct )
        MolData(molNum)%ridgid = .false.
        write(nout,*) "Molecule uses Linear CBMC"

      case default
        lineStat = -1
        write(0,*) "Invalid Regrowth Type Specified in Molecule Definition!"
        write(0,*) "Molecule Type: ", molNum
        write(0,*) "Command:", line
        stop
    end select
!    call MolData(molNum)%molConstruct%Constructor(molNum)
    call MolData(molNum)%molConstruct%SetMolType(molNum)
    call MolData(molNum)%molConstruct%ProcessIO(line, linestat)

  end subroutine
!================================================================================
end module
!================================================================================
