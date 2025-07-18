!========================================================            
#define __StdErr__ 0
!========================================================            
    module ScriptInput
      use Input_Format, only: LoadFile,GetXCommand, LowerCaseLine, FindCommandBlock, maxLineLen
    contains
!========================================================            
      subroutine Script_ReadParameters(infile)
      use ClassyConstants
      use Input_Forcefield
      use Input_AnalysisType, only: Script_AnalysisType
      use Input_Sampling, only: Script_SamplingType
      use Input_NeighType, only: Script_NeighType
      use Input_Initialize, only: Script_Initialize
#ifdef EMBPYTHON
      use SimPython, only: RunPythonTest
#endif
      use SimControl, only: TimeStart, TimeEnd
      use SimMonteCarlo, only: RunMonteCarlo
      use SimMinimize, only: RunMinimize
      use ParallelVar
      use Units
      implicit none
      integer :: i, ii, j, nArgs, boxnum
      integer :: iLine, lineStat, AllocateStat
      integer :: nLines, nForceLines, lineBuffer
      integer, allocatable :: lineNumber(:)
      character(len=*), intent(in), optional :: infile


      character(len=maxLineLen), allocatable :: lineStore(:)
      character(len=30) :: command, command2, dummy
      character(len=50) :: fileName
      character(len=50) :: forcefieldFile

    
      nLines = 0
      if(present(infile)) then
         filename = trim(adjustl(infile))
         call LoadFile(lineStore, nLines, lineNumber, fileName)
      else
         !If no filename was passed, get the filename from the command line. 
          nArgs = command_argument_count()
          if(nArgs > 0) then
            call get_command_argument(1, fileName)
            call LoadFile(lineStore, nLines, lineNumber, fileName)
            if(nLines == 0) then
              write(__StdErr__,*) "ERROR! Input file is empty or could not be read!"
              stop
            else
              write(nout, *) "File successfully loaded!"
            endif
          elseif(nArgs == 0) then
            write(__StdErr__,*) "ERROR! No Input File has been given!"
            stop
          endif
      endif

!      This block counts the number of lines in the input file to determine how large the lineStorage array needs to be.

      do iLine = 1, nLines
        call LowerCaseLine(lineStore(iLine))
      enddo
!      call setDefaults(seed, screenEcho)

      lineBuffer = 0
      do iLine = 1, nLines
        if(lineBuffer > 0) then
          lineBuffer = lineBuffer - 1
          cycle
        endif
        lineStat = 0        
        call GetCommand(lineStore(iLine), command, lineStat)
!         If line is empty or commented, move to the next line.         
        if(lineStat .eq. 1) then
          cycle
        endif 

        select case(trim(adjustl( command )))
          case("create")
            call createCommand(iLine, linestore, lineBuffer, lineStat)

          case("forcefield")
            call GetXCommand(lineStore(iLine), filename, 2, lineStat)  
            call Script_ReadFieldFile(filename, lineStat)

          case("modify")
            call modifyCommand( lineStore(iLine), lineStat )

          case("neighlist")
            call Script_NeighType( lineStore(iLine), lineStat)

          case("samplingtype")
            call Script_SamplingType(iLine, lineStore, lineStat)

          case("set")
            call setCommand( lineStore(iLine), lineStat )

          case("minimize")

            call CPU_TIME(TimeStart)
            call Script_Initialize
            call GetXCommand(lineStore(iLine), command2, 2, lineStat)  
            read(command2,*) boxnum
            call RunMinimize(boxnum)
            call CPU_TIME(TimeEnd)


          case("run")

            call CPU_TIME(TimeStart)
            call Script_Initialize
            call RunMonteCarlo
            call CPU_TIME(TimeEnd)
#ifdef EMBPYTHON
          case("testpython")

            call CPU_TIME(TimeStart)
            call GetXCommand(lineStore(iLine), filename, 2, lineStat)  
            call Script_Initialize
            call RunPythonTest(filename)
            call CPU_TIME(TimeEnd)
#endif
           
          case default
            write(__StdErr__,"(A,2x,I10)") "ERROR! Unknown Command on Line", lineNumber(iLine)
            write(__StdErr__,*) trim(adjustl(lineStore(iLine)))
            error stop 
        end select

        ! Ensure that the called processes exited properly.
        if(lineStat .eq. -1) then
          write(__StdErr__,"(A,1x,I10)") "ERROR! Command on line:", lineNumber(iLine)
          write(__StdErr__, "(A)") "could not be understood. Please check command for accuracy and try again."
          write(__StdErr__,*) trim(adjustl(lineStore(iLine)))
          error stop 
        endif
        
      enddo
!      write(*,*) "Finished Reading Input Script."      
      if(allocated(lineStore)) then
        deallocate(lineStore)
      endif


      end subroutine
!========================================================            
      subroutine setCommand(line, lineStat)
      use Common_MolInfo
      use Common_NeighData
      use ParallelVar
      use RandomGen, only: initSeed
      use SimControl, only: nMoves, nCycles, screenFreq, energyCheck, &
                            Etol, Forcetol, lrate, configfreq
      use Units, only: outEngUnit, outLenUnit, outAngUnit,outPressUnit,  &
                       inEngUnit, inLenUnit, inAngUnit,inPressUnit, &
                       FindEngUnit, FindLengthUnit, FindAngularUnit,  &
                       FindPressureUnit, engStr
      use VarPrecision
      implicit none
      character(len=maxLineLen), intent(in) :: line      
      integer, intent(out) :: lineStat

      character(len=50) :: command, command2
      logical :: logicValue
      integer :: intValue
      real(dp) :: realValue
      


      lineStat  = 0

      call GetXCommand(line, command, 2, lineStat)
      select case(trim(adjustl(command)))
        case("cycles")
          call GetXCommand(line, command, 3, lineStat)
          read(command, *) realValue  
          nCycles = nint(realValue, kind=8)

        case("energycheck")
          call GetXCommand(line, command, 3, lineStat)
          read(command, *) intValue  
          energyCheck = intValue

        case("moves")
          call GetXCommand(line, command, 3, lineStat)
          read(command, *) realValue  
          nMoves = nint(realValue, kind=8)

        case("screenecho")
          call GetXCommand(line, command, 3, lineStat)
          read(command, *) logicValue
          if(logicValue) then
            if(myid == 0) then
              nout = 6
            endif
          else
            if(myid == 0) then
              nout = 100
            endif
          endif

        case("rng_seed")
          call GetXCommand(line, command2, 3, lineStat)
          read(command2,*) intValue
          initSeed = intValue

        case("angleunits")
          call GetXCommand(line, command2, 3, lineStat)
          outAngUnit = FindAngularUnit(command2)

        case("distunits")
          call GetXCommand(line, command2, 3, lineStat)
          outLenUnit = FindLengthUnit(command2)

        case("energyunits")
          call GetXCommand(line, command2, 3, lineStat)
          outEngUnit = FindEngUnit(command2)
          engStr = trim(adjustl(command2))

        case("energytol")
          call GetXCommand(line, command2, 3, lineStat)
          read(command2, *) realValue
          Etol = realValue

        case("forcetol")
          call GetXCommand(line, command2, 3, lineStat)
          read(command2, *) realValue
          forcetol = realValue

        case("learnrate")
          call GetXCommand(line, command2, 3, lineStat)
          read(command2, *) realValue
          lrate = realValue

        case("pressureunits")
          call GetXCommand(line, command2, 3, lineStat)
          inPressUnit = FindPressureUnit(command2)
          outPressUnit = FindPressureUnit(command2)

        case("neighskin")
          call GetXCommand(line, command2, 3, lineStat)
          read(command2, *) realValue
          neighSkin = realValue

        case("configfrequency")
          call GetXCommand(line, command2, 3, lineStat)
          read(command2, *) realValue
          configFreq = nint(realValue)

        case("screenfrequency")
          call GetXCommand(line, command2, 3, lineStat)
          read(command2, *) realValue
          screenFreq = nint(realValue)

        case default
          lineStat = -1
      end select

     
      end subroutine
!========================================================            
      subroutine createCommand(iLine, linestore, lineBuffer, lineStat)
      use AnalysisData, only: AnalysisArray, analyCommon
      use BoxData, only: BoxArray
      use BoxPresets, only: Preset
      use TrajData, only: TrajArray
      use MCMoveData, only: Moves, MoveProb
      use ForcefieldData, only: EnergyCalculator, nForceFields

      use Input_SimBoxes, only: Script_BoxType
      use Input_Forcefield, only: Script_FieldType
      use Input_AnalysisType, only: Script_AnalysisType
      use Input_Constraint, only: Script_Constraint
      use Input_TrajType, only: Script_TrajType
      use Input_Moves, only: Script_MCMoves
      use Input_LoadCoords, only: Script_ReadCoordFile

      use VarPrecision
      use Units
      implicit none
      integer, intent(in) :: iLine
      character(len=maxLineLen), intent(in) :: linestore(:) 
      integer, intent(out) :: lineStat, lineBuffer

      character(len=30) :: dummy, command, command2
      character(len=50) :: fileName
      logical :: logicValue
      integer :: i, curLine
      integer :: intValue, AllocateStat, nItems
      real(dp) :: realValue
      

      lineStat  = 0
      AllocateStat = 0

      read(linestore(iLine),*) dummy, command
      call FindCommandBlock(iLine, lineStore, "end_create", lineBuffer)
      nItems = lineBuffer - 1
!      write(*,*) nItems

      select case(trim(adjustl(command)))
        case("analysis") 
           if( .not. allocated(AnalysisArray) ) then
             allocate(AnalysisArray(1:nItems), stat = AllocateStat)
             allocate(analyCommon(1:nItems), stat = AllocateStat)
             do i = 1, nItems
               curLine = iLine + i
               call Script_AnalysisType(linestore(curLine), i, lineStat)
             enddo   
             do i = 1, nItems
               call AnalysisArray(i) % func % CastCommonType(analyCommon(i)%val)
             enddo   

           else
             write(0,*) "ERROR! The create analysis command has already been used and can not be called twice"
             stop
           endif

         !-------------------------------------------------------------------------------------
        case("boxes")
           if( .not. allocated(BoxArray) ) then
             allocate(BoxArray(1:nItems), stat = AllocateStat)
             do i = 1, nItems
               curLine = iLine + i
               call GetXCommand(lineStore(curLine), command2, 1, lineStat)
               if( trim(adjustl(command2)) == "fromfile") then
                   call GetXCommand(lineStore(curLine), command2, 2, lineStat)
                   fileName = ""
                   fileName(1:30) = command2(1:30)
                   call Script_ReadCoordFile(fileName, i, lineStat)

               elseif( trim(adjustl(command2)) == "preset") then
                   call Preset(curLine, lineStore, i, lineStat)

               else
                   call Script_BoxType(linestore(curLine), i, lineStat)
               endif
               BoxArray(i)%box%boxID = i 
               BoxArray(i)%box%screenIO = .true.
             enddo             
           else
             write(0,*) "ERROR! The create box command has already been used and can not be called twice"
             stop
           endif

         !-------------------------------------------------------------------------------------
        case("constraint")
           if( .not. allocated(BoxArray) ) then
             error stop "Box array not allocated!"
           endif
           call GetXCommand(lineStore(iLine), command2, 3, lineStat)
           read(command2,*) intValue
           call Script_Constraint(lineStore, iLine, intValue, lineBuffer, lineStat)


!        case("energyfunctions")
!           if( .not. allocated(EnergyCalculator) ) then
!             nForceFields = nItems
!             allocate(EnergyCalculator(1:nItems), stat = AllocateStat)
!             do i = 1, nItems
!               curLine = iLine + i
!               call Script_FieldType(linestore(curLine), i, lineStat)
!               call EnergyCalculator(i)%Method%Constructor
!             enddo   
!           else
!             write(*,*) "ERROR! The create energycalculators command has already been used and can not be called twice"
!             stop
!           endif
!
         !-------------------------------------------------------------------------------------
        case("moves") 
           if( .not. allocated(Moves) ) then
             nForceFields = nItems
             allocate(Moves(1:nItems), stat = AllocateStat)
             allocate(MoveProb(1:nItems), stat = AllocateStat)
             do i = 1, nItems
               curLine = iLine + i
               call Script_MCMoves(linestore(curLine), i, lineStat)
             enddo   
           else
             write(*,*) "ERROR! The create moves command has already been used and can not be called twice"
             error stop
           endif

         !-------------------------------------------------------------------------------------
        case("trajectory") 
           if( .not. allocated(TrajArray) ) then
             allocate(TrajArray(1:nItems), stat = AllocateStat)
             do i = 1, nItems
               curLine = iLine + i
               call Script_TrajType(linestore(curLine), i, lineStat)
             enddo   
           else
             write(*,*) "ERROR! The create trajectory command has already been used and can not be called twice"
             error stop
           endif
 
        case default
          write(*,*) command
          lineStat = -1
      end select

      IF (AllocateStat /= 0) then
        write(*,*) AllocateStat
        error STOP "Allocation Error in Create Command"
      endif
     
      end subroutine   
!========================================================            
    subroutine modifyCommand(line, lineStat)
      use AnalysisData, only: AnalysisArray, analyCommon
      use BoxData, only: BoxArray
      use Template_SimBox
      use ForcefieldData, only: EnergyCalculator, nForceFields
      use MCMoveData, only: Moves
      use CommonSampling, only: sampling
      use VarPrecision
      use Units
      implicit none
      character(len=maxLineLen), intent(in) :: line      
      integer, intent(out) :: lineStat

      character(len=30) :: dummy, command, command2
      logical :: logicValue
      integer :: intValue, AllocateStat
      real(dp) :: realValue
      

      lineStat  = 0
      AllocateStat = 0
!      read(line,*) dummy, command
      call GetXCommand(line, command, 2, lineStat)
      call LowerCaseLine(command)
      select case(adjustl(trim(command)))
        case("analysis")
           if(.not. allocated(AnalysisArray)) then
             write(__StdErr__,*) "ERROR! You are trying to modify a analysis object before they have been defined!"
             linestat = -1
             return
           endif

           call GetXCommand(line, command2, 3, lineStat)
           read(command2, *) intValue

           if(intvalue > size(AnalysisArray)) then
             write(__StdErr__,*) "Invalid object number given to modify."
             write(__StdErr__,*) "Analysis number does not exit!"
             write(__StdErr__,*) line
             return
           endif
           call AnalysisArray(intValue) % func % ModifyIO(line, lineStat)

        case("box")
           if(.not. allocated(BoxArray)) then
             write(__StdErr__,*) "ERROR! You are trying to modify a box before boxes have been defined!"
             linestat = -1
             return
           endif

           call GetXCommand(line, command2, 3, lineStat)
           read(command2, *) intValue
           if(intvalue > size(BoxArray) .or. intvalue < 1) then
             write(__StdErr__,*) "Invalid object number given to modify."
             write(__StdErr__,*) "Box Object does not exist!"
             write(__StdErr__,*) line
             linestat = -1
             return
           endif
           call BoxArray(intValue) % box % ProcessIO(line, lineStat)

        case("move")
           if(.not. allocated(Moves)) then
             write(__StdErr__,*) "ERROR! You are trying to modify a move before they have been defined!"
             linestat = -1
             return
           endif

           call GetXCommand(line, command2, 3, lineStat)
           read(command2, *) intValue
           if(intvalue > size(Moves) .or. intvalue < 1) then
             write(__StdErr__,*) "Invalid object number given to modify."
             write(__StdErr__,*) "Move Number does not exist!."
             write(__StdErr__,*) line
             return
           endif
           call Moves(intValue) % Move % ProcessIO(line, lineStat)


        case("sampling")
           if(.not. allocated(Sampling)) then
             write(__StdErr__,*) "ERROR! You are trying to modify the sampling object before it has been defined!"
             linestat = -1
             return
           endif
           call Sampling % ProcessIO(line, lineStat)

        case default
           lineStat = -1
      end select

      IF (AllocateStat /= 0) error STOP "CRITICAL ERROR! Allocation Error in the Modify Command"
     
      end subroutine
!========================================================            
      subroutine IO_ErrorControl(iLine, lineNumber, lineStore, lineStat)
      implicit none
      integer, intent(in) :: iLine, lineStat
      integer, intent(in) :: lineNumber(:)
      character(len=maxLineLen), intent(in) :: lineStore(:)

      if(lineStat .eq. -1) then
        write(*,"(A,2x,I10)") "ERROR! Unknown Variable Name on Line", lineNumber(iLine)
        write(*,*) lineStore(iLine)
        error stop 
      endif

      end subroutine
!========================================================            
      end module
