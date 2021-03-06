!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Originally created by James McDonagh, Stuart Davie and Rob Coates (Contributions from Dale Stuchfield)!
! at the University of Manchester 2015 in the Popelier group.                                           !
! Acknowledge the use of strip spaces from Jauch 27 November and sort_asend from Rossetta code          !
! Licensed under MIT License                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! File written by JMcD 3/4/16                                                                           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! A set of program and routines cobined into one program. This is master source which controls execution 
! of the rest. A set of python and R scripts are provided for pre and post processing of some of the data
! TASK 1 : Determines the structural deviations over a DL POLY MD trajectory or normal mode distortions 
!          from Tyche. This is sampling statistcs and information. Pre-processing required HISTORY-XYZ.py
!          or prep_Tychegjf_for_structure.py. Post-processing MDplot.R over the directories of csv files
!          that are output from the program
! TASK 2 : Cut out a micro-solvation shell around one atom from a DL POLY MD trajectory. Produces an RDF
!          From the atom of interest and gives a stacked XYZ files for visual tools as output. 
!          Post-processing XYZ-to-GJF-sort.py to turn the stacked XYZ file in to a set of GJFs ordered by
!          Total number of surrounding water molecules. Limited to water solvent.
! TASK 3 : Cut out the first solvation shell around one atom from a DL POLY MD trajectory. Produces an 
!          RDF from each solute atom of interest and gives a stacked XYZ files for visual tools as  
!          output. Post-processing XYZ-to-GJF-sort.py to turn the stacked XYZ file in to a set of GJFs 
!          ordered by. Total number of surrounding water molecules. Limited to water solvent.
! TASK 4 : Calculate prediction and model statistics from a Ferebus Kriging task. These are useful to 
!          evaluate how good of model is made. This can be done for IQA (self and interaction) and 
!          multipole models from FEREBUS kriging engine.
! TASK 5 : Place nodes at high solvent density position in a system. Order the features such that they 
!          are labelled relative to a fixed node position not labelled as free distingushable molecules.  
!          Reduce the features by analysing the imporatnce. Cut doen the feature list, and keep a backup  
!          of the original feature list. Cut out the first solvation shell around one atom from a
!          DL POLY MD trajectory.
!
! TASK 1 JMcD, TASK 2 and TASK 3 R.C., JMcD and S.D. TASK 4 S.D.(D.S. contributed to discussion of TASK1) 
! TASK 5 S.D (contributions from JMcD)
! 
! GitHub website https://github.com/Jammyzx1/Hermes/tree/FINAL
! GitHub alternative website https://github.com/Jammyzx1/Hermes
!
! Version 1.3
! CHANGE LOG 
! Version 1  : Separate programs with some interlinking pre and post processing
! Version 1.1: Modular format and incorporation in the Hermes program.
! Version 1.2: Added functionality for reduction and sorting of training data.
! Version 1.3: Stuart has added the ability for the user to enable a fixed number of solvent molecules 
!              around a solute
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Program Hermes

  use structurechanges
  use HISTORYRDFmicro
  use HISTORYrdf
  use Stats
  use krigingcoefficientiqa
  use den_nodes
  use wat_order
  use analysis_of_varience
  use nodorder
  use solorder
  
  implicit none
  
  integer              :: ierr, narg, i, j, k, l, tasklen, noptions, AtomsBefore, NoWater, o_lap, natoms, Nbins 
  integer              :: istatus, MinWaters
  double precision     :: res, fvt
  character (len=100)  :: fname, dummy, dummy2, dummy3, dummy4, strucfname, solutefname, connect, strucatom
  character (len=300)  :: buffer 
  character (len=9)    :: dum, task, Omicro, solutemicro, chargemicro, atomnummicro
  character (len=72)   :: line
  character (len=50)   :: TSname
  logical              :: file_ex, file_ex2

  print*,'------------------------------------ HERMES Version 1.2 ----------------------'
  print*,' TASK 1 : Determines the structural deviations over a DL POLY MD trajectory or' 
  print*,'          normal mode distortions from Tyche. This is sampling statistcs and '
  print*,'          information. Pre-processing required HISTORY-XYZ.py'
  print*,'          or prep_Tychegjf_for_structure.py. Post-processing MDplot.R over the'
  print*,'          directories of csv files that are output from the program'
  print*,' TASK 2 : Cut out a micro-solvation shell around one atom from a DL POLY MD' 
  print*,'          trajectory. Produces an RDF From the atom of interest and gives a '
  print*,'          stacked XYZ files for visual tools as output. Post-processing '
  print*,'          XYZ-to-GJF-sort.py to turn the stacked XYZ file in to a set of GJFs'
  print*,'          ordered by. Total number of surrounding water molecules. Limited to'
  print*,'          water solvent.'
  print*,' TASK 3 : Cut out the first solvation shell around one atom from a DL POLY ' 
  print*,'          MD trajectory. Produces an RDF from each solute atom of interest '  
  print*,'          and gives a stacked XYZ files for visual tools as output. ' 
  print*,'          Post-processing XYZ-to-GJF-sort.py to turn the stacked XYZ file'
  print*,'          in to a set of GJFs ordered by. Total number of surrounding water'
  print*,'          molecules. Limited to water solvent.'
  print*,' TASK 4 : Calculate prediction and model statistics from a Ferebus Kriging ' 
  print*,'          task. These are useful to evaluate how good of model is made. '
  print*,'          This can be done for IQA (self and interaction) and multipole'
  print*,'          models from FEREBUS kriging engine.'
  print*,' TASK 5 : Place nodes at high solvent density position in a system. Order the' 
  print*,'          features such that they are labelled relative to a fixed node '  
  print*,'          position not labelled as free distingushable molecules. Reduce the'    
  print*,'          features by analysing the imporatnce. Cut down the feature list,' 
  print*,'          and keep a backup of the original feature list.' 
  print*,'------------------------------------------------------------------------------'

  open (status="replace", unit=500, file="Hermes.log", action="write",iostat=ierr)                                ! Unit 40 new output file BONDS
  if (ierr.ne.0) then 
     stop 'ERROR - Unable to create Hermes.log.'
  else  
     print*,'SUCCESS - Hermes.log has been created'
  end if

  write(500,'(A)'), 'Hermes Version 1.1, authors J. L. McDonagh, S. Davie and R. Coates (contribution from Dale Stuchfield) '
  write(500,'(A)'), 'Acknowledge the use of adapted routines : '
  write(500,'(A)'), 'strip spaces from Jauch 27 November'
  write(500,'(A)'), 'sort_asend from Rossetta code'
  write(500,'(A)')'------------------------------------ HERMES Version 1.2 -----------------------------------------------'
  write(500,'(A)')' TASK 1 : Determines the structural deviations over a DL POLY MD trajectory or normal mode distortions '
  write(500,'(A)')'          from Tyche. This is sampling statistcs and information. Pre-processing required HISTORY-XYZ.py'
  write(500,'(A)')'          or prep_Tychegjf_for_structure.py. Post-processing MDplot.R over the directories of csv files'
  write(500,'(A)')'          that are output from the program'
  write(500,'(A)')' TASK 2 : Cut out a micro-solvation shell around one atom from a DL POLY MD trajectory. Produces an RDF'
  write(500,'(A)')'          From the atom of interest and gives a stacked XYZ files for visual tools as output.'
  write(500,'(A)')'          Post-processing XYZ-to-GJF-sort.py to turn the stacked XYZ file in to a set of GJFs ordered by'
  write(500,'(A)')'          Total number of surrounding water molecules. Limited to water solvent.'
  write(500,'(A)')' TASK 3 : Cut out the first solvation shell around one atom from a DL POLY MD trajectory. Produces an' 
  write(500,'(A)')'          RDF from each solute atom of interest and gives a stacked XYZ files for visual tools as'  
  write(500,'(A)')'          output. Post-processing XYZ-to-GJF-sort.py to turn the stacked XYZ file in to a set of GJFs' 
  write(500,'(A)')'          ordered by. Total number of surrounding water molecules. Limited to water solvent.'
  write(500,'(A)')' TASK 4 : Calculate prediction and model statistics from a Ferebus Kriging task. These are useful to' 
  write(500,'(A)')'          evaluate how good of model is made. This can be done for IQA (self and interaction) and'
  write(500,'(A)')'          multipole models from FEREBUS kriging engine.'
  write(500,'(A)')' TASK 5 : Place nodes at high solvent density position in a system. Order the features such that ' 
  write(500,'(A)')'          they are labelled relative to a fixed node position not labelled as free distingushable '  
  write(500,'(A)')'          molecules. Reduce the features by analysing the imporatnce. Cut down the feature list, and'    
  write(500,'(A)')'          keep a backup of the original feature list.' 
  write(500,'(A)')'-------------------------------------------------------------------------------------------------------'
!
!----------------------------------------- Initalisation -----------------------------------------------------------------
!
 i = 0
 j = 0
 k = 0
 l = 0
 narg = 0
 noptions = 0
! 
!----------------------------------------- Get input ---------------------------------------------------------------------
! 
 narg = IARGC()
 i = 1

 if(narg.lt.1) then
    do while (i==1)
       print*,'INTERACTIVE'
       write(*,'(A)', ADVANCE="NO") 'Please enter the input file name - ' 
       read*,fname
       write(500,'(A)') 'Filename ', fname
       fname = trim(adjustl(fname))
       open (status="old", unit=700, file=fname, action="read", iostat=ierr)                             ! Unit 700  input file
       if (ierr.ne.0) then 
          stop 'ERROR - Unable to open input file.'
       else  
          print*,'SUCCESS - input file opened'
          i = 0
       end if
    end do

 else
    call GETARG(1,fname)
    write(500,'(A)') 'Filename ', fname
    open (status="old", unit=700, file=fname, action="read",iostat=ierr)                                ! Unit 700  input file
    if (ierr.ne.0) then 
       stop 'ERROR - Unable to open input file.'
    else  
       print*,'SUCCESS - input file opened'
    end if
 end if
! 
!----------------------------------------- read input run the appropriate calculations -------------------------------
!
 do while(ierr.eq.0)
    read(700,'(A)') line
    if(index(line,'task').ne.0) then
       print*, 'task if loop entered'
       read(line,*), dum, task
       print*,'dummy ', dum
       print*,'Task ', task
       ierr = 1
    else
       print*, 'ERROR: Task was not given on the first input line of input file please check. Input file - ', fname
    end if
 end do
 write(500,'(A)')'Task ', trim(adjustl(task))
 rewind 700
 ierr = 0
 tasklen = len_trim(adjustl(trim(task)))
!
!- - - - - - - - - Task 1 - - - - - - - - - 
!
 if (index(task,'structure').ne.0) then
    do while(ierr.eq.0)
        read(700,'(A)') line
        if(index(line,'filename').ne.0) then         ! File of bonding, angle and dihedral information or auto to make it itself
            read(line,*), dummy, strucfname
            strucfname=trim(adjustl(strucfname))
            noptions = noptions + 1                  ! The noptions variable counts the number of valid options found in the input file too few or too many gives an error
        else if(index(line,'atom').ne.0) then        ! The atom label to centre the viewing file on
            read(line,*), dummy, strucatom
            strucatom=trim(adjustl(strucatom))
            noptions = noptions + 1
        else if(index(line,'END').ne.0) then
            print*, 'Input file read'
            write(500,'(A)') 'Input file read'
            if (noptions.lt.2) then
               write(500,'(A)') 'ERROR: Too few arguments passed for task structure', noptions
               stop 'ERROR: Too few arguments passed for task structure'
            else if (noptions.gt.2) then
               write(500,'(A)') 'ERROR: Too many arguments passed for task structure', noptions
               stop 'ERROR: Too many arguments passed for task structure'
            end if
            ierr=1
        end if
    end do
    call stucture(strucfname, strucatom)
    write(500,'(A)') 'Structural analysis ended normally'
!
 else if (index(task,'1').ne.0) then
    do while(ierr.eq.0)
        read(700,'(A)') line
        if(index(line,'filename').ne.0) then         ! File of bonding, angle and dihedral information or auto to make it itself
            read(line,*), dummy, strucfname
            strucfname=trim(adjustl(strucfname))
            noptions = noptions + 1
        else if(index(line,'atom').ne.0) then        ! The atom label to centre the viewing file on
            read(line,*), dummy, strucatom
            strucatom=trim(adjustl(strucatom))
            noptions = noptions + 1
        else if(index(line,'END').ne.0) then
            print*, 'Input file read'
            write(500,'(A)') 'Input file read'
            if (noptions.lt.2) then
               write(500,'(A)') 'ERROR: Too few arguments passed for task 1'
               stop 'ERROR: Too few arguments passed for task 1'
            else if (noptions.gt.2) then 
               write(500,'(A)') 'ERROR: Too many arguments passed for task 1'
               stop 'ERROR: Too many arguments passed for task 1'
            end if
            ierr=1
        end if
    end do
    call stucture(strucfname, strucatom)
    write(500,'(A)') 'Structural analysis ended normally'
!
!- - - - - - - - - Task 2 - - - - - - - - - 
!    
 else if(index(task,'microsolv').ne.0) then
     do while(ierr.eq.0)
        read(700,'(A)') line
        if(index(line,'Otype').ne.0) then            ! Water O label eg OTP
            read(line,*), dummy, Omicro
            noptions = noptions + 1
        else if(index(line,'solute name').ne.0) then      ! The atom label of the solute to calculate the RDF from
            read(line,*), dummy, solutemicro
            noptions = noptions + 1
        else if(index(line,'solvent charge sites').ne.0) then      ! The number of solvent charge sites
            read(line,*), dummy, dummy2, dummy3, chargemicro
            noptions = noptions + 1
        else if(index(line,'solute number').ne.0) then      ! The number of solute atoms
           print*,'2',line
            read(line,*), dummy, dummy2, atomnummicro
            noptions = noptions + 1
        else if(index(line,'END').ne.0) then
            print*, 'Input file read'
            write(500,'(A)') 'Input file read'
            if (noptions.lt.4) then
               write(500,'(A)') 'ERROR: Too few arguments passed for task microsolv'
               stop 'ERROR: Too few arguments passed for task microsolv'
            else if (noptions.gt.4) then
               write(500,'(A)') 'ERROR: Too many arguments passed for task microsolv'
               stop 'ERROR: Too many arguments passed for task micosolv'
            end if
            ierr=1
        end if
     end do
    connect = 'GEN'
    call RDFmicro(Omicro, solutemicro, chargemicro, atomnummicro, connect)
    write(500,'(A)') 'Micro-solvation ended normally'
!
 else if(index(task,'2').ne.0) then
     do while(ierr.eq.0)
        read(700,'(A)') line
        if(index(line,'Otype').ne.0) then            ! Water O label eg OTP
            read(line,*), dummy, Omicro
            noptions = noptions + 1
        else if(index(line,'solute name').ne.0) then      ! The atom label of the solute to calculate the RDF from
            read(line,*), dummy, solutemicro
            noptions = noptions + 1
        else if(index(line,'solvent charge sites').ne.0) then      ! The number of solvent charge sites
            read(line,*), dummy, dummy2, dummy3, chargemicro
            noptions = noptions + 1
         else if(index(line,'solute number').ne.0) then      ! The number of solute atoms 
            read(line,*), dummy, dummy2, atomnummicro
            noptions = noptions + 1
        else if(index(line,'END').ne.0) then
            print*, 'Input file read'
            write(500,'(A)') 'Input file read'
            if (noptions.lt.4) then
               write(500,'(A)') 'ERROR: Too few arguments passed for task 2'
               stop 'ERROR: Too few arguments passed for task 2'
            else if (noptions.gt.4) then
               write(500,'(A)') 'ERROR: Too many arguments passed for task 2'
               stop 'ERROR: Too many arguments passed for task 2'
            end if
            ierr=1
        end if
     end do
    connect = 'GEN'
    print*,'Solute atom count is ',atomnummicro
    call RDFmicro(Omicro, solutemicro, chargemicro, atomnummicro, connect)
    write(500,'(A)') 'Micro-solvation ended normally'
!
!- - - - - - - - - Task 3 - - - - - - - - - 
!    
 else if(index(task,'solvshell').ne.0) then
    do while(ierr.eq.0)
        read(700,'(A)') line
        if(index(line,'Otype').ne.0) then                  ! Water O label eg OTP
            read(line,*), dummy, Omicro
            noptions = noptions + 1
        else if(index(line,'solute filename').ne.0) then      ! filename of a file which has the atom labels of the solute to calculate the RDF from
            read(line,*), dummy, dummy2, solutefname
            noptions = noptions + 1
        else if(index(line,'solvent charge sites').ne.0) then      ! The number of solvent charge sites 
            read(line,*), dummy, dummy2, dummy3, chargemicro
            noptions = noptions + 1
        else if(index(line,'solute number').ne.0) then      ! The number of solute atoms 
            read(line,*), dummy, dummy2, atomnummicro
            noptions = noptions + 1
        else if(index(line,'END').ne.0) then
            print*, 'Input file read'
            write(500,'(A)') 'Input file read'
            if (noptions.lt.4) then
               write(500,'(A)') 'ERROR: Too few arguments passed for task solvshell'
               stop 'ERROR: Too few arguments passed for task solvshelll'
            else if (noptions.gt.4) then
               write(500,'(A)') 'ERROR: Too many arguments passed for task solvshell'
               stop 'ERROR: Too many arguments passed for task solvshell'
            end if
            ierr = 1
        end if
    end do
    call rdf(Omicro, solutefname, chargemicro, atomnummicro)
    write(500,'(A)') 'First solvation ended normally'
!
 else if(index(task,'3').ne.0) then
    do while(ierr.eq.0)
        read(700,'(A)') line
        if(index(line,'Otype').ne.0) then                  ! Water O label eg OTP
            read(line,*), dummy, Omicro
            noptions = noptions + 1
        else if(index(line,'solute filename').ne.0) then      ! filename of a file which has the atom labels of the solute to calculate the RDF from
            read(line,*), dummy, dummy2, solutefname
            noptions = noptions + 1
        else if(index(line,'solvent charge sites').ne.0) then      ! The number of solvent charge sites 
            read(line,*), dummy, dummy2, dummy3, chargemicro
            noptions = noptions + 1
        else if(index(line,'solute number').ne.0) then      ! The number of solute atoms 
            read(line,*), dummy, dummy2, atomnummicro
            noptions = noptions + 1
        else if(index(line,'END').ne.0) then
            print*, 'Input file read'
            write(500,'(A)') 'Input file read'
            if (noptions.lt.4) then
               write(500,'(A)') 'ERROR: Too few arguments passed for task 3'
               stop 'ERROR: Too few arguments passed for task 3'
            else if (noptions.gt.4) then
               write(500,'(A)') 'ERROR: Too many arguments passed for task 3'
               stop 'ERROR: Too many arguments passed for task 3'
            end if
            ierr = 1
        end if
    end do
    call rdf(Omicro, solutefname, chargemicro, atomnummicro)
    write(500,'(A)') 'First solvation ended normally'
!
!- - - - - - - - - Task 4 - - - - - - - - -     
!
 else if(index(task,'statskrig').ne.0) then
     do while(ierr.eq.0)
         read(700,'(A)') line
        if(index(line,'multipole').ne.0) then                  ! For multipole moments
            noptions = noptions + 1
            call krigcoefficent
        else if(index(line,'IQA').ne.0) then                   ! For IQA
            noptions = noptions + 1
            call krigcoefficentiqa 
        else if(index(line,'END').ne.0) then
            print*, 'Input file read'
            write(500,'(A)') 'Input file read'
            if (noptions.eq.0) then
               write(500,'(A)') 'ERROR: Too few arguments passed for task statskrig'
               stop 'ERROR: Too few arguments passed for task statskrig'
            end if
            ierr = 1
        end if
     end do
!
 else if(index(task,'4').ne.0) then
     do while(ierr.eq.0)
         read(700,'(A)') line
        if(index(line,'multipole').ne.0) then                  ! For multipole moments
            noptions = noptions + 1
            call krigcoefficent
        else if(index(line,'IQA').ne.0) then                   ! For IQA
            noptions = noptions + 1
            call krigcoefficentiqa
        else if(index(line,'END').ne.0) then
            print*, 'Input file read'
            write(500,'(A)') 'Input file read'
            if (noptions.eq.0) then
               write(500,'(A)') 'ERROR: Too few arguments passed for task 4'
               stop 'ERROR: Too few arguments passed for task 4'
            end if
            ierr = 1
        end if
     end do

!
!- - - - - - - - - Task 5 - - - - - - - - -     
!
 else if(index(task,'featorder').ne.0) then
    do while(ierr.eq.0)
        read(700,'(A)') line
        if(index(line,'Training_set_name').ne.0) then         
            read(line,*), dummy, TSname
            TSname=trim(adjustl(TSname))
            noptions = noptions + 1                  
        else if(index(line,'No_Solute_atoms').ne.0) then        
            read(line,*), dummy, AtomsBefore
            noptions = noptions + 1
        else if(index(line,'No_waters').ne.0) then        
            read(line,*), dummy, NoWater
            noptions = noptions + 1
        else if(index(line,'Resolution').ne.0) then        
            read(line,*), dummy, res
            noptions = noptions + 1
        else if(index(line,'Overlap').ne.0) then        
            read(line,*), dummy, O_lap
            noptions = noptions + 1
        else if(index(line,'Feat_var_thresh').ne.0) then        
            read(line,*), dummy, fvt
            noptions = noptions + 1
        else if(index(line,'No_bins').ne.0) then        
            read(line,*), dummy, Nbins
            noptions = noptions + 1
        else if(index(line,'Min_waters').ne.0) then
            read(line,*), dummy, MinWaters
            noptions = noptions + 1
        else if(index(line,'END').ne.0) then
            print*, 'Input file read'
            write(500,'(A)') 'Input file read'
            if (noptions.lt.8) then
               write(500,'(A)') 'ERROR: Too few arguments passed for task structure', noptions
               stop 'ERROR: Too few arguments passed for task structure'
            else if (noptions.gt.8) then
               write(500,'(A)') 'ERROR: Too many arguments passed for task structure', noptions
               stop 'ERROR: Too many arguments passed for task structure'
            end if
            ierr=1
        end if
    end do
    call density_nodes(TSname, AtomsBefore, NoWater, res, O_lap)
    call waterorder(TSname, AtomsBefore, NoWater)
    call anova(TSname, AtomsBefore, NoWater, fvt, Nbins, MinWaters)
    call node_order(natoms, AtomsBefore, NoWater)
    call solvent_order(TSname, AtomsBefore, NoWater)
 
    inquire (file='FINPUT-original.txt', exist=file_ex)
    inquire (file='FINPUT.txt', exist=file_ex2)
    istatus = 0 
    if (file_ex.eqv..FALSE..and.file_ex2.eqv..TRUE.) then
       open(status="new", unit=706, file='FINPUT-original.txt', action="write",iostat=ierr)
       if(ierr > 0) then
          print*, 'ERROR - opening new file FINPUT-original.txt'
       end if
       ierr=0
       open(status="old", unit=705, file='FINPUT.txt', action="read",iostat=ierr)
       if(ierr > 0) then
          print*, 'ERROR - opening FINPUT.txt'
       end if
       ierr=0
       
       do while (ierr.eq.0)
          read(705,'(A)',iostat=ierr) buffer
          write(706,'(A)') buffer
       end do
       if(ierr > 0) then
          print*, 'ERROR - Reading FINPUT.txt to copy to FINPUT-original.txt'
       end if
       close(705)
       close(706)
       
       ierr = 0
       open(status="replace", unit=706, file='FINPUT.txt', action="write",iostat=ierr)
       if(ierr > 0) then
          print*, 'ERROR - opening new file FINPUT.txt to overwrite'
       end if
       ierr = 0
       open(status="old", unit=705, file='FINPUT-original.txt', action="read",iostat=ierr)
       if(ierr > 0) then
          print*, 'ERROR - opening FINPUT-original.txt'
       end if
       ierr = 0

       do while (ierr.eq.0)
          read(705,'(A)',iostat=ierr) buffer
          if(index(buffer,'natoms').ne.0) then
             write(706,'(A6,1X,I5)') 'natoms',natoms
          else
             write(706,'(A)') buffer
          end if
       end do
    else
       print*, 'Cannot correct the FINPUT.txt file automatically this should be done by hand'
       print*, 'Change the values of natoms in the FINPUT.txt file to the value in '
       print*, 'Change_natoms_in_FINPUT_to_value_in_this_file.txt'
    end if
          write(500,'(A)') 'Feature ordering finished normally'
!
 else if (index(task,'5').ne.0) then
    do while(ierr.eq.0)
        read(700,'(A)') line
        if(index(line,'Training_set_name').ne.0) then         
            read(line,*), dummy, TSname
            TSname=trim(adjustl(TSname))
            noptions = noptions + 1                  
        else if(index(line,'No_Solute_atoms').ne.0) then        
            read(line,*), dummy, AtomsBefore
            noptions = noptions + 1
        else if(index(line,'No_waters').ne.0) then        
            read(line,*), dummy, NoWater
            noptions = noptions + 1
        else if(index(line,'Resolution').ne.0) then        
            read(line,*), dummy, res
            noptions = noptions + 1
        else if(index(line,'Overlap').ne.0) then        
            read(line,*), dummy, O_lap
            noptions = noptions + 1
        else if(index(line,'Feat_var_thresh').ne.0) then        
            read(line,*), dummy, fvt
            noptions = noptions + 1
        else if(index(line,'No_bins').ne.0) then        
            read(line,*), dummy, Nbins
            noptions = noptions + 1
        else if(index(line,'Min_waters').ne.0) then
            read(line,*), dummy, MinWaters
            noptions = noptions + 1
        else if(index(line,'END').ne.0) then
            print*, 'Input file read'
            write(500,'(A)') 'Input file read'
            if (noptions.lt.8) then
               write(500,'(A)') 'ERROR: Too few arguments passed for task structure', noptions
               stop 'ERROR: Too few arguments passed for task structure'
            else if (noptions.gt.8) then
               write(500,'(A)') 'ERROR: Too many arguments passed for task structure', noptions
               stop 'ERROR: Too many arguments passed for task structure'
            end if
            ierr=1
        end if
    end do
    call density_nodes(TSname, AtomsBefore, NoWater, res, O_lap)
    call waterorder(TSname, AtomsBefore, NoWater)
    call anova(TSname, AtomsBefore, NoWater, fvt, Nbins, MinWaters)
    call node_order(natoms, AtomsBefore, NoWater)
    call solvent_order(TSname, AtomsBefore, NoWater)
 
    inquire (file='FINPUT-original.txt', exist=file_ex)
    inquire (file='FINPUT.txt', exist=file_ex2)
    istatus = 0 
    if (file_ex.eqv..FALSE..and.file_ex2.eqv..TRUE.) then
       open(status="new", unit=706, file='FINPUT-original.txt', action="write",iostat=ierr)
       if(ierr > 0) then
          print*, 'ERROR - opening new file FINPUT-original.txt'
       end if
       ierr=0
       open(status="old", unit=705, file='FINPUT.txt', action="read",iostat=ierr)
       if(ierr > 0) then
          print*, 'ERROR - opening FINPUT.txt'
       end if
       ierr=0
       
       do while (ierr.eq.0)
          read(705,'(A)',iostat=ierr) buffer
          write(706,'(A)') buffer
       end do
       if(ierr > 0) then
          print*, 'ERROR - Reading FINPUT.txt to copy to FINPUT-original.txt'
       end if
       close(705)
       close(706)
       
       ierr = 0
       open(status="replace", unit=706, file='FINPUT.txt', action="write",iostat=ierr)
       if(ierr > 0) then
          print*, 'ERROR - opening new file FINPUT.txt to overwrite'
       end if
       ierr = 0
       open(status="old", unit=705, file='FINPUT-original.txt', action="read",iostat=ierr)
       if(ierr > 0) then
          print*, 'ERROR - opening FINPUT-original.txt'
       end if
       ierr = 0

       do while (ierr.eq.0)
          read(705,'(A)',iostat=ierr) buffer
          if(index(buffer,'natoms').ne.0) then
             write(706,'(A6,1X,I5)') 'natoms',natoms
          else
             write(706,'(A)') buffer
          end if
       end do
    else
       print*, 'Cannot correct the FINPUT.txt file automatically this should be done by hand'
       print*, 'Change the values of natoms in the FINPUT.txt file to the value in '
       print*, 'Change_natoms_in_FINPUT_to_value_in_this_file.txt'
    end if
          write(500,'(A)') 'Feature ordering finished normally'
!
!- - - - - - - - No input - - - - - - - - -
!
 else
    stop 'No recognised input task in input file'
 end if
! 
!----------------------------------------- Close files -------------------------------------------------------------------
!
 close(700)
 write(500,'(A)') 'Hermes ended normally from ', task, ', which has completed'
 close(500)
 print'(A)', 'Hermes ended normally after, ', task, ' was completed'

end Program Hermes
