!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Program author Robert Coates, James L. McDonagh and Stuart Davie Created 2015 in the            !
! Popelier group, Univeristy of Manchester                                                        !
! Licensed under Attribution-ShareAlike 2.5 Generic (CC BY-SA 2.5)                                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Originally created by Rob with oversight from S.D. Extensive edit and expansion by JMcD 
! Version 1.1
! CHANGE LOG 
! Version 1  : Used to make the micro solvated shell around N atom of amino acid only
! Version 1.1: Expanded to cut out the 1st solvation shell. Modular format and incorporation 
!              in the Hermes program. Still resitricted to H2O as the only solvent
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    module HISTORYrdf

     use module_subroutines
     use module_functions

     implicit none

    contains

    subroutine rdf(solvent_atom_type, solute_file, csolvent_charge_sites, cno_solute_atoms)
!
! - - - - - - - - - - - - - - - - - - - - - - Variables - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!	
    integer :: numatms, atms, ierror, stringlen, i, tstep, stepcount, total, counter, narg, solvent_charge_sites, no_solute_atoms
    integer :: charges, sites, no_solvent_molecules, tstep25, tstep50, tstep75, solute_count, iostat, linecount, num_atom_type
    integer :: soluatom, label_len, j, unique_count, l, m, n, k, icharlensolu, icharlensolv
    integer, dimension(:), allocatable :: timeindex, atom_type_count
    character (len=256) :: lineread, lrd, filename, newfilename
    character (len=100) :: coordfile, solute_file
    character (len=115) :: CMD
    character (len=4) :: element
    character (len=9) :: solvent_atom_type, csolvent_charge_sites, cno_solute_atoms
    character (len=8) :: labelled_element
    character (len=6), dimension(100) :: solutes_atom_type
    character (len=6), dimension(:,:), allocatable :: atom_types
    character (len=10), dimension(:), allocatable :: unique_atom_types
    character (len=17) :: dummyn
    real, dimension(:), allocatable :: nx, ny, nz, shells
    real, dimension(3,200) :: waterstore
    real, dimension(:,:), allocatable :: solutexyz
    real :: x, y, z, atomdist
    real :: bl, shell1, shell_size
    character (len=8) :: search = 'timestep'				! Search HISTORY file for steps
    logical :: fileexist, rdfexist, done, solute, waterinout					! Check the file exists
	
! - - - - - - - - - - - - - - - - - - - - - - Read inputs - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    print*, '------------------------------------------------------------------------------'
    print*, '|* RDF-CREATE * VERSION 1.2 * By R.Coates, J. L. McDonagh and S. Davie       |'
    print*, '|****************************************************************************|'
    print*, '|* THIS * PROGRAM * WILL * TRANSFORM * HISTORY * FILES * TO * .XYZ * FILES * |'
    print*, '|****************************************************************************|'
    print*, '|PLEASE ENTER; H2O O ATOM TYPE, SOLUTE ATOM TYPE AS THE RDF SOURCE, # SOLVENT|'
    print*, '|CHARGE SITES PER MOLECULE, NAME OF FILE CONTAINING ONE ATOM TYPE PER LINE TO|'
    print*, '|CALCULATE THE RDF FROM AND THE NAME OF AN XYZ FILE OR GEN TO GENERATE       |'
    print*, '|****************************************************************************|'
  
    filename = 'HISTORY'							! Reading HISTORY
    
    stringlen = len(trim(filename))					! Name .xyz file using HISTORY file name
    newfilename = filename(1:stringlen)
    newfilename(stringlen + 1:stringlen + 4) = '.xyz'
    stringlen = len(trim(newfilename))
    solvent_charge_sites = 3

    num_atom_type = 15                              ! Change to the total number of atom type definitions the program has
    allocate(atom_types(num_atom_type,2))
    allocate(atom_type_count(num_atom_type))
    atom_type_count = 0
    atom_types(1,1) = 'NH3'                         ! Atom type definintion add more for other systems. If you add one you must add the length of the atom type string
    atom_types(1,2) = '3'                           ! The length of the atom type string is very important these correspond to the atom types
    atom_types(2,1) = 'HC'                            
    atom_types(2,2) = '2'
    atom_types(3,1) = 'CT1'
    atom_types(3,2) = '3'
    atom_types(4,1) = 'HB'
    atom_types(4,2) = '2'
    atom_types(5,1) = 'CC'
    atom_types(5,2) = '2'
    atom_types(6,1) = 'HA'
    atom_types(6,2) = '2'
    atom_types(7,1) = 'CT3'
    atom_types(7,2) = '3'
    atom_types(8,1) = 'OC'
    atom_types(8,2) = '2'
    atom_types(9,1) = 'C '
    atom_types(9,2) = '2'
    atom_types(10,1) = 'O '
    atom_types(10,2) = '2'
    atom_types(11,1) = 'NH1'
    atom_types(11,2) = '3'
    atom_types(12,1) = 'H '
    atom_types(12,2) = '2'
    atom_types(13,1) = 'CT2'
    atom_types(13,2) = '3'
    atom_types(14,1) = 'CA'
    atom_types(14,2) = '2'
    atom_types(15,1) = 'HP'
    atom_types(15,2) = '2'
!    
    solvent_atom_type = trim(adjustl(solvent_atom_type))
    icharlensolv = len(trim(adjustl(solvent_atom_type)))
!
    open(status="old", unit=40, file=solute_file, action="read", iostat=ierror, access="sequential")	! Open solute atom type containing file for reading. This file should be a list one per line of the atom types of the solute atoms.
    solute_count = 0
    if (ierror.gt.0) then
       print('(A)'), '|***** ERROR - There is a problem opening the file containing the solutes atom type ',&
            trim(adjustl(solute_file))
       stop '|****************************************************************************|'
    end if

    do while (ierror.eq.0)
       solute_count = solute_count + 1
       
       if (solute_count.gt.100) then
          stop 'ERROR - Number of solute atoms too high, set solute_count in main-md-to-xyz.f90 high, maximum default is 100'
       end if
        
       read(40,'(A)', iostat=ierror) solutes_atom_type(solute_count)
       if(index('END',trim(adjustl(solutes_atom_type(solute_count)))).gt.0)then
          ierror = 1
          exit
       else
          solutes_atom_type(solute_count) = trim(adjustl(solutes_atom_type(solute_count)))
          print*, solutes_atom_type(solute_count)
       end if
!
       if(ierror.gt.0) then
          print('(A)'), 'ERROR - There is a problem reading the values in the file &
               ',trim(adjustl(solute_file)), ' . Please check this file.'
          stop
       else if(ierror.lt.0) then
          print('(A)'), 'WARNING - Solute atom types have been found but file did not end with END.&
                  Number of solute atoms to evaluate rdfs at - ',solute_count
       end if
    end do
    close(40)
!
! - - - - - - - - - - - - - - - - - - - - - - O P E N - F I L E S - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! 	
    open (status="old", unit=60, action="read", access="sequential", file=filename, iostat=ierror)	 ! Open the HISTORY, to read only
    if (ierror.ne.0) then
       print*, '|*********** ERROR * OPENING * HISTORY * FILE *******************************|'
       stop '|****************************************************************************|'
    else
       print*, '|*********** HISTORY * OPENED * SUCCESSFULLY ********************************|'
       print*, '|****************************************************************************|'
    end if
!  
    print*, '|*********** OPENING * NEW * FILES ******************************************|'
    print*, '|****************************************************************************|'
!    
    open(status="replace", unit=50, file="history.xyz", action="write", iostat=ierror)			! Open new history file for writing
    if (ierror.ne.0) then
      print*, '|*********** there was an error creating the new history.xyz file ***********|'
      stop '|****************************************************************************|'
    else
      print*, '|*********** the history.xyz file was sucessfully created *******************|'
               
    end if
!    
    open(status="replace", unit=51, file="historyxyz.rdf", action="readwrite", iostat=ierror)	        ! Open new historyxyz file for writing
    if (ierror.ne.0) then
      print*, '|*********** there was an error creating the new historyxyz.rdf file ********|'
      stop '|****************************************************************************|'
    else
      print*, '|*********** the historyxyz.rdf file was sucessfully created ****************|'
    end if
!	
    open(status="replace", unit=52, file="rdf", action="write", iostat=ierror) ! Open new rdf file for writing
    if (ierror.ne.0) then
        print*, '|*********** there was an error creating the new rdf file *******************|'
        stop '|****************************************************************************|'
    else
        print*, '|*********** the rdf file was sucessfully created ***************************|'
    end if 
!   
    open(status="replace", unit=53, file="temp-his-cut.xyz", action="readwrite", iostat=ierror)			
    if (ierror.ne.0) then
        print*, '|*********** there was an error creating the new temp-his-cut.xyz file ******|'
        stop '|****************************************************************************|'
    else
        print*, '|*********** the temp-his-cut.xyz file was sucessfully created **************|'
    end if
!
    open(status="replace", unit=54, file="HISTORYCUT.xyz", action="write", iostat=ierror)			
    if (ierror.ne.0) then
        print*, '|*********** there was an error creating the new HISORYCUT.xyz file *********|'
        stop '|****************************************************************************|'
    else
        print*, '|*********** the HISTORYCUT.xyz file was sucessfully created ****************|'
     end if
!
     open(status="replace", unit=55, file="Hermes-rdf.log", action="write", iostat=ierror)			
    if (ierror.ne.0) then
        print*, '|*********** there was an error creating the new Hermes-rdf.log file ********|'
        stop '|****************************************************************************|'
    else
        print*, '|*********** the Hermes-rdf.log file was sucessfully created ***************|'
    end if
!	
    print*, '|*********** FILES * OPENED *************************************************|'
    print*, '|****************************************************************************|'
!	
! - - - - - - - - - - - - - - - - - - - - - - I N I T I A L I S A T I O N - - - - - - - - - - - - - - - - - - - -
!	
    counter = 0
    shell1 = 0
    bl = 0
    linecount = 0 
    i = 0
    j = 0
    k = 0
    l= 0
    m = 0
    n = 0
    tstep = 0
    stepcount = 0
    numatms = 0
    atms = 0
    total = 0
    ierror = 0
    atom_type_count = 0
    unique_count = 0
!	
! - - - - - - - - - - - - - - - - - - - - - - C O U N T I N G - T H E - T I M E S T E P S / A T O M S - - - - - - - - - - - - - -
!
    print*, '|*********** COUNTING * ATOMS/FRAMES ****************************************|'
    print*, '|****************************************************************************|'
!    
    do while (ierror.eq.0)
       read(60, '(A)', iostat=ierror) lineread		! Read HISTORY line by line
       read(60, '(A)', iostat=ierror) lineread
       lrd = adjustl(lineread(23:35))
       read(lrd, '(I6)') numatms
       lrd = adjustl(lineread(46:53))
       read(lrd, '(I5)') tstep
       read(60, '(A)', iostat=ierror) lineread		
       read(60, '(A)', iostat=ierror) lineread
       lrd = adjustl(lineread(6:22))
       read(lrd, "(ES17.10)") bl
       exit
    end do
!	
! - - - - - - - - - - - - - - - - - - - - - - P R I N T - A T O M S / S N A P S - - - - - - - - - - - - - - - - - - - - - - - - -
!
    print*, tstep, 'Frames found.'
    print*, numatms, 'Atoms found, including water charges.'
    if (solvent_charge_sites.le.3) then
       atms = numatms              !total number of atoms taking away charges on water ! JMcD TIP5P numatms - int((numatms-13)*(0.4))
       print*, atms, 'Atoms found. No additional water charges found.'
    else
       no_solvent_molecules = (numatms - no_solute_atoms) ! JMcD The number of solvent molecules 
       atms = (no_solvent_molecules * 3) + no_solute_atoms!total number of atoms taking away charges on water
       print*, atms, 'Atoms found not including additional water charges.'
    end if
!    
    allocate(unique_atom_types(int(atms/3)))
!   
    print*, numatms-no_solute_atoms
    print*, solvent_charge_sites
    total = (numatms-no_solute_atoms)/solvent_charge_sites!total number of water molecules! JMcD TIP5P int(((numatms-13)*(0.2)))
    total = int(total)
    print*, total, 'Solvent water molecules found.'
! 
! - - - - - - - - - - - - - - - - - - - - - - R E W I N D - A N D - R E I N I T I A L I S E - - - - - - - - - - - - - - - - - - -
!
    linecount = 0
    rewind(60)
    ierror = 0
!    
! - - - - - - - - - - - - - - - - Work through atom types and solvent atoms to determine which remain - - - - - - - - - - - - - -	
!
    tstep25 = int(tstep*0.25)
    tstep50 = int(tstep*0.5)
    tstep75 = int(tstep*0.75)
    print*, '|*********** Creating xyz files *********************************************|'

    do while (ierror.eq.0)
          read(60, '(A)', iostat=ierror) lineread
          if (index(lineread, 'timestep').gt.0) then	! Starts at the timestep on the HISTORY
             atom_type_count = 0 
             unique_count = 0
             stepcount = stepcount + 1
             if (stepcount.eq.tstep25) then
                print*, '|*********** 25% ************************************************************|'
             else if (stepcount.eq.tstep50) then
                print*, '|*********** 50% ************************************************************|'
             else if (stepcount.eq.tstep75) then
                print*, '|*********** 75% ************************************************************|'
             else if (stepcount.eq.tstep) then
                print*, '|*********** 100% ***********************************************************|'
             end if
             write(50, "(I20)") atms
             write(50, "(A, I7)") 'Timestep: ', stepcount
             write(51, "(I20)") atms
             write(51, "(A, I7)") 'Timestep: ', stepcount
!            
             do i = 1,3
                 read(60, '(A)', iostat=ierror) lineread
             end do
!            
             do i = 1, numatms
                solute = .FALSE.
                read(60, '(A)', iostat=ierror) lineread
                element = trim(adjustl(lineread(1:4)))
                read(60, '(A)', iostat=ierror) lineread
                read(lineread(4:21), "(ES17.10)") x	
                read(lineread(24:41), "(ES17.10)") y
                read(lineread(44:61), "(ES17.10)") z
!               
                if (index(element, 'Q').gt.0) then	! Skips the tp5p charges
                   cycle
                end if
!

                do j = 1, num_atom_type, 1 ! JMcD this loop uniquely labels atoms by counting the number of times an atom type is present in a Timestep 
                   if (index(element(1:1),atom_types(j,1)(1:1)).eq.1.and.index(element(2:2),atom_types(j,1)(2:2)).eq.1) then ! JMcD Only matches the first two characters of the atom type unique labels are given even if two are the same.
                      atom_type_count(j) = atom_type_count(j) + 1
                      unique_count = unique_count + 1
                      write(labelled_element,'(A4,I4)') element, atom_type_count(j)
                      call StripSpaces(labelled_element)
                      write(50, "(A2, 4X, ES17.10, 2X, ES17.10, 2X, ES17.10)") element(1:1), x, y, z
                      write(51, "(A, 2X, ES17.10, 2X, ES17.10, 2X, ES17.10)") trim(labelled_element), x, y, z
                      unique_atom_types(unique_count) = labelled_element
                      solute = .TRUE.
                      exit
                   end if
                end do
!               
                if (solute.eqv..FALSE.) then
                    write(50, "(A2, 4X, ES17.10, 2X, ES17.10, 2X, ES17.10)") element(1:1), x, y, z
                    write(51, "(A, 2X, ES17.10, 2X, ES17.10, 2X, ES17.10)") element, x, y, z
                end if
!                
             end do
          end if
    end do
    close(50)
    rewind(51)
    allocate(shells(unique_count))
    allocate(solutexyz(3,unique_count))
    deallocate(atom_type_count)
    print*, '|*********** history.xyz made ***********************************************|'
!    
    allocate(nx(tstep), ny(tstep), nz(tstep)) !an array containing the coordinates of the target atom at each frame
!    
    do i = 1, unique_count
       rewind(51)
       call rdfcreate(tstep, bl, total, shell1, nx, ny, nz, solvent_atom_type, unique_atom_types(i))
       shells(i) = shell1
       shell_size = maxval(shells)
       print*, '|*********** ', unique_atom_types(i), shell1,' ***********************************|'
    end do
    print*, shell_size
    rewind(51)
    ierror = 0
!
! - - - - - - - - - - - - - - - - - - - - - - C U T - C O O R D I N A T E S - - - - - - - - - - - - - - - - - - - - - - - - - -
! WRITES COORDINATES TO A TEMORARY FILE, THE NUMBER OF ATOMS IS COUNTED LATER BY READING THIS FILE
    print*, '|*********** CUTTING * OUT * 2ND * SOLV * SHELL * UPWARDS *******************|'
    write(55,*) 'Solute,          x,           y,            z,      Solvent,       x       y,        z,    Separating Distance'
!   '(A6,X13,A1,2(X17,A1),X5,A7,X8,A1,2(X17,A1),X5,A20)')
    do while (ierror.eq.0)
        read(51, '(A)', iostat=ierror) lineread               ! 51 is historyxyz.rdf
        if (index(lineread, 'Timestep').gt.0) then	      ! Print timestep between each set of coordinates
            i = i + 1
            if (i.eq.tstep25) then
                print*, '|*********** 25% ************************************************************|'
            else if (i.eq.tstep50) then
                print*, '|*********** 50% ************************************************************|'
            else if (i.eq.tstep75) then
                print*, '|*********** 75% ************************************************************|'
            else if (i.eq.tstep) then
                print*, '|*********** 100% ***********************************************************|'
            end if
            write(53, '(A)') lineread                         ! 53 is temp-his-cut.xyz
        end if
                
        do n = 1, unique_count, 1
            icharlensolu = len(trim(adjustl(unique_atom_types(n))))
            if (index(lineread, unique_atom_types(n)(1:icharlensolu)).gt.0)then	 ! alanine elements, labels dependent on model
                read(lineread(1:2), "(A1)") element
                read(lineread,*), dummyn, x, y, z
                write(53, "(A1, 4X, ES17.10, 2X, ES17.10, 2X, ES17.10)") element, x, y, z
                solutexyz(1,n) = x 
                solutexyz(2,n) = y
                solutexyz(3,n) = z
            else if (index(lineread, solvent_atom_type(1:icharlensolv)).gt.0) then      		 ! tp5p O label, could vary OT5E TIP5P OTP TIP3P
                read(lineread(1:2), "(A1)") element
                read(lineread,*), dummyn, x, y, z
                do while(waterinout.eqv..false.)
                    do l = 1, unique_count
                        atomdist = dist(solutexyz(1,l), solutexyz(2,l), solutexyz(3,l), x, y, z)
                        if (atomdist.le.shell_size) then 			! O is within first solvation shell
                            write(55,*) unique_atom_types(l), solutexyz(1,l), solutexyz(2,l), solutexyz(3,l), &
                                 element, x, y, z, atomdist
                            write(53, "(A1, 4X, ES17.10, 2X, ES17.10, 2X, ES17.10)") element, x, y, z	! The water protons are on the
                            read(51, '(A)', iostat=ierror) lineread					! next 2 lines of the file
                            read(lineread(1:2), "(A1)") element						! produced by dl poly ! therefore can be read and written directly
                            read(lineread,*), dummyn, x, y, z
                            write(53, "(A1, 4X, ES17.10, 2X, ES17.10, 2X, ES17.10)") element, x, y, z
                            read(51, '(A)', iostat=ierror) lineread
                            read(lineread(1:2), "(A1)") element
                            read(lineread,*), dummyn, x, y, z
                            write(53, "(A1, 4X, ES17.10, 2X, ES17.10, 2X, ES17.10)") element, x, y, z
                            waterinout = .true.
                            exit
                        else if(l.eq.unique_count) then
                            waterinout = .true.
                        end if 
                    end do
                end do
                waterinout = .false.
            end if
        end do
    end do
    print*, '|*********** DONE ***********************************************************|'
    print*, '|****************************************************************************|'
    close(51)					! Delete the temp file
    deallocate(nx, ny, nz)
    deallocate(solutexyz, shells)
!    
! - - - - - - - - - - - - - - - - - - - - - - C O U N T - T H E - N O - A T O M S - A T - E A C H - S T E P - - - - - - - - - -
!
    i = 0					        ! counts lines in the HISCUTtemp file
    counter = 0						! and stores the line number at
    ierror = 0						! 'timestep' lines' in the timeindex
    rewind(53)						! array for use calculating No. atoms
    allocate(timeindex(tstep+1))
    do while (ierror.eq.0)
        read(53, '(A)', iostat=ierror) lineread	
        counter = counter + 1
        if (index(lineread, 'Timestep').GT.0) then
            i = i + 1
            timeindex(i) = counter
        end if
    end do
    timeindex(tstep+1) = counter
!
! - - - - - - - - - - - - - - - - - - - - - - W R I T E - F I N A L - C U T - X Y Z - F I L E - - - - - - - - - - - - - - - - -		
!
     print*, '|*********** WRITING * HISCUT.XYZ *******************************************|'
     print*, '|****************************************************************************|'
     i = 0
     ierror = 0	
     rewind(53)
     do while (ierror.EQ.0)
         read(53, '(A)', iostat=ierror) lineread
         if (index(lineread, 'Timestep').GT.0) then
             i = i + 1
             if (i.EQ.tstep25) then
                 print*, '|*********** 25% ************************************************************|'
             else if (i.EQ.tstep50) then
                 print*, '|*********** 50% ************************************************************|'
             else if (i.EQ.tstep75) then
                 print*, '|*********** 75% ************************************************************|'
             else if (i.EQ.tstep) then
                 print*, '|*********** 100% ***********************************************************|'
             end if
             write(54, '(I20)') (timeindex(i+1) - timeindex(i) - 1) !No. atoms, the difference between two line numbers - 1
             write(54, '(A)') lineread
         else if (ierror.EQ.0) then
             write(54, '(A)') lineread	!write coordinates
         end if
     end do
     print*, '|*********** PROGRAM * COMPLETE *********************************************|'
     print*, '------------------------------------------------------------------------------'
!
! - - - - - - - - - - - - - - - - - - - - - - C L E A N - U P - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
     close(53)
     close(54)
     deallocate(atom_types)
     deallocate(unique_atom_types)
     deallocate(timeindex)
!
!- - - - - - - - - - - - - - - - - - - - - - close program - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
   end subroutine  RDF

 end module HISTORYrdf
