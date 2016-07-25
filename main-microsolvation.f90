!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Program author Robert Coates, Created 2015 in the Popelier group, Univeristy of Manchester      !
! Adapted and edited by James L. McDonagh and Stuart Davie, Univeristy of Manchester 2015         !
! Licensed under Attribution-ShareAlike 2.5 Generic (CC BY-SA 2.5)                                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This is the version originally created by R.C with oversight from S.D. and minor edits from JMcD 
! Version 1.1
! CHANGE LOG 
! Version 1  : Used to make the micro solvated shell around N atom of amino acid only
! Version 1.1: Expanded to run over any atom type the user enters. Modular format and incorporation 
!              in the Hermes program. Still resitricted to H2O as the only solvent
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  module HISTORYRDFmicro
  
    use module_subroutines
    use module_functions

    implicit none
    
    contains

      subroutine RDFmicro(solvent_atom_type, solutes_atom_type, csolvent_charge_sites, cno_solute_atoms, coordfile)
	
! - - - - - - - - - - - - - - - - - - - - - - V A R I A B L E S - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	
        integer :: numatms, atms, ierror, stringlen, i, tstep, stepcount, total, counter, narg, solvent_charge_sites 
        integer :: charges, sites, no_solvent_molecules, tstep25, tstep50, tstep75, E, S, Q, linecount, endcount, k 
        integer :: no_solute_atoms
        integer, dimension(:), allocatable :: timeindex
        character (len=256) :: lineread, lrd, filename, newfilename
        character (len=100) :: coordfile
        character (len=115) :: CMD
        character (len=4) :: element
        character (len=9) :: solvent_atom_type, solutes_atom_type, csolvent_charge_sites, cno_solute_atoms
        real, dimension(:), allocatable :: nx, ny, nz
        real, dimension(3) :: vector
        real :: x, y, z, dist, atomdist
        real :: bl, shell1, distance2
        character (len=8) :: search = 'timestep'				! Search HISTORY file for steps
        logical :: fileexist, rdfexist 					! Check the file exists
	
! - - - - - - - - - - - - - - - - - - - - - - R E A D - I N P U T S - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        print*, '------------------------------------------------------------------------------'
        print*, '|* MICROSOLVATION * VERSION 1.1 * By R. Coates, S. Davie and J. L. McDonagh  |'
        print*, '|****************************************************************************|'
        print*, '|* THIS * PROGRAM * WILL * TRANSFORM * HISTORY * FILES * TO * .XYZ * FILES * |'
        print*, '|****************************************************************************|'
        print*, '|PLEASE ENTER; H2O O ATOM TYPE, SOLUTE ATOM TYPE AS THE RDF SOURCE, # SOLVENT|'
        print*, '|CHARGE SITES PER MOLECULE, # OF SOLUTE ATOMS AND XYZ FILE OR GEN TO GENERATE|'
        print*, '|****************************************************************************|'
        
        filename = 'HISTORY'							! Reading HISTORY
        
        stringlen = len(trim(filename))					! Name .xyz file using HISTORY file name
        newfilename = filename(1:stringlen)
        newfilename(stringlen + 1:stringlen + 4) = '.xyz'
        stringlen = len(trim(newfilename))
!
        read(csolvent_charge_sites,'(I6)'), solvent_charge_sites
        print*,cno_solute_atoms
!
        read(cno_solute_atoms,'(I6)'), no_solute_atoms
! - - - - - - - - - - - - - - - - - - - - - - O P E N - F I L E S - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 	
        open (status="old", unit=60, action="read", access="sequential", file=filename, iostat=ierror)	!Open the HISTORY, to read only
        if (ierror.NE.0) then
           print*, '|*********** ERROR * OPENING * HISTORY * FILE *******************************|'
           print*, '|****************************************************************************|'
        else
           print*, '|*********** HISTORY * OPENED * SUCCESSFULLY ********************************|'
           print*, '|****************************************************************************|'
        end if
  
        print*, '|*********** OPENING * NEW * FILES ******************************************|'
        print*, '|****************************************************************************|'
  
        inquire (file=newfilename, exist=fileexist)			!Check if the files already exist (prevents overwriting data)
        inquire (file="rdf", exist=rdfexist)
        if(fileexist.EQV..TRUE..OR.rdfexist.EQV..TRUE.) then
           print*, '|*********** RDF * AND/OR * HISTORY.XYZ * FILES * ALREADY * EXISTS **********|'
           print*, '|****************************************************************************|'
           stop
        else
           if(coordfile.eq.'GEN') then
              open(status="new", unit=50, file=newfilename, action="write", iostat=ierror)		!Open new xyz file for writing
              if (ierror.NE.0) then
                 print*, 'there was an error creating the new xyz file'
              else
                 print*, 'the xyz file was successfully created'
              end if
           else
              open(status="old", unit=50, file=coordfile, action="read", iostat=ierror)		!Open new xyz file for writing
              if (ierror.NE.0) then
                 print*, 'there was an error creating the new xyz file'
              else
                 print*, 'the xyz file was successfully created'
              end if

              open(status="old", unit=51, file="xyztemp", action="readwrite", iostat=ierror)		!Open temporary file for writing
              if (ierror.NE.0) then
                 print*, 'there was an error creating the new temporary file'
              else
                 print*, 'the temorary file was successfully created'
              end if
        
              do while (ierror.eq.0)
                 read(50,'(A)',iostat=ierror) lineread
                 write(51,'(A)') lineread
              end do

           end if
		
           open(status="new", unit=51, file="xyztemp", action="readwrite", iostat=ierror)		!Open temporary file for writing
           if (ierror.NE.0) then
              print*, 'there was an error creating the new temporary file'
           else
              print*, 'the temorary file was successfully created'
           end if
    
           open(status="new", unit=52, file="rdf", action="write", iostat=ierror)			!Open new rdf file for writing
           if (ierror.NE.0) then
              print*, 'there was an error creating the new rdf file'
           else
              print*, 'the rdf file was sucessfully created'
           end if
    
           open(status="new", unit=53, file="HISCUTtemp", action="readwrite", iostat=ierror)		
           if (ierror.NE.0) then
              print*, 'there was an error creating the new HISCUTtemp file'
           else
              print*, 'the temporary file was sucessfully created'
           end if
           open(status="new", unit=54, file="HISCUT.xyz", action="write", iostat=ierror)			
           if (ierror.NE.0) then
              print*, 'there was an error creating the new HISCUT.xyz file'
           else
              print*, 'the HISCUT.xyz file was sucessfully created'
           end if
        end if
	
        print*, '|*********** FILES * OPENED *************************************************|'
        print*, '|****************************************************************************|'
	
! - - - - - - - - - - - - - - - - - - - - - - I N I T I A L I S A T I O N - - - - - - - - - - - - - - - - - - - -
	
        counter = 0
        shell1 = 0
        bl = 0
        linecount = 0 
        endcount = 0
        i = 0
        S = 0
        E = 0
        Q = 0
        tstep = 0
        stepcount = 0
 
! - - - - - - - - - - - - - - - - - - - - - - C O U N T I N G - T H E - T I M E S T E P S / A T O M S - - - - - - - - - - - - - -

        print*, '|*********** COUNTING * ATOMS/FRAMES ****************************************|'
        print*, '|****************************************************************************|'
        
        do while (ierror.EQ.0)
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

	
! - - - - - - - - - - - - - - - - - - - - - - P R I N T - A T O M S / S N A P S - - - - - - - - - - - - - - - - - - - - - - - - -
 
        print*, tstep, 'Frames found.'
        print*, numatms, 'Atoms found, including water charges.'
        if (solvent_charge_sites.le.3) then
           atms = numatms              !total number of atoms taking away charges on water ! JMcD TIP5P numatms - int((numatms-13)*(0.4))
        else
           no_solvent_molecules = (numatms - no_solute_atoms)  ! JMcD The number of solvent molecules 
           atms = (no_solvent_molecules * 3) + no_solute_atoms !total number of atoms taking away charges on water
        end if
        
        total = (numatms-no_solute_atoms)/solvent_charge_sites!total number of water molecules! JMcD TIP5P int(((numatms-13)*(0.2)))
        total = int(total)
 
! - - - - - - - - - - - - - - - - - - - - - - R E W I N D - A N D - R E I N I T I A L I S E - - - - - - - - - - - - - - - - - - -

        linecount = 0
        rewind(60)
        ierror = 0
	
! - - - - - - - - - - - - - - - - - - - - - - P R O D U C I N G - X Y Z - F I L E S - - - - - - - - - - - - - - - - - - - - - - -

        tstep25 = int(tstep*0.25)
        tstep50 = int(tstep*0.5)
        tstep75 = int(tstep*0.75)
 
        if (coordfile.eq.'GEN')then
           print*, '|*********** WRITING * XYZ * FILE *******************************************|'
           
           do while (ierror.EQ.0)
              read(60, '(A)', iostat=ierror) lineread
    
              If (index(lineread, 'timestep').GT.0) then	! Starts at the timestep on the HISTORY
                 stepcount = stepcount + 1
                 if (stepcount.EQ.tstep25) then
                    print*, '|*********** 25% ************************************************************|'
                 else if (stepcount.EQ.tstep50) then
                    print*, '|*********** 50% ************************************************************|'
                 else if (stepcount.EQ.tstep75) then
                    print*, '|*********** 75% ************************************************************|'
                 else if (stepcount.EQ.tstep) then
                    print*, '|*********** 100% ***********************************************************|'
                 end if
                 write(50, "(I20)") atms
                 write(50, "(A, I7)") 'Timestep: ', stepcount	! writes the coordinates to an xyz and a temp file for use calculating the rdf
                 write(51, "(A, I7)") 'Timestep: ', stepcount
                 read(60, '(A)', iostat=ierror) lineread
                 read(60, '(A)', iostat=ierror) lineread
                 read(60, '(A)', iostat=ierror) lineread
       
                 do i = 1, numatms
                    read(60, '(A)', iostat=ierror) lineread
                    element = adjustl(lineread(1:4))
                    read(60, '(A)', iostat=ierror) lineread
                    read(lineread(4:21), "(ES17.10)") x	
                    read(lineread(24:41), "(ES17.10)") y
                    read(lineread(44:61), "(ES17.10)") z
          
                    if (index(element, 'Q1').GT.0.OR.index(element, 'Q2').GT.0.) then	!Skips the tp5p charges
                       cycle
                    else
                       write(50, "(A2, 4X, ES17.10, 2X, ES17.10, 2X, ES17.10)") element(1:1), x, y, z
                       write(51, "(A4, 2X, ES17.10, 2X, ES17.10, 2X, ES17.10)") element, x, y, z
                    end if
                 end do
       
              end if
           end do
           print*, '|*********** XYZ * FILE * PRODUCED ******************************************|'
           print*, '|****************************************************************************|'
        else
           print*,'The following file will be read for the snapshot geometries, ', trim(coordfile)
        end if
        
        close(60)
        close(50)

        rewind(51) !rewind for reading by the rdfcreate subroutine
	
        allocate(nx(tstep), ny(tstep), nz(tstep)) !an array containing the coordinates of the target atom at each frame
        call rdfcreatemicro(tstep, bl, total, shell1, nx, ny, nz, solvent_atom_type, solutes_atom_type)	

        !gives subroutine the number of frames(tstep), the boxlength, total water molecules
	!subroutine returns 1st solv shell distance(shell1 and array of target atom coordinates

        close(52)	!Close RDF file
	
        rewind(51)
        ierror = 0
        i = 0
	
! - - - - - - - - - - - - - - - - - - - - - - C U T - C O O R D I N A T E S - - - - - - - - - - - - - - - - - - - - - - - - - -
	
! WRITES COORDINATES TO A TEMORARY FILE, THE NUMBER OF ATOMS IS COUNTED LATER BY READING THIS FILE
        print*, '|*********** CUTTING * OUT * 2ND * SOLV * SHELL * UPWARDS *******************|'
        do while (ierror.EQ.0)
           read(51, '(A)', iostat=ierror) lineread
           if (index(lineread, 'Timestep').GT.0) then		! Print timestep between each set of coordinates
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
              write(53, '(A)') lineread
           else if (index(lineread, 'NH3').GT.0.OR.index(lineread, 'HC').GT.0.OR.index(lineread, 'CT1').GT.0&
                .OR.index(lineread, 'HB').GT.0.OR.index(lineread, 'CC').GT.0.OR.index(lineread, 'HA').GT.0&
                .OR.index(lineread, 'CT3').GT.0.OR.index(lineread, 'OC').GT.0.OR.index(lineread, 'C ').GT.0&
                .OR.index(lineread, 'O ').GT.0.OR.index(lineread, 'NH1').GT.0.&
                .OR.index(lineread, 'H ').GT.0) then	 ! solute elements, labels dependent on model
              read(lineread(1:2), "(A1)") element
              read(lineread(7:24), "(ES17.10)") x	
              read(lineread(26:43), "(ES17.10)") y
              read(lineread(45:62), "(ES17.10)") z
              write(53, "(A1, 4X, ES17.10, 2X, ES17.10, 2X, ES17.10)") element, x, y, z
           else if (index(lineread, solvent_atom_type).GT.0) then      		 ! tp5p O label, could vary OT5E TIP5P OTP TIP3P
              read(lineread(1:2), "(A1)") element
              read(lineread(7:24), "(ES17.10)") x	
              read(lineread(26:43), "(ES17.10)") y
              read(lineread(45:62), "(ES17.10)") z
              vector(1) = nx(i) - x
              vector(2) = ny(i) - y
              vector(3) = nz(i) - z
              do k = 1,3,1
                distance2 = distance2 + (vector(k)**2.0)
              end do
              atomdist = sqrt(distance2)
              if (atomdist.LE.(shell1)) then 			! O is within first solvation shell
                 write(53, "(A1, 4X, ES17.10, 2X, ES17.10, 2X, ES17.10)") element, x, y, z	! The water protons are on the
                 read(51, '(A)', iostat=ierror) lineread					! next 2 lines of the file
                 read(lineread(1:2), "(A1)") element						! produced by dl poly
                 read(lineread(7:24), "(ES17.10)") x						! therefore can be read and written directly
                 read(lineread(26:43), "(ES17.10)") y
                 read(lineread(45:62), "(ES17.10)") z
                 write(53, "(A1, 4X, ES17.10, 2X, ES17.10, 2X, ES17.10)") element, x, y, z
                 read(51, '(A)', iostat=ierror) lineread
                 read(lineread(1:2), "(A1)") element
                 read(lineread(7:24), "(ES17.10)") x	
                 read(lineread(26:43), "(ES17.10)") y
                 read(lineread(45:62), "(ES17.10)") z
                 write(53, "(A1, 4X, ES17.10, 2X, ES17.10, 2X, ES17.10)") element, x, y, z
              end if
           end if
        end do
        print*, '|*********** DONE ***********************************************************|'
        print*, '|****************************************************************************|'
        close(51)					! Delete the temp file
        deallocate(nx, ny, nz)

! - - - - - - - - - - - - - - - - - - - - - - C O U N T - T H E - N O - A T O M S - A T - E A C H - S T E P - - - - - - - - - -
	
        i = 0							! counts lines in the HISCUTtemp file
        counter = 0						! and stores the line number at
        ierror = 0						! 'timestep' lines' in the timeindex
        rewind(53)						! array for use calculating No. atoms
        allocate(timeindex(tstep+1))
        do while (ierror.EQ.0)
           read(53, '(A)', iostat=ierror) lineread		
           counter = counter + 1
           if (index(lineread, 'Timestep').GT.0) then 
              i = i + 1
              timeindex(i) = counter
           end if
        end do
 
        timeindex(tstep+1) = counter
			
! - - - - - - - - - - - - - - - - - - - - - - W R I T E - F I N A L - C U T - X Y Z - F I L E - - - - - - - - - - - - - - - - -		
	
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
! - - - - - - - - - - - - - - - - - - - - - - C L E A N - U P - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        close(53)
        close(54)
	
! - - - - - - - - - - - - - - - - - - - - - - C L O S E - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      end subroutine   RDFmicro

    end module HISTORYRDFmicro


