!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Originally created by James McDonagh at the University of Manchester 2015, in the Popelier group      !
! Contributions from Dale Stuchfield                                                                    !
! Licensed under Attribution-ShareAlike 2.5 Generic (CC BY-SA 2.5)                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Originally created by James with contributions from Dale stuchfield.
! Version 1.1
! CHANGE LOG 
! Version 1  : Calculates the molecular deviations over a DL POLY trajectory or Tyche distortions. 
!              Statistics about the sampling are given.
! Version 1.1: Modular format and incorporation in the Hermes program.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module structurechanges

  use module_subroutines
  use module_bonding

  implicit none

 contains 
  subroutine stucture(filename, write_norm_xyz)

  integer :: ierr, status, line_no, geometry, atom_count, start_atc, end_atc, string_len, struc, stat, alloc_err
  integer :: column_no, ichar, location, location1, location2, location3, location4, status1, icharxyz, i, j, k, li, ichar_title
  integer :: icharAL1, icharAL2, icharAL3, icharAL4, unit70line_no, status2, status3, narg 
  real :: blength, bangle, dangle, bangledeg, dangledeg, possensedihed, bsum, bond_av, bond_sq_sum           ! in out terms for the subroutines
  real :: Bond_RMSD, asum, angle_av, angle_av_sum, angle_sq_sum, Angle_RMSD, dsum, dihed_av, dihed_av_sum, equil
  real :: dihed_sq_sum, Dihed_RMSD, bmin, bmax, amin, amax, Dmin, Dmax, bbias, bsd, abias, asd, Dbias, Dsd, alt
  real, dimension(3) :: atom_coord_1, atom_coord_2, atom_coord_3, atom_coord_4                               ! to read the atom coordinates into
  real, allocatable :: brmsd(:), armsd(:), drmsd(:), sorted(:)
  character (len=256) :: longstring, longstringxyz, title                                                    ! line of the read in file
  character (len=100) :: filename, new_filename, write_norm_xyz, fil_nam                                     ! filenames
  character (len=10) :: atco_search = "Atom Count"                                                           ! regular expression to find the atom count line
  character (len=1) :: calc, bond = "B", angle = "A", dihed = "D"                                            ! letter to define caluclation bondlength,bondangle or dihedral
  character (len=5) :: atomlabel1,atomlabel2, atomlabel3, atomlabel4, dummy                                  ! atom labels to find
  character (len=11) :: search = "-----------"                                                               ! string to find for each new geometry
  character (len=20) :: atom_xyz_1, atom_xyz_2, atom_xyz_3, atom_xyz_4 
  character(len=255) :: cwd
  character(len=20) :: names
  logical :: file_ex, FLB1,FLB2,FLB3,FLB4
!
! ---------------- Get User Input ---------------------------------------------------------------------------------------------------
! THE SECOND PART OF THIS LOOP SETS THE OUTPUT FILE NAME AS <FILENAME>.CSV
  print*, 'STRUCTURE VERSION 1'
  print*, 'This code will calculate all bond angle and dihedral deviation over a MD Trajectory'
  print*, 'The input is a list of all of the bonds angles and dihedrals you want this of'
  print*, 'The code can find all unique bonds, angles and dihedrals (note including redundent ones)'
  print*, 'For the code to do this enter auto to the filename' 
  print*, 'WARNING ALL PREVIOUS SUMMARY FILES WILL BE OVER WRITTEN'

  call system('mkdir -p BONDS_plot' )
  call system('mkdir -p ANGLES_plot' )
  call system('mkdir -p DIHEDRAL_plot' )
  call getcwd(cwd)
  cwd = trim(adjustl(cwd))
  write_norm_xyz = 'yes'

  open (status="replace", unit=200, file='structure_input.log', action="write",iostat=ierr)                     ! Unit 200 the structure log file
     if ( ierr.ne.0 ) then
        write(*,'(A)')'Module bonds : ERROR - unable to create structure_input.log please check file'
        stop
     else
        print*,'Module bonds : SUCCESS - structure_input.log - has been created successfully'
        write(200,'(A)'), 'STRUCTURE VERSION 1.1'
     end if
!
    i = 1   
    if(filename.eq.'auto') then
        call bonding ! subroutine in module-binding.f90
        filename = 'structure_input.txt'
        write(200,'(A)') 'filename of structure input ', trim(filename)
        print*, 'Do you wish to continue or do you want to look at the input file the taht has been autiomatically created &
            ? (c carry on : e exit ) '
        read*, dummy
        if(dummy.eq.'e') then
            stop
        end if
    end if
    
    inquire (file=filename, exist=file_ex)
    if (file_ex.eqv..FALSE.) then
        print*,'File name is ', filename
        stop 'The input file name does not exist please enter a useable filename'
    end if
   

  print*,'Output files are BONDS, ANGLES, DIHEDRALS'
!
! ---------------- Open The Read Files ----------------------------------------------------------------------------------------------
!     
  open (status="old", unit=70, action="read",access="sequential", file=filename, iostat=ierr)               ! Unit 70 user input
  if ( ierr .ne. 0 ) then 
     write(*,'(A)') 'Main : Unable to open ', trim(adjustl(filename)), ', please check file'
     write(200,'(A)') 'Main : Unable to open ', trim(adjustl(filename)), ', please check file'
  else  
     write(*,'(A)') 'Main : ',trim(adjustl(filename)),'has been opened successfully'
     write(200,'(A)') 'Main : ',trim(adjustl(filename)),'has been opened successfully'
  end if

  open (status="old", unit=50, action="read",access="sequential", file="labelled.xyz", iostat=ierr)         ! Unit 50 labelled.xyz
  if ( ierr .ne. 0 ) then 
     write(*,'(A)') 'Main : Unable to open labelled.xyz please check file'
     write(200,'(A)') 'Main : Unable to open labelled.xyz please check file'     
  else  
     write(*,'(A)') 'Main : Labelled.xyz has been opened successfully'
     write(200,'(A)') 'Main : Labelled.xyz has been opened successfully'
  end if
!
! ---------------- Open The Ouput File ----------------------------------------------------------------------------------------------
!
  inquire (file="BONDS.csv", exist=file_ex)
  if (file_ex.eqv..TRUE.) then
     print*, 'Main : WARNING : FILE  BONDS  ALREADY EXISTS.'
     write(200,'(A)') 'Main : WARNING : FILE  BONDS  ALREADY EXISTS. OUTPUT WRITTEN TO fort.40'
!     stop
  else if (file_ex.eqv..FALSE.) then
     open (status="new", unit=40, file="BONDS.csv", action="write",iostat=ierr)                                ! Unit 40 new output file BONDS
     if ( ierr /= 0 ) then 
        write(*,'(A)') 'Unable to create BONDS please check file'
        write(200,'(A)') 'Unable to create BONDS please check file'
     else  
        print*,'Main : BONDS has been created successfully'
        write(200,'(A)') 'Main : BONDS has been created successfully'
     end if
     open (status="replace", unit=44, file="summary_BONDS.csv", action="write",iostat=ierr)                        ! Unit 44 new output file summary BONDS 
     if ( ierr /= 0 ) then 
        write(*,'(A)') 'Unable to create summary_BONDS please check file'
        write(200,'(A)') 'Unable to create summary_BONDS please check file'
     else  
        print*,'Main : summary_BONDS has been created successfully'
        write(200,'(A)') 'Main : summary_BONDS has been created successfully'
     end if
  end if

  inquire (file="ANGLES.csv", exist=file_ex)
  if (file_ex.eqv..TRUE.) then
     print*, 'Main : WARNING : FILE  ANGLES  ALREADY EXISTS.'
     write(200,'(A)') 'Main : WARNING : FILE  ANGLES  ALREADY EXISTS .OUTPUT WRITTEN TO fort.20'
!     stop
  else if (file_ex.eqv..FALSE.) then
     open (status="new", unit=20, file="ANGLES.csv", action="write",iostat=ierr)                               ! Unit 20 new output file ANGLES
     if ( ierr /= 0 ) then 
        write(*,'(A)') 'Main : Unable to create ANGLES please check file'
        write(200,'(A)') 'Main : Unable to create ANGLES please check file'
     else  
        print*,'Main : ANGLES has been created successfully'
        write(200,'(A)') 'Main : ANGLES has been created successfully'
     end if
     open (status="replace", unit=24, file="summary_ANGLES.csv", action="write",iostat=ierr)                        ! Unit 24 new output file summary ANGLES.csv
     if ( ierr /= 0 ) then 
        write(*,'(A)') 'Main : Unable to create summary_ANGLES please check file'
        write(200,'(A)') 'Main : Unable to create summary_ANGLES please check file'
     else  
        print*,'Main : summary_ANGLES has been created successfully'
        write(200,'(A)') 'Main : summary_ANGLES has been created successfully'
     end if
  end if

  inquire (file="DIHEDRALS.csv", exist=file_ex)
  if (file_ex.eqv..TRUE.) then
     print*, 'Main : WARNING : FILE  DIHEDRALS  ALREADY EXISTS.'
     write(200,'(A)') 'Main : WARNING : FILE  DIHEDRALS  ALREADY EXISTS. OUTPUT WRITTEN TO fort.60'
!     stop
  else if (file_ex.eqv..FALSE.) then
     open (status="new", unit=60, file="DIHEDRALS.csv", action="write",iostat=ierr)                            ! Unit 60 new output file DIHEDRALS
     if ( ierr /= 0 ) then 
        write(*,'(A)') 'Unable to create DIHEDRALS please check file'
        write(200,'(A)') 'Unable to create DIHEDRALS please check file'
     else  
        print*,'Main : DIHEDRALS has been created successfully'
        write(200,'(A)') 'Main : DIHEDRALS has been created successfully'
     end if
     open (status="replace", unit=64, file="summary_DIHEDRALS.csv", action="write",iostat=ierr)                    ! Unit 64 new output file summary DIHEDRALS
     if ( ierr /= 0 ) then 
        write(*,'(A)') 'Unable to create summary_DIHEDRALS please check file'
        write(200,'(A)') 'Unable to create summary_DIHEDRALS please check file'
     else  
        print*,'Main : summary_DIHEDRALS has been created successfully'
        write(200,'(A)') 'Main : summary_DIHEDRALS has been created successfully'
     end if

  end if
!
! ---------------- Initialise Number Variables -------------------------------------------------------------------------------------
!
  line_no = 0
  struc = 0
  geometry = 0
  atom_count = 0
  start_atc = 0
  end_atc = 0
  column_no = 0
  status = 0
  status1 = 0
  status2 = 0
  status3 = 0
  FLB1 = .FALSE.
  FLB2 = .FALSE.
  FLB3 = .FALSE.
  FLB4 = .FALSE.
  do i = 1, 3
     atom_coord_1(i) = 0.0
     atom_coord_2(i) = 0.0
     atom_coord_3(i) = 0.0
     atom_coord_4(i) = 0.0
  end do
  bsum = 0
  bond_sq_sum = 0
  Bond_RMSD = 0
  bond_av = 0
  asum = 0
  angle_sq_sum = 0
  Angle_RMSD = 0
  angle_av = 0
  dsum = 0
  dihed_sq_sum = 0
  Dihed_RMSD = 0
  dihed_av = 0
  bmax = 0.0
  bmin = 100
  amax = 0.0
  amin = 360
  dmax = 0.0
  dmin = 360
  bbias = 0
  abias = 0
  Dbias = 0
  bsd = 0
  asd = 0
  Dsd = 0
!
! ---------------- Get Read In Values ---------------------------------------------------------------------------------------------
! THIS LOOP STORES THE XYZ FILE TITLE AND FINDS THE ATOM COUNT
  
  do while (status.eq.0)
     read(50,'(A)',iostat=status) longstring
     ichar=len_trim(longstring)
     line_no=line_no+1
     column_no = index(longstring(1:ichar),atco_search)
     location = index(longstring(1:ichar),search)
     
     if (location.ne.0) then
        struc = struc + 1    ! Counts the number of structures in the stacked xyz file
     end if

     if (line_no.eq.1) then
        title = longstring(1:ichar) ! Stores the title
        ichar_title = ichar
        print*,title
        write(40,"(A)"),title
        write(20,"(A)"),title
        write(60,"(A)"),title
     end if

     if (column_no.ne.0) then ! Print these column heading at the start of each new molecules being found 
        write(40,"(A)"),longstring(1:ichar)
        write(40,"(A)"),'Term, Atoms                   ,,,, Geom     ,       Value'
        write(20,"(A)"),longstring(1:ichar)
        write(20,"(A)"),'Term, Atoms                   ,,,, Geom     ,       Value'
        write(60,"(A)"),longstring(1:ichar)
        write(60,"(A)"),'Term, Atoms                   ,,,, Geom     ,       Value,          positive_sense_angle'
        read(longstring(14:19),"(I5)") atom_count
        print*,'Atom Count ,',atom_count
     end if
  end do
  Rewind(50)
  print*,'Number of time structures - ', struc

  allocate(armsd(1:struc),brmsd(1:struc),drmsd(1:struc),stat=alloc_err) ! allocate memory of the appropriate size
  if (stat.eq.41.or.stat.eq.179.or.stat.eq.632.or.stat.eq.718.or.stat.eq.727) then
     write(200,'(A)')'Main : ERROR : Not enough memory for allocating armsd, brmsd and drmsd.'
     stop 'Main : ERROR : Not enough memory for allocating armsd, brmsd and drmsd.'
  end if
!
! ---------------- Main Loops ------------------------------------------------------------------------------------------------------
!
  line_no = 0
  status = 0
  location = 0
  j = 0
  do while (status.eq.0)
     read(70,'(A)',iostat=status) longstring 

     if (status.gt.0) then
        print*,'Main : ERROR : something went wrong when reading the user input file.'
        write(200,'(A)') 'Main : ERROR : something went wrong when reading the user input file.'
        stop
     else if (status.lt.0) then
        print*,'Main : WARNING : end of user input file reached ensure END is the last line of the file.'
        write(200,'(A)') 'Main : WARNING : end of user input file reached ensure END is the last line of the file.'
        stop
     else
        unit70line_no = unit70line_no +1 ! counts the lines of the users input file this is the atom labels to calculate the BL, BA or DA between
     end if

     ichar=len_trim(longstring)

     if (longstring(1:ichar).eq.'END') then
        print*,'End of input defined by user'
        exit
     end if
     
     calc = adjustl(longstring(:1))
     print*,calc
!
! - - - - - - - - - - - - - - - - - - - - - - - - BONDS - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
     if (calc.eq.bond) then 

        read(longstring(2:),*)  atomlabel1, atomlabel2
        write(44,'(2(A4,A1,1X))') atomlabel1,',', atomlabel2, ','
        icharAL1 = len_trim(atomlabel1)
        icharAL2 = len_trim(atomlabel2)
        print*, atomlabel1, atomlabel2
        blength = 0.0
        
        names = trim(adjustl(trim(adjustl(atomlabel1))//'-'//trim(adjustl(atomlabel2))//'.csv')) ! New ploting file for each bond

        open (status="replace", unit=45, file=trim(adjustl(cwd))//"/BONDS_plot/"//names, action="write", iostat=ierr)                                     ! Unit 45 new output file BONDS plot
        if ( ierr /= 0 ) then
           print*, ierr
           write(*,'(A)') 'Main : WARNING : Unable to create Bonds plot file for ', trim(adjustl(atomlabel1))//&
                trim(adjustl(atomlabel2))
           write(200,'(A)') 'Main : WARNING : Unable to create Bonds plot file for ', trim(adjustl(atomlabel1))//&
                trim(adjustl(atomlabel2))
        else
           write(45,"(A)"),'Term, Atom1, Atom2, Geom, Value'
           write(200,'(A)') 'Main : Main : Bonds plot file for ', trim(adjustl(atomlabel1))//&
                trim(adjustl(atomlabel2)),' created successfully'
        end if

        do while (status1.eq.0)
           line_no = line_no + 1
           read(50,'(A)',iostat=status1) longstringxyz 

           if (status1.gt.0) then
              print*,'Main : ERROR : something went wrong when reading labelled.xyz.'
              write(200,'(A)') 'Main : ERROR : something went wrong when reading labelled.xyz.'
              stop
           else if (status1.lt.0) then
              print*,'Main : WARNING : end of labelled.xyz file reached ensure END is the last line of the file.'
              write(200,'(A)') 'Main : WARNING : end of labelled.xyz file reached ensure END is the last line of the file.'
              stop
           else
              icharxyz=len_trim(longstringxyz)
           end if

           if (longstringxyz(1:icharxyz).eq.'END') then
              print*,'All geometries cycled'
              exit
           end if

           location = index(longstringxyz(1:icharxyz),search)

           if (location.ne.0) then
              start_atc = line_no + 1
              end_atc = start_atc + atom_count
              geometry = geometry + 1
           end if
        
           if (line_no.ge.start_atc.and.line_no.le.end_atc) then ! Look for the atom type between specified line number defining the start of a molecule and the end of a molecule
              location1 = 0
              location1 = index(longstringxyz(1:ichar),atomlabel1(1:icharAL1))
!
! The if loops below store the atoms and the coordinates 
!
              if (location1.ne.0) then
                 FLB1=.TRUE.
                 read(longstringxyz,*) dummy, atom_xyz_1, atom_xyz_2, atom_xyz_3
                 read(atom_xyz_1,'(f20.8)') atom_coord_1(1)
                 read(atom_xyz_2,'(f20.8)') atom_coord_1(2)
                 read(atom_xyz_3,'(f20.8)') atom_coord_1(3)
              end if
              
              location2 = 0
              location2 = index(longstringxyz(1:ichar),atomlabel2(1:icharAL2))
 
              if (location2.ne.0) then
                 FLB2=.TRUE.
                 read(longstringxyz,*) dummy, atom_xyz_1, atom_xyz_2, atom_xyz_3
                 read(atom_xyz_1,'(f20.8)') atom_coord_2(1)
                 read(atom_xyz_2,'(f20.8)') atom_coord_2(2)
                 read(atom_xyz_3,'(f20.8)') atom_coord_2(3)
              end if
              
              if (FLB1.and.FLB2.eqv..TRUE.) then

                 call bondlength(atom_coord_1,atom_coord_2,blength) ! Routine in module-subroutines
                 write(40,'(3(A6),A3,I7,A1,f20.8)'),'Bond, ', atomlabel1//',', atomlabel2,',,,',geometry,',', blength
                 write(45,'(3(A6),A3,I7,A1,f20.8)'),'Bond, ', atomlabel1//',', atomlabel2,',',geometry,',', blength
                 brmsd(geometry) = blength
                 
                 if (blength.gt.bmax)then
                    bmax = blength ! store the largest bond length
                 end if

                 if(blength.lt.bmin) then
                    bmin = blength ! store the smalled bond length
                 end if
                 
                 bsum = bsum + blength ! continuous sum of the bond lengths calculated
                 do i = 1, 3
                    atom_coord_1(i) = 0.0
                    atom_coord_2(i) = 0.0
                 end do

                 FLB1 = .FALSE.
                 FLB2 = .FALSE.
                 blength = 0.0

              end if

           end if
        end do
        
        geometry = 0
        rewind(50)
        status1 = 0

        bond_av = bsum/struc ! Average bond length
        do i = 1,struc
           bbias = bbias + (brmsd(i) - bond_av) ! bias in the bnd lengths
           bsd = bsd + (brmsd(i) - bond_av)**2  ! standard deviation of the bond lengths
           bond_sq_sum = bond_sq_sum + (brmsd(i) - brmsd(1))**2 
        end do
        bbias = bbias/struc
        bsd = sqrt(bsd/(struc - 1))
        Bond_RMSD = sqrt(bond_sq_sum/struc) ! RMSD of the bond lengths relative to the first example
        
        write(40,'(A)')'All units are the same as the input units'
        write(40,'(A16,4X,f20.8)')'RMSD(t;initial),',Bond_RMSD 
        write(40,'(A16,4X,f20.8,1X,A10)')'Normalised RMSD,',Bond_RMSD/(bmax - bmin),',unitless'
        write(40,'(A13,7X,f20.8,1X,A10)')'RMSD/average,',Bond_RMSD/bond_av,',unitless'
        write(40,'(A8,12X,f20.8)')'Average,',bond_av
        write(40,'(A6,14X,3(f20.8,1X,A1,1X))')'Range,',bmin,',',bmax,',', bmax - bmin,''
        write(40,'(A5,15X,f20.8)')'Bias,',bbias
        write(40,'(A19,1X,f20.8)')'Standard deviation,',bsd
        write(40,"(A)"),'------------------- New Bond --------------------'
        write(44,'(A)')'All units are the same as the input units'
        write(44,'(A16,4X,f20.8)')'RMSD(t;initial),',Bond_RMSD 
        write(44,'(A16,4X,f20.8,1X,A10)')'Normalised RMSD,',Bond_RMSD/(bmax - bmin),'unitless'
        write(44,'(A13,7X,f20.8,1X,A10)')'RMSD/average,',Bond_RMSD/bond_av,',unitless'
        write(44,'(A8,12X,f20.8)')'Average,',bond_av
        write(44,'(A6,14X,2(f20.8,1X,A1,1X)A11,1X,f20.8)')'Range,',bmin,',',bmax,',','Difference,', bmax - bmin
        write(44,'(A5,15X,f20.8)')'Bias,',bbias
        write(44,'(A19,1X,f20.8)')'Standard deviation,',bsd
        write(44,'(20X)')
        bond_av = 0.0
        bsum = 0.0
        blength = 0.0
        bond_sq_sum = 0.0
        brmsd = 0.0
        Bond_RMSD = 0.0
        bmax = 0.0
        bmin = 100
        bbias = 0
        bsd = 0
        brmsd = 0.0
        close(45)
!
! - - - - - - - - - - - - - - - - - - - - - - - - ANGLES - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!           
     else if (calc.eq.angle) then
        
        read(longstring(2:),*)  atomlabel1, atomlabel2, atomlabel3
        write(24,'(3(A4,A1,1X))') atomlabel1,',', atomlabel2, ',',atomlabel3, ','
        icharAL1 = len_trim(atomlabel1)
        icharAL2 = len_trim(atomlabel2)
        icharAL3 = len_trim(atomlabel3)
        print*, atomlabel1, atomlabel2, atomlabel3

        names = trim(adjustl(trim(adjustl(atomlabel1))//'-'//trim(adjustl(atomlabel2))//'-'//trim(adjustl(atomlabel3))//'.csv')) ! New ploting file for each angle

        open (status="replace", unit=25, file=trim(adjustl(cwd))//"/ANGLES_plot/"//names, action="write", iostat=ierr)                                     ! Unit 25 new output file ANGLES plot
        if ( ierr /= 0 ) then
           print*, ierr
           write(*,'(A)') 'Main : WARNING : Unable to create Bonds plot file for ', trim(adjustl(atomlabel1))//&
                trim(adjustl(atomlabel2))//trim(adjustl(atomlabel3))
           write(200,'(A)') 'Main : WARNING : Unable to create Bonds plot file for ', trim(adjustl(atomlabel1))//&
                trim(adjustl(atomlabel2))//trim(adjustl(atomlabel3))
        else
           write(25,"(A)"),'Term, Atom1, Atom2, Atom3, Geom, Value'
           write(200,'(A)') 'Main : Main : Bonds plot file for ', trim(adjustl(atomlabel1))//&
                trim(adjustl(atomlabel2))//trim(adjustl(atomlabel3)),' created successfully'
        end if

        do while (status2.eq.0)
           line_no = line_no + 1
           read(50,'(A)',iostat=status2) longstringxyz

           if (status2.gt.0) then
              print*,'Main : ERROR : something went wrong when reading labelled.xyz.'
              write(200,'(A)') 'Main : ERROR : something went wrong when reading labelled.xyz at the start of angles.'
              stop
           else if (status2.lt.0) then
              print*,'Main : WARNING : end of labelled.xyz file reached ensure END is the last line of the file.'
              write(200,'(A)') 'Main : WARNING : end of labelled.xyz file reached ensure END is the last line of the file.'
              stop
           else
              icharxyz=len_trim(longstringxyz)
           end if

           if (longstringxyz.eq.'END') then
              print*,'All geometries cycled'
              exit
           end if

           location = index(longstringxyz(1:icharxyz),search)

           if (location.ne.0) then
              start_atc = line_no + 1
              end_atc = start_atc + atom_count
              geometry = geometry + 1
           end if
!
! The if loops below find and store the atom types specified and the coordinates of the atoms 
!
           if (line_no.ge.start_atc.and.line_no.le.end_atc) then
              location1 = 0
              location1 = index(longstringxyz(1:ichar),atomlabel1(1:icharAL1))

              if (location1.ne.0) then
                 FLB1=.TRUE.
                 read(longstringxyz,*) dummy, atom_xyz_1, atom_xyz_2, atom_xyz_3
                 read(atom_xyz_1,'(f20.8)') atom_coord_1(1)
                 read(atom_xyz_2,'(f20.8)') atom_coord_1(2)
                 read(atom_xyz_3,'(f20.8)') atom_coord_1(3)
              end if
              
              location2 = 0
              location2 = index(longstringxyz(1:ichar),atomlabel2(1:icharAL2))                       
              
              if (location2.ne.0) then
                 FLB2=.TRUE.
                 read(longstringxyz,*) dummy, atom_xyz_1, atom_xyz_2, atom_xyz_3
                 read(atom_xyz_1,'(f20.8)') atom_coord_2(1)
                 read(atom_xyz_2,'(f20.8)') atom_coord_2(2)
                 read(atom_xyz_3,'(f20.8)') atom_coord_2(3)
              end if
              
              location3 = 0
              location3 = index(longstringxyz(1:ichar),atomlabel3(1:icharAL3))

              if (location3.ne.0) then
                 FLB3=.TRUE.
                 read(longstringxyz,*) dummy, atom_xyz_1, atom_xyz_2, atom_xyz_3
                 read(atom_xyz_1,'(f20.8)') atom_coord_3(1)
                 read(atom_xyz_2,'(f20.8)') atom_coord_3(2)
                 read(atom_xyz_3,'(f20.8)') atom_coord_3(3)
              end if
              
              if (FLB1.and.FLB2.and.FLB3.eqv..TRUE.) then

                 call bondangle(atom_coord_1,atom_coord_2,atom_coord_3,bangle,bangledeg) ! subroutine in module-subroutines.f90
                 if (bangledeg.gt.amax)then
                    amax = bangledeg ! biggest angle
                 end if

                 if(bangledeg.lt.amin) then
                    amin = bangledeg ! smallest angle
                 end if
                 
                 write(20,'(4(A6),A2,I7,A1,f20.8)'),'angle, ', atomlabel1//',', atomlabel2//',',atomlabel3,',,',geometry,','&
                      , bangledeg
                 write(25,'(4(A6),A2,I7,A1,f20.8)'),'angle, ', atomlabel1//',', atomlabel2//',',atomlabel3,',',geometry,','&
                      , bangledeg
                 armsd(geometry) = bangledeg 
                 asum = asum + bangledeg
                 do i = 1, 3
                    atom_coord_1(i) = 0.0
                    atom_coord_2(i) = 0.0
                    atom_coord_3(i) = 0.0
                 end do

                 FLB1 = .FALSE.
                 FLB2 = .FALSE.
                 FLB3 = .FALSE.
                 bangle = 0.0
                 bangledeg = 0.0

              end if

           end if
        end do

        rewind(50)
        geometry = 0
        status2 = 0

        angle_av = asum/struc ! mean angle
        do i = 1,struc
           abias = abias + (armsd(i) - angle_av) ! angle bias
           asd = asd + (armsd(i) - angle_av)**2 ! angle standard deviation
           angle_sq_sum = angle_sq_sum + (armsd(i) - armsd(1))**2 
        end do
        abias = abias/struc
        asd = sqrt(asd/struc)
        Angle_RMSD = sqrt(angle_sq_sum/struc) ! angle RMSD

        write(20,'(A)')'All units are the same as the input units,'
        write(20,'(A16,4X,f20.8,1X)')'RMSD(t;initial),',Angle_RMSD 
        write(20,'(A16,4X,f20.8,1X,A10)')'Normalised RMSD,',Angle_RMSD/(amax - amin),',unitless,'
        write(20,'(A13,7X,f20.8,1X,A10)')'RMSD/average,',Angle_RMSD/angle_av,',unitless,'
        write(20,'(A8,12X,f20.8)')'Average,',angle_av
        write(20,'(A6,14X,3(f20.8,1X,A1,1X))')'Range,',amin,',',amax,',', amax - amin, ','
        write(20,'(A5,15X,f20.8)')'Bias,',abias
        write(20,'(A19,1X,f20.8)')'Standard deviation,',asd
        write(20,"(A)"),'------------------- New Angle --------------------'
        write(24,'(A)')'All units are the same as the input units,'
        write(24,'(A16,4X,f20.8,1X)')'RMSD(t;initial),',Angle_RMSD 
        write(24,'(A16,4X,f20.8,1X,A10)')'Normalised RMSD,',Angle_RMSD/(amax - amin),',unitless,'
        write(24,'(A13,7X,f20.8,1X,A10)')'RMSD/average,',Angle_RMSD/angle_av,',unitless,'
        write(24,'(A8,12X,f20.8)')'Average,',angle_av
        write(24,'(A6,14X,2(f20.8,1X,A1,1X)A11,1X,f20.8)')'Range,',amin,',',amax,',','Difference,', amax - amin
        write(24,'(A5,15X,f20.8)')'Bias,',abias
        write(24,'(A19,1X,f20.8)')'Standard deviation,',asd
        write(24,'(20X)')
        angle_av = 0.0
        asum = 0.0
        bangledeg = 0.0
        angle_sq_sum = 0.0
        armsd = 0.0
        Angle_RMSD = 0.0
        amax = 0
        amin = 360
        asd = 0
        abias = 0
        armsd = 0
        close(25)
!
! - - - - - - - - - - - - - - - - - - - - - - - - DIHEDRALS - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
!
     else if (calc.eq.dihed) then
        read(longstring(2:),*)  atomlabel1, atomlabel2, atomlabel3, atomlabel4
        write(64,'(4(A4,A1,1X))') atomlabel1,',', atomlabel2, ',',atomlabel3, ',',atomlabel4, ','
        icharAL1 = len_trim(atomlabel1)
        icharAL2 = len_trim(atomlabel2)
        icharAL3 = len_trim(atomlabel3)
        icharAL4 = len_trim(atomlabel4)
        print*, atomlabel1, atomlabel2, atomlabel3, atomlabel4

        names = trim(adjustl(trim(adjustl(atomlabel1))//'-'//trim(adjustl(atomlabel2))//'-'//trim(adjustl(atomlabel3))//'-'&
             //trim(adjustl(atomlabel4))//'.csv'))! New ploting file for each dihedral

        open (status="replace", unit=65, file=trim(adjustl(cwd))//"/DIHEDRAL_plot/"//names, action="write", iostat=ierr)                                     ! Unit 65 new output file DIHEDRAL plot
        if ( ierr /= 0 ) then
           print*, ierr
           write(*,'(A)') 'Main : WARNING : Unable to create Bonds plot file for ', trim(adjustl(atomlabel1))//&
                trim(adjustl(atomlabel2))//trim(adjustl(atomlabel3))//trim(adjustl(atomlabel4))
           write(200,'(A)') 'Main : WARNING : Unable to create Bonds plot file for ', trim(adjustl(atomlabel1))//&
                trim(adjustl(atomlabel2))//trim(adjustl(atomlabel3))//trim(adjustl(atomlabel4))
        else
           write(65,"(A)"),'Term, Atom1, Atom2, Atom3, Atom3, Geom, Value, Positive_sense_angle'
           write(200,'(A)') 'Main : Main : Bonds plot file for ', trim(adjustl(atomlabel1))//&
                trim(adjustl(atomlabel2))//trim(adjustl(atomlabel3))//trim(adjustl(atomlabel4)),' created successfully'
        end if

        do while (status3.eq.0)
           line_no = line_no + 1
           read(50,'(A)',iostat=status3) longstringxyz

           if (status3.gt.0) then
              print*,'Main : ERROR : something went wrong when reading labelled.xyz.'
              write(200,'(A)')'Main : ERROR : something went wrong when reading labelled.xyz at the start of the dihedrals.'
              stop
           else if (status3.lt.0) then
              print*,'Main : WARNING : end of labelled.xyz file reached ensure END is the last line of the file.'
              write(200,'(A)') 'Main : WARNING : end of labelled.xyz file reached ensure END is the last line of the file.'
              stop
           else
              icharxyz=len_trim(longstringxyz)
           end if

           if (longstringxyz.eq.'END') then
              print*,'All geometries cycled'
              exit
           end if

           location = index(longstringxyz(1:icharxyz),search)

           if (location.ne.0) then
              start_atc = line_no + 1
              end_atc = start_atc + atom_count
              geometry = geometry + 1
           end if
        
           if (line_no.ge.start_atc.and.line_no.le.end_atc) then
              location1 = 0
              location1 = index(longstringxyz(1:ichar),atomlabel1(1:icharAL1))
!
! The if loops below store the atoms and its coordinates when they are found
!
              if (location1.ne.0) then
                 FLB1=.TRUE.
                 read(longstringxyz,*) dummy, atom_xyz_1, atom_xyz_2, atom_xyz_3
                 read(atom_xyz_1,'(f20.8)') atom_coord_1(1)
                 read(atom_xyz_2,'(f20.8)') atom_coord_1(2)
                 read(atom_xyz_3,'(f20.8)') atom_coord_1(3)
              end if
              
              location2 = 0
              location2 = index(longstringxyz(1:ichar),atomlabel2(1:icharAL2))                    
              
              if (location2.ne.0) then
                 FLB2=.TRUE.
                 read(longstringxyz,*) dummy, atom_xyz_1, atom_xyz_2, atom_xyz_3
                 read(longstringxyz,*) dummy, atom_xyz_1, atom_xyz_2, atom_xyz_3
                 read(atom_xyz_1,'(f20.8)') atom_coord_2(1)
                 read(atom_xyz_2,'(f20.8)') atom_coord_2(2)
                 read(atom_xyz_3,'(f20.8)') atom_coord_2(3)
              end if
              
              location3 = 0
              location3 = index(longstringxyz(1:ichar),atomlabel3(1:icharAL3))

              if (location3.ne.0) then
                 FLB3=.TRUE.
                 read(longstringxyz,*) dummy, atom_xyz_1, atom_xyz_2, atom_xyz_3
                 read(atom_xyz_1,'(f20.8)') atom_coord_3(1)
                 read(atom_xyz_2,'(f20.8)') atom_coord_3(2)
                 read(atom_xyz_3,'(f20.8)') atom_coord_3(3)
              end if

              location4 = 0
              location4 = index(longstringxyz(1:ichar),atomlabel4(1:icharAL4))

              if (location4.ne.0) then
                 FLB4=.TRUE.
                 read(longstringxyz,*) dummy, atom_xyz_1, atom_xyz_2, atom_xyz_3
                 read(atom_xyz_1,'(f20.8)') atom_coord_4(1)
                 read(atom_xyz_2,'(f20.8)') atom_coord_4(2)
                 read(atom_xyz_3,'(f20.8)') atom_coord_4(3)
              end if
              
              if (FLB1.and.FLB2.and.FLB3.and.FLB4.eqv..TRUE.) then

                 call dihedral(atom_coord_1,atom_coord_2,atom_coord_3,atom_coord_4,dangle,dangledeg) ! subroutine in module-subroutines.f90
!
! Determines the Dihedral in a consistent frame/direction                   
!
                 if (dangledeg.lt.0) then
                       possensedihed = 360 - abs(dangledeg)
                 else 
                       possensedihed = dangledeg
                 end if
 
                 if (possensedihed.gt.Dmax)then
                    Dmax = possensedihed ! biggest dihedral 
                 end if

                 if(possensedihed.lt.Dmin) then
                    Dmin = possensedihed ! smallest dihedral
                 end if

                 write(60,'(5(A6),A1,I6,A1,f20.8,A1,f20.8)'),'Dihed, ', atomlabel1//',', atomlabel2//',',atomlabel3//',',&
                 atomlabel4,',' ,geometry,',', dangledeg, ',', possensedihed
                 write(65,'(5(A6),A1,I6,A1,f20.8,A1,f20.8)'),'Dihed, ', atomlabel1//',', atomlabel2//',',atomlabel3//',',&
                 atomlabel4,',' ,geometry,',', dangledeg, ',', possensedihed
                 drmsd(geometry) = possensedihed
                 dsum = dsum + possensedihed
                 do i = 1, 3
                    atom_coord_1(i) = 0.0
                    atom_coord_2(i) = 0.0
                    atom_coord_3(i) = 0.0
                    atom_coord_4(i) = 0.0
                 end do

                 FLB1 = .FALSE.
                 FLB2 = .FALSE.
                 FLB3 = .FALSE.
                 FLB4 = .FALSE.
                 dangle = 0.0
                 dangledeg = 0.0

              end if
           end if
        end do

        rewind(50)
        geometry = 0
        status3 = 0

        dihed_av = dsum/struc ! mean dihedral
        do i = 1,struc
           Dbias = Dbias + (drmsd(i) - dihed_av) ! bias dihedral
           Dsd = Dsd + (drmsd(i) - dihed_av)**2 ! standard deviation dihedral
           dihed_sq_sum = dihed_sq_sum + (drmsd(i) - drmsd(1))**2
        end do
        Dbias = Dbias/struc 
        Dsd = sqrt(Dsd/struc)
        Dihed_RMSD = sqrt(dihed_sq_sum/struc) ! RMSD dihedral

        write(60,'(A)')'All units are the same as the input units,'
        write(60,'(A16,4X,f20.8)'),'RMSD(t;inital),',Dihed_RMSD 
        write(60,'(A16,4X,f20.8,1X,A10)')'Normalised RMSD,',Dihed_RMSD/(Dmax - Dmin),',unitless,'
        write(60,'(A13,7X,f20.8,1X,A10)')'RMSD/average,',Dihed_RMSD/dihed_av,',unitless,'
        write(60,'(A8,12X,f20.8)')'Average,',dihed_av
        write(60,'(A6,14X,2(f20.8,1X,A1,1X))')'Range,',Dmin,',',Dmax,',', Dmax - Dmin, ','
        write(60,'(A5,15X,f20.8)')'Bias,',Dbias
        write(60,'(A19,1X,f20.8)')'Standard deviation,',Dsd
        write(60,"(A)"),'------------------- New Dihedrl --------------------'
        write(64,'(A)')'All units are the same as the input units,'
        write(64,'(A16,4X,f20.8)'),'RMSD(t;inital),',Dihed_RMSD 
        write(64,'(A16,4X,f20.8,1X,A10)')'Normalised RMSD,',Dihed_RMSD/(Dmax - Dmin),',unitless,'
        write(64,'(A13,7X,f20.8,1X,A10)')'RMSD/average,',Dihed_RMSD/dihed_av,',unitless,'
        write(64,'(A8,12X,f20.8)')'Average,',dihed_av
        write(64,'(A6,14X,2(f20.8,1X,A1,1X)A11,1X,f20.8)')'Range,',Dmin,',',Dmax,',','Difference,', Dmax - Dmin
        write(64,'(A5,15X,f20.8)')'Bias,',Dbias
        write(64,'(A19,1X,f20.8)')'Standard deviation,',Dsd
        write(64,'(5X)')

        dihed_av = 0.0
        dsum = 0.0
        possensedihed = 0.0
        dihed_sq_sum = 0.0
        drmsd = 0.0
        Dihed_RMSD = 0.0
        Dmin = 360
        Dmax = 0
        Dbias = 0
        Dsd = 0
        drmsd = 0
        close(65)
!
! - - - - - - - - - - - - - - - - - - - - - - - - ERROR - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
     else
        Print*, 'Main : WARNING : unrecognised calculation requested'
        write(200,'(A)') 'Main : WARNING : unrecognised calculation requested'
        print*, longstring(2:)
        cycle
     end if

  end do
!
  close(50)
  close(40)
  close(70)
  close(20)
  close(60)
  close(24)
  close(44)
  close(64)

  if(write_norm_xyz.eq."no") then
     print*, 'DONE'
  else if(write_norm_xyz.ne."no") then
     call create_norm_xyz(atom_count,fil_nam,write_norm_xyz,title,ichar_title) ! subroutine in module-subroutines.f90 creates a jmol plotable image 
  else
     print*,'Main : ERROR : no recognisable option to use with writing normalised xyz'
     write(200,'(A)') 'Main : ERROR : no recognisable option to use with writing normalised xyz'
  end if

  close(200)

end subroutine  stucture

end module structurechanges
  
  

  
  
