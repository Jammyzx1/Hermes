!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Originally created by James McDonagh at the University of Manchester 2015, in the Popelier group      ! 
! Licensed under Attribution-ShareAlike 2.5 Generic (CC BY-SA 2.5)                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Bonding is assumed within a distance cutoff so this routine is not fool proof. Heavy atom bond if
! distance > 0.4 A but < 2.0 A. If it H then bonded if the distance is  > 0.4 A but < 1.2 A
! (A = Angstroms).
!
! Originally created by James.
! Version 1.1
! CHANGE LOG 
! Version 1  : Subroutine for determining atoms which are bonded, in an angle and in a dihedral angle
! Version 1.1: Modular format and incorporation in the Hermes program.
! Version 2  : This subroutine replaces the messy module_bonding.f90 for improved clarity. This code
!              does not attempt to write only unique bond, angle and dihedral information due to concerns
!              lost information from doing so with module_bonding.f90.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module module_structure_determine
  
  use module_subroutines

  implicit none

contains

  subroutine bonding_determine
    integer :: limit = 500                                                                                       ! used to size the array looking for replicated angles and dihedrals
    character (len = 10) :: atco_search = "Atom Count"                                                           ! regular expression to find the atom count line
    character (len = 11) :: search = "-----------"                                                               ! string to find for each new geometry
    character (len = 5), dimension(:), allocatable  :: atoms
    character (len = 30), dimension(:), allocatable :: angle_list                                                ! array to store previously found angles and dihedra
    character (len = 256) :: longstring                                                                          
    character (len = 19) :: output_filename = "structure_input.txt"
    character (len = 30) :: angle, dihedrals, at
    real :: coord_1(3), coord_2(3), coord_3(3), bond_length, bond_length1
    real, allocatable :: coords(:,:)
    real, allocatable :: bonds(:,:,:)
    real, dimension(5) :: temp
    real, allocatable :: structure(:,:)
    character (len = 20) :: atom_xyz_1, atom_xyz_2, atom_xyz_3
    integer :: ichar, line_no, atom_count, atom_count_line, stc, edc, j, i, n, m, k, l, x, o, ierr, start_structure, status, store
    integer :: hostatom, storepoint, hostatom2, storepoint2, hydrogen, hydrogen2, label, a, b, ilabel, i2label, angl, replic
    integer :: angle_list_line, a1, hydrogen3
    integer :: a2, a3, a4, reang, a123, a1234, var
    logical :: file_ex, swap 
!    
    print*,'This subroutine automatically assigns bonds on the basis of distance, this is by no means fool proof.'
    print*,'CRITERIA FOR A BOND: BOND LENGTH > 0.4 ANGSTROMS < 2.0 ANGSTROMS. IF H THEN BOND LENGTH > 0.4 ANGSTROMS < 1.2 ANGSTROMS'
    print*,'Based on this list all possible unique angles and dihedrals are defined.'
    print*,'Please check the file input.txt.' 
    print*,'If it goes wrong dont blame me.'
    print*,'                                                    '
!
!---------------------- Initalisation ------------------------------------------
!
    i = 0
    j = 0
    start_structure = 0
    status = 0
    line_no = 0
    ichar = 0
    atom_count = 0
    stc = 0
    edc = 0
    n = 0
    m = 0
    k = 0
    l = 0
    angle_list_line = 0
    reang = 0
    a123 = 0
!    
!--------------------------- Open labelled.xyz ---------------------------------
!
    open (status="old", unit=50, action="read",access="sequential", file="labelled.xyz", iostat=ierr)         ! Unit 50 labelled.xyz
    if ( ierr .ne. 0 ) then
       write(*,'(A)')' Module bonds : ERROR - Unable to open labelled.xyz please check file'
       write(200,'(A)')' Module bonds : ERROR - Unable to open labelled.xyz please check file'
    else
       write(*,'(A)')' Module bonds : SUCCESS - Labelled.xyz has been opened successfully'
       write(200,'(A)')' Module bonds : SUCCESS - Labelled.xyz has been opened successfully'
    end if
!
!------------------------- Open Output -----------------------------------------
!
    open (status="replace", unit=160, file=output_filename, action="write",iostat=ierr)                       ! Unit 160 new output file output_filename defualt structure_input.txt
    if ( ierr.ne.0 ) then
       write(200,'(A)')'Module bonds : ERROR - unable to create ',output_filename,' please check file'
       write(*,'(A)')'Module bonds : ERROR - unable to create ',output_filename,' please check file'
       stop 
    end if
 
!   
!------------------- Read the first geometry -----------------------------------
!   
  do while (status.eq.0)
     read(50,'(A)',iostat=status) longstring
     ichar=len_trim(longstring)
     line_no=line_no+1
     atom_count_line = index(longstring(1:ichar),atco_search)
!
     if(atom_count_line.ne.0) then
        read(longstring(14:19),"(I5)") atom_count
        print*, atom_count
        allocate(atoms(1:atom_count), structure(1:atom_count,1:(6+atom_count)))
        allocate(angle_list(1:limit))
        structure = 0.0
        print*, 'allocated'
     end if
!
     start_structure = index(longstring(1:ichar),search) ! Search for regex that specifies the start of a structure
     if (start_structure.ne.0) then                      ! if the start of a structure is found 
        i = i + 1                                        ! Count the structures found 
        if (i.gt.1) then                                 ! If more than one i.e. you go past the first exit all deterination are made based on the first strutcure
           exit
        end if
        print*, atom_count
        stc = line_no
        edc = line_no+atom_count
        print*,'Start - ', stc                           ! Print and write the expected start and end lines 
        write(200,'(A19,1X,I2)') 'Start read line - ', stc
        print*,'End line - ', edc
        write(200,'(A19,1X,I2)') 'End read line - ', edc
     end if
!     
     if (line_no.gt.stc.and.line_no.le.edc) then         ! If in the molecular geometry is being specified then store the atomic coordinates
        j = j + 1
        if (j.gt.atom_count) then
           print*, 'Module bonds : WARNING - coordinates exceed the atom count, number of atoms mismatch'
           write(200,'(A)') 'Module bonds : WARNING - coordinates exceed the atom count, number of atoms mismatch'
        end if
        read(longstring(1:ichar),*) atoms(j), atom_xyz_1, atom_xyz_2, atom_xyz_3
        structure(j,1) = j                               ! Structure as an array holds the atom index x y z count of bonded parteners list of bonded parteners. 
        read(atom_xyz_1,'(f20.8)') structure(j,2)        ! Atoms has the atom labels stored in the order they are read
        read(atom_xyz_2,'(f20.8)') structure(j,3)
        read(atom_xyz_3,'(f20.8)') structure(j,4) 
     end if
  end do
  print*, 'Atoms array'
  print*, atoms
!  
!------------------------- Calculate bonded atoms  -----------------------------
!
  do i = 1, atom_count
     structure(i,2) = coord_1(1) ! Coord is a vector of x y z for the first atom (reference atom) 
     structure(i,3) = coord_1(2)
     structure(i,4) = coord_1(3)
     n = int(structure(i,1))
     write(200,'(A)') '----- Atom ', atoms(n),' -----'
     do j = 1, atom_count
        structure(j,2) = coord_2(1) ! Coord is a vector of x y z for the second atom (test atom) 
        structure(j,3) = coord_2(2)
        structure(j,4) = coord_2(3)
        m = int(structure(j,1))

        hydrogen = index(atoms(n),'H')
        hydrogen2 = index(atoms(j),'H')
        
        call bondlength(coord_1, coord_2, bond_length)

         if(hydrogen.ne.0.or.hydrogen2.ne.0) then
            if (bond_length.ge.0.4.and.bond_length.le.1.2) then
               structure(i,5) = structure(i,5) + 1
               write(200,'(A)') 'Atoms ', atoms(int(structure(i,1))), ' and ', atoms(int(structure(j,1))), ' Bond'
            end if
            write(200,'(A)') 'Atoms ', atoms(int(structure(i,1))), ' and ', atoms(int(structure(j,1))), ' do not bond'
         else 
            if (bond_length.ge.0.4.and.bond_length.le.2.0) then
               structure(i,5) = structure(i,5) + 1
               write(200, '(A)') 'Atoms ', atoms(int(structure(i,1))), ' and ', atoms(int(structure(j,1))), ' Bond'
            end if
            write(200, '(A)') 'Atoms ', atoms(int(structure(i,1))), ' and ', atoms(int(structure(j,1))), ' do not bond'
         endif

      end do
   end do
!  
!------------------------- Sort the list in order of bonds  --------------------
!
   print*, 'Structure array contents: '
   do i = 1, atom_count
      print*, atoms(i), structure(i,1), structure(i,2), structure(i,3), structure(i,4), structure(i,5)
   end do

   do i = 1, atom_count 
      print*,'Before sort: Atom number ',structure(i,1),' Atom label ', atoms(i),' No. of bonded atoms ', structure(i,5)
   end do
!
!- - - - - Sort 1 Num of bonded atoms - - - - - Bubble sort quick enough for small molecules change for quick sort to improve 
!
   do j = atom_count-1, 1, -1
      swap = .FALSE.
      do i = 1, j
         if (structure(i,5).gt.structure(i + 1,5)) then
            
            temp(1) = structure(i,1)
            temp(2) = structure(i,2)
            temp(3) = structure(i,3)
            temp(4) = structure(i,4)
            temp(5) = structure(i,5)

            structure(i,1) = structure(i+1,1)
            structure(i,2) = structure(i+1,2)
            structure(i,3) = structure(i+1,3)
            structure(i,4) = structure(i+1,4)
            structure(i,5) = structure(i+1,5)

            structure(i+1,1) = temp(1)
            structure(i+1,2) = temp(2)
            structure(i+1,3) = temp(3)
            structure(i+1,4) = temp(4)
            structure(i+1,5) = temp(5)
            
            swap = .TRUE.
            
         end if
      end do
      if (.NOT. swap) exit
   end do
!
!- - - - - Sort 2 Hydrogen first - - - - - Bubble sort quick enough for small molecules change for quick sort to improve
!
   do j = atom_count-1, 1, -1
      swap = .FALSE.
      do i = 1, j
         if (index(atoms(i),'H').eq.0.and.index(atoms(i+1),'H').ne.0) then
!                        
            temp(1) = structure(i,1)
            temp(2) = structure(i,2)
            temp(3) = structure(i,3)
            temp(4) = structure(i,4)
            temp(5) = structure(i,5)
!
            structure(i,1) = structure(i+1,1)
            structure(i,2) = structure(i+1,2)
            structure(i,3) = structure(i+1,3)
            structure(i,4) = structure(i+1,4)
            structure(i,5) = structure(i+1,5)
!
            structure(i+1,1) = temp(1)
            structure(i+1,2) = temp(2)
            structure(i+1,3) = temp(3)
            structure(i+1,4) = temp(4)
            structure(i+1,5) = temp(5)
!
            swap = .TRUE.
!            
         end if
      end do
      if (.NOT. swap) exit
   end do
   temp = 0.0
!
   print*,'Sorted list'
   do i = 1, atom_count
      label = int(bonds(i,1,4))  
      print*,'After sort: Atom number ',bonds(i,1,4),' Atom label ', atoms(label),' No. of bonded atoms ', bonds(i,1,5)
   end do
!  
!------------------------- Calculate bonded atoms  -----------------------------
!
  do i = 1, atom_count
     coord_1(1) = structure(i,2) ! Coord is a vector of x y z for the first atom (reference atom) 
     coord_1(2) = structure(i,3)
     coord_1(3) =  structure(i,4)
     n = int(structure(i,1))
     do j = 1, atom_count
        coord_2(1) = structure(j,2) ! Coord is a vector of x y z for the second atom (test atom) 
        coord_2(2) = structure(j,3)
        coord_2(3) = structure(j,4)
        m = int(structure(j,1))
!
        hydrogen = index(atoms(n),'H')
        hydrogen2 = index(atoms(j),'H')
!        
        call bondlength(coord_1, coord_2, bond_length)
!
        ! In these loops we store the unique atom bonding combinations to do this if the atoms meet the bonding rules then the index of the test atom is stored in the array structure
        ! after the reference atom data is stored (after the second structure array index >= 7). In principle this data is all we need to work out bonds, angles and dihedrals. Below
        ! this is done in the second line of the if loops with the second structure index given by adding one to the current bonding atom count from the reference atom (index 6)
        ! which is set in the first line after the if loops here. This atom count can be compared with the first atom count to check the sorting hasn't made a change.
!
         if(hydrogen.ne.0.or.hydrogen2.ne.0) then
            if (bond_length.ge.0.4.and.bond_length.le.1.2) then
               structure(i,6) = structure(i,6) + 1
               structure(i,(int(structure(i,6)) + 1)) = j
            end if
         else
            if (bond_length.ge.0.4.and.bond_length.le.2.0) then
               structure(i,6) = structure(i,6) + 1
               structure(i,(int(structure(i,6)) + 1)) = j
            end if
         endif
      end do
   end do
!  
!------------------------- bond length to output -------------------------------
!
   do i = 1, atom_count ! write to output files warning if there is a bonded atom mismatch and a bonded atom list to the structure input file
      if(structure(i,5).ne.structure(i,6)) then
         write(200,'(A)') 'WARNING - atom ', atoms(i) , ' has a mismatch in number of bonding atoms after sort please check the &
              output is as expected for that atom'
      end if
      do j = 7,(6 + int(structure(i,6)))
         write(160,'(A)')'B ', atoms(int(structure(i,1))), atoms(int(structure(i,j)))
      end do
   end do
!  
!------------------------- Write molecular angle -------------------------------
!
   do i = 1, atom_count ! write to output files warning if there is a bonded atom mismatch and a bonded atom list to the structure input file
      if(structure(i,5).ne.structure(i,6)) then
         write(200,'(A)') 'WARNING - atom ', atoms(i) , ' has a mismatch in number of bonding atoms after sort please check the &
              output is as expected for that atom'
      end if
      do j = 7,(6 + int(structure(i,6)))
         do k = 7, (6 + int(structure(j,6)))
            write(160,'(A)')'A ', atoms(int(structure(i,1))), atoms(int(structure(i,j))), atoms(int(structure(j,k)))
         end do
      end do
   end do
!  
!------------------------- Write molecular dihedral ----------------------------
!
   do i = 1, atom_count ! write to output files warning if there is a bonded atom mismatch and a bonded atom list to the structure input file
      if(structure(i,5).ne.structure(i,6)) then
         write(200,'(A)') 'WARNING - atom ', atoms(i) , ' has a mismatch in number of bonding atoms after sort please check the &
              output is as expected for that atom'
      end if
      do j = 7,(6 + int(structure(i,6)))
         do k = 7, (6 + int(structure(j,6)))
            do l = 7, (6 + int(structure(k,6)))
               write(160,'(A)')'D ', atoms(int(structure(i,1))), atoms(int(structure(i,j))), atoms(int(structure(j,k))), &
                    atoms(int(structure(k,l)))
            end do
         end do
      end do
   end do
!  
!------------------------- Close the subroutine --------------------------------
!
 end subroutine bonding_determine
   

 
         
