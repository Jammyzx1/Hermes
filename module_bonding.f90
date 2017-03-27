!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Originally created by James McDonagh at the University of Manchester 2015, in the Popelier group      ! 
! Licensed under MIT License                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This is a messy subroutine that could be improved, but it is quick. Bonding is assumed within a
! distance cutoff so this routine is not fool proof. Heavy atom bond if distance > 0.4 A but < 2.0 A. 
! If it H then bonded if the distance is  > 0.4 A but < 1.2 A (A = Angstroms).
!
! Originally created by James.
! Version 1.1
! CHANGE LOG 
! Version 1  : Subroutine for determining atoms which are bonded, in an angle and in a dihedral angle
! Version 1.1: Modular format and incorporation in the Hermes program.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module module_bonding
  
  use module_subroutines

  implicit none

contains

  subroutine bonding
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
    open (status="replace", unit=160, file=output_filename, action="write",iostat=ierr)                          ! Unit 160 new output file output_filename defualt structure_input.txt
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
        allocate(atoms(1:atom_count),coords(1:atom_count,1:3),bonds(1:atom_count,1:atom_count,1:5))
        allocate(angle_list(1:limit))
        print*, 'allocated'
     end if
!
     start_structure = index(longstring(1:ichar),search)
     if (start_structure.ne.0) then
        i = i + 1
        if (i.gt.1) then
           exit
        end if
        print*, atom_count
        stc = line_no
        edc = line_no+atom_count
        print*,'Start - ', stc
        write(200,'(A19,1X,I2)') 'Start read line - ', stc
        print*,'End line - ', edc
        write(200,'(A19,1X,I2)') 'End read line - ', edc
     end if
!     
     if (line_no.gt.stc.and.line_no.le.edc) then
        j = j + 1
        if (j.gt.atom_count) then
           print*, 'Module bonds : WARNING - coordinates exceed the atom count, number of atoms mismatch'
           write(200,'(A)') 'Module bonds : WARNING - coordinates exceed the atom count, number of atoms mismatch'
        end if
        read(longstring(1:ichar),*) atoms(j), atom_xyz_1, atom_xyz_2, atom_xyz_3
        read(atom_xyz_1,'(f20.8)') coords(j,1)
        read(atom_xyz_2,'(f20.8)') coords(j,2)
        read(atom_xyz_3,'(f20.8)') coords(j,3)
     end if
  end do
  print*, 'Atoms array'
  print*, atoms
!  
!------------------------- Initalise bonds array -------------------------------
!
  bonds = 0.0
!
!------------------------ Calculate the atom atom distances --------------------
!
  hostatom = 1
  storepoint = 1
  do i = 1,atom_count
!        
     if(hostatom.ne.i) then
        write(200,'(A3,A17,I3)') atoms(i-1),'No. bonded atoms ',storepoint - 1
        bonds(hostatom,1,5) = storepoint - 1
        hostatom = hostatom + 1                        ! JMcD Count of atoms held constant while other atom is varied
        storepoint = 1                                 ! JMcD The array index to store the label of bonded atoms to the host atom in storepoint - 1 is the number of bonded atoms to the host atom
     end if

     do n = 1,3
        coord_1(n) = coords(i,n)
        bonds(i,1,n) = coords(i,n)
        bonds(i,1,4) = i
     end do
!
     do j = 1,atom_count
!
        do m = 1,3
           coord_2(m) = coords(j,m)
           store = j
        end do
!        
        call bondlength(coord_1,coord_2,bond_length)
!        
        n = int(bonds(i,1,4))

        hydrogen = index(atoms(n),'H')
        hydrogen2 = index(atoms(j),'H')
!
        if(hydrogen.ne.0.or.hydrogen2.ne.0) then
           if (bond_length.ge.0.4.and.bond_length.le.1.2) then
              storepoint = storepoint + 1
              n = 0
           else
              n = 0
           end if
        else
           if (bond_length.ge.0.4.and.bond_length.le.2.0) then
              storepoint = storepoint + 1
              n = 0
           else
              n = 0
           end if
        end if
     end do                                                ! JMcD i loop
  end do                                                   ! JMcD j loop
  write(200,'(A15,I3)') 'No. bonded atoms ', storepoint - 1
  n = 0
!    
!------------------------- Print information to log file -------------------------------------------
!
  l = 1
  k = 4
  m = 2
  n = 3
  write(200,'(20x)')
  write(200,'(A)')'----------- Array Information -----------'
  write(200,'(A)')' Atom Type          Cells               Labels            Coordinates (xyz)'
  write(200,'(A)')'--------------------------------------------------------------------------------------------'
  do i = 1, atom_count
     do j = 1, atom_count
        a = int(bonds(i,1,4))
        b = int(bonds(i,j,4))
        if (a.ne.0) then
           write(200,'(A15,3(I3,1X),f20.8,1X,A3,5X,3(f10.7,1X))') 'Host atom      ', i, l, k, bonds(i,1,4), atoms(a),& 
                bonds(i,1,l), bonds(i,1,m), bonds(i,1,n)
        else
           write(200,'(A15,3(I3,1X),f20.8,9X,3(f10.7,1X))') 'Host atom      ', i, l, k, bonds(i,1,4),& 
                bonds(i,1,l), bonds(i,1,m), bonds(i,1,n)
        end if
!        
        if (b.ne.0) then
           write(200,'(A15,3(I3,1X),f20.8,1X,A3,5X,3(f10.7,1X))') 'Bonded atom    ', i, j, k, bonds(i,j,4), atoms(b),& 
                bonds(i,j,l), bonds(i,j,m), bonds(i,j,n)
        else
           write(200,'(A15,3(I3,1X),f20.8,9X,3(f10.7,1X))') 'Bonded atom    ', i, j, k, bonds(i,j,4),& 
                bonds(i,j,l), bonds(i,j,m), bonds(i,j,n)
        end if
     end do
  end do
  l = 0
  k = 0
  m = 0
  n = 0
  
!
!-------------------------- Make a sorted list -----------------------------------------------------
!
  do i = 1, atom_count
     label = int(bonds(i,1,4))  
     print*,'Before sort: Atom number ',bonds(i,1,4),' Atom label ', atoms(label),' No. of bonded atoms ', bonds(i,1,5)
  end do
!
!- - - - - Sort 1 Num of bonded atoms - - - - -
!
  do j = atom_count-1, 1, -1
     swap = .FALSE.
     do i = 1, j
        if (bonds(i,1,5).gt.bonds(i + 1,1,5)) then
           do k = 1, 5
              temp(k) = bonds(i,1,k)
              bonds(i,1,k) = bonds(i + 1,1,k)
              bonds(i + 1,1,k) = temp(k)
           end do
           swap = .TRUE.
        end if
     end do
     if (.NOT. swap) exit
  end do
!
!- - - - - Sort 2 Hydrogen first - - - - -
!
  do j = atom_count-1, 1, -1
     swap = .FALSE.
     do i = 1, j
        ilabel = bonds(i,1,4)
        i2label = bonds(i+1,1,4)
        if (index(atoms(ilabel),'H').eq.0.and.index(atoms(i2label),'H').ne.0) then
           do k = 1, 5
              temp(k) = bonds(i,1,k)
              bonds(i,1,k) = bonds(i + 1,1,k)
              bonds(i + 1,1,k) = temp(k)
           end do
           swap = .TRUE.
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
!---------------------- Find the unique bonding combinations ---------------------------------------
!   
  i = 1
  j = 1
  hostatom = 1
  storepoint = 1
  do i = 1,atom_count
!        
     if(hostatom.ne.i) then
        write(200,'(A17,I3)') 'No. bonded atoms ',storepoint - 1
        bonds(hostatom,1,5) = storepoint - 1
        hostatom = hostatom + 1                        ! JMcD Count of atoms held constant while other atom is varied
        storepoint = 1                                 ! JMcD The array index to store the label of bonded atoms to the host atom in storepoint - 1 is the number of bonded atoms to the host atom
     end if
!
     write(200,'(20X)')
     write(200,'(A)')'Unique bonding pairs'
     write(200,'(A)')'------------ Next Host Atom ------------'
     write(200,'(A11,4X,A5,4X,A12)')'Bonding/Not','Atoms','Bond length'
     write(200,'(A)')'-------------------------------------'
!     
     do n = 1,3
        coord_1(n) = bonds(i,1,n)
     end do
!
     do j = i,atom_count
!
        do m = 1,3
           coord_2(m) = bonds(j,1,m)
           store = j
        end do
!        
        call bondlength(coord_1,coord_2,bond_length)
!        
        n = int(bonds(i,1,4))
        m = int(bonds(j,1,4))
        hydrogen = index(atoms(n),'H')
        hydrogen2 = index(atoms(m),'H')
!
        if(hydrogen.ne.0.or.hydrogen2.ne.0) then
           if (bond_length.ge.0.4.and.bond_length.le.1.2) then
              storepoint = storepoint + 1      
              write(200,'(A14,2(1XA3),f10.5)')'Bonding       ', atoms(n), atoms(m), bond_length ! JMcD Print to log file
              write(160,'(A1,1X,2(A3,1X))') 'B',atoms(n), atoms(m)                              ! JMcD Print to structure input
              bonds(hostatom,storepoint,4) = m
              m = 0
              do m = 1,3
                 bonds(hostatom,storepoint,m) = coord_2(m)
              end do
              n = 0
           else
              write(200,'(A14,2(1XA3),f10.5)')'Not bonding   ', atoms(n), atoms(m), bond_length ! JMcD log file
              n = 0
           end if

        else
           if (bond_length.ge.0.4.and.bond_length.le.2.0) then
              storepoint = storepoint + 1      
              write(200,'(A14,2(1XA3),f10.5)')'Bonding       ', atoms(n), atoms(m), bond_length ! JMcD Print to log file
              write(160,'(A1,1X,2(A3,1X))') 'B',atoms(n), atoms(m)                              ! JMcD Print to structure input
              bonds(hostatom,storepoint,4) = m
              do m = 1,3
                 bonds(hostatom,storepoint,m) = coord_2(m)
              end do
              n = 0
           else
              write(200,'(A14,2(1XA3),f10.5)')'Not bonding   ', atoms(n), atoms(m), bond_length ! JMcD log file
              n = 0
           end if
        end if
     end do                                                ! JMcD j loop
  end do                                                   ! JMcD i loop
  write(200,'(A15,I3)') 'No. bonded atoms ', storepoint - 1
  n = 0
  write(200,'(20X)')
!- - - - - - - - - Log - - - - - - - - - - - - - - 
  write(200,'(A)') 'Bonded list'
  write(200,'(13(A5,1X))') atoms 
!  do var = 1, atom_count  
!     write(200,'(13(I5,1X))') int(bonds(1,var,4)), int(bonds(2,var,4)), int(bonds(3,var,4)), &
!          int(bonds(4,var,4)), int(bonds(5,var,4)), int(bonds(6,var,4)) ,int(bonds(7,var,4)), &
!          int(bonds(8,var,4)), int(bonds(9,var,4)), int(bonds(10,var,4)), int(bonds(11,var,4)) &
!          , int(bonds(12,var,4)), int(bonds(13,var,4))      ! JMcD DIAGNOS
!  end do
  write(200,'(20X)')
!
!--------------- Re-zero bonds --------------------------------------------------
!
!  bonds = 0.0
!
!-------------------- Get all unique angles -------------------------------------
!
  i = 1
  j = 1
  hostatom = 1
  storepoint = 1
  do i = 1,atom_count
!        
     if(hostatom.ne.i) then
        bonds(hostatom,1,5) = storepoint - 1
        hostatom = hostatom + 1                        ! JMcD Count of atoms held constant while other atom is varied
        storepoint = 1                                 ! JMcD The array index to store the label of bonded atoms to the host atom in storepoint - 1 is the number of bonded atoms to the host atom
     end if
!     
     do n = 1,3
        coord_1(n) = bonds(i,1,n)
     end do
!
     do j = 1,atom_count
!
        do m = 1,3
           coord_2(m) = bonds(j,1,m)
           store = j
        end do
!
        call bondlength(coord_1,coord_2,bond_length)
!        
        n = int(bonds(i,1,4))
	!print*, 'n = ', n, ' bonds(i,1,4) = ', bonds(i,1,4)
        m = int(bonds(j,1,4))
        hydrogen = index(atoms(n),'H')
        hydrogen2 = index(atoms(m),'H')
!
        if(hydrogen.ne.0.or.hydrogen2.ne.0) then
           if (bond_length.ge.0.4.and.bond_length.le.1.2) then
              storepoint = storepoint + 1      
              bonds(hostatom,storepoint,4) = m
              m = 0
              do m = 1,3
                 bonds(hostatom,storepoint,m) = coord_2(m)
              end do
              n = 0
           else
              n = 0
           end if

        else
           if (bond_length.ge.0.4.and.bond_length.le.2.0) then
              storepoint = storepoint + 1
              bonds(hostatom,storepoint,4) = m
              m = 0
              do m = 1,3
                 bonds(hostatom,storepoint,m) = coord_2(m)
              end do
              n = 0
           else
              n = 0
           end if
        end if
     end do                                                ! JMcD j loop
  end do                                                   ! JMcD i loop
  n = 0
  print*, atoms
  write(200,'(20X)')
!
!-------------------- Find the unique angles ------------------------------  
!
  write(200,'(20X)')
  write(200,'(A)') 'Unique atom Angle trimers'
  i = 1
  j = 1
  hostatom = 1
  storepoint = 1
  do i = 1,atom_count
!        
     if(hostatom.ne.i) then
        hostatom = hostatom + 1                       
        storepoint = 1                                 
     end if
!
     write(200,'(20X)')
     write(200,'(A)')'Unique Angle trimers'
     write(200,'(A)')'------------ Next Host Atom ------------'
     write(200,'(A9,20X,A5,4X)')'Angle/Not','Atoms'
     write(200,'(A)')'----------------------------------------'
!     
     do n = 1,3
        coord_1(n) = bonds(i,1,n)
     end do
!
     do j = 1,atom_count
!
        do m = 1,3
           coord_2(m) = bonds(j,1,m)
           store = j
        end do
!        
        call bondlength(coord_1,coord_2,bond_length)
!        
        n = int(bonds(i,1,4))
        m = int(bonds(j,1,4))
        hydrogen = index(atoms(n),'H')
        hydrogen2 = index(atoms(m),'H')
!
        if(hydrogen2.ne.0) then
           cycle
        else if(hydrogen.ne.0.or.hydrogen2.ne.0) then
           if (bond_length.ge.0.4.and.bond_length.le.1.2) then
              storepoint = storepoint + 1      
              do k = 2, atom_count
                 if(bonds(j,k,4).eq.0) then
                    exit
                 end if
!                 
                 label = int(bonds(j,k,4))
                 if(n.eq.label) then
                    exit
                 end if
!                 
                 do l = 1,limit                    
                    angle = 'A ' // atoms(n) // atoms(m) // atoms(label)
                    a1 = index(angle_list(l),atoms(n))
                    a2 = index(angle_list(l),atoms(m))
                    a3 = index(angle_list(l),atoms(label))
                    angl = verify(angle,angle_list(l))
                    hydrogen = index(atoms(label),'H')                   
                    if(a1.ne.0.and.a2.ne.0.and.a3.ne.0) then
                       a123 = 1
                    end if
!                
                    if(angl.eq.0) then
                       reang = 1
                    end if
                 end do
!
                 if(a123.eq.0.and.reang.eq.1) then
                    write(200,'(38X,A104,3(1X,A3))') 'WARNING : module bonding : The code is unclear if this entry is a &
                         replicate it has been included anyway',atoms(n), atoms(m), atoms(label)
                    reang = 0
                 end if
!
                 if(reang.eq.0.and.a123.eq.0) then
                    angle_list_line = angle_list_line + 1
                    write(160,'(A1,1X,3(A3,1X))') 'A', atoms(n), atoms(m), atoms(label)
                    write(200,'(A24,1X,3(A3,1X),f10.5)') 'Valid angle term       ', atoms(n), atoms(m), atoms(label),bond_length
                    angle_list(angle_list_line) = trim(angle)
                 end if
                 reang = 0
                 a123 = 0
              end do
              n = 0
           else
              write(200,'(A24,3(1XA3))')'Not valid angle term   ', atoms(n), atoms(m), atoms(label) ! JMcD log file
              n = 0
           end if
        else
           if (bond_length.ge.0.4.and.bond_length.le.2.0) then
              storepoint = storepoint + 1      
              do k = 2, atom_count
                 if(bonds(j,k,4).eq.0) then
                    exit
                 end if
!                 
                 label = int(bonds(j,k,4))
                 if(n.eq.label) then
                    exit
                 end if
!                 
                 do l = 1,limit                    
                    angle = 'A ' // atoms(n) // atoms(m) // atoms(label)
                    a1 = index(angle_list(l),atoms(n))
                    a2 = index(angle_list(l),atoms(m))
                    a3 = index(angle_list(l),atoms(label))
                    angl = verify(angle,angle_list(l))
                    hydrogen = index(atoms(label),'H')                   
!
                    if(a1.ne.0.and.a2.ne.0.and.a3.ne.0) then
                       a123 = 1
                    end if
!
                    if(angl.eq.0) then
                       reang = 1
                    end if
                 end do
!
                 if(a123.eq.0.and.reang.eq.1) then
                    write(200,'(38X,A103,3(1X,A3))') 'WARNING : module bonding : The code is unclear if this entry is a  &
                         replicate it has been included anyway',atoms(n), atoms(m), atoms(label)
                    reang = 0
                 end if
!
                 if(reang.eq.0.and.a123.eq.0) then
                    angle_list_line = angle_list_line + 1
                    write(160,'(A1,1X,3(A3,1X))') 'A', atoms(n), atoms(m), atoms(label)
                    write(200,'(A24,1X,3(A3,1X),f10.5)') 'Valid angle term       ', atoms(n), atoms(m), atoms(label),bond_length
                    angle_list(angle_list_line) = 'A ' // atoms(n) // atoms(m) // atoms(label)
                 end if
                 reang = 0
                 a123 = 0
              end do
              n = 0
           else
              write(200,'(A24,3(1XA3))')'Not valid angle term   ', atoms(n), atoms(m), atoms(label) ! JMcD log file
              n = 0
           end if
        end if                                             ! JMcD hydrogen2
     end do                                                ! JMcD j loop
  end do                                                   ! JMcD i loop  
  n = 0
  angle_list_line = 0
  storepoint = 0
  hostatom = 0
  write(200,'(20X)')
!- - - - - - - - - Log - - - - - - - - - - - - - - 
  write(200,'(A)') 'Bonded list'
  write(200,'(13(A5,1X))') atoms 
!  do var = 1, atom_count
!     write(200,'(13(I5,1X))') int(bonds(1,var,4)), int(bonds(2,var,4)), int(bonds(3,var,4)), &
!          int(bonds(4,var,4)), int(bonds(5,var,4)), int(bonds(6,var,4)) ,int(bonds(7,var,4)), &
!          int(bonds(8,var,4)), int(bonds(9,var,4)), int(bonds(10,var,4)), int(bonds(11,var,4)) &
!          , int(bonds(12,var,4)), int(bonds(13,var,4))      ! JMcD DIAGNOS
!  end do
  write(200,'(20X)')
!  print*, atoms                                                                                                                               ! JMcD DIAGNOS
!  write(*,'(13(I5,1X))') int(bonds(1,1,4)), int(bonds(2,1,4)), int(bonds(3,1,4)), int(bonds(4,1,4)), int(bonds(5,1,4)), int(bonds(6,1,4)) &   ! JMcD DIAGNOS
!,int(bonds(7,1,4)), int(bonds(8,1,4)), int(bonds(9,1,4)), int(bonds(10,1,4)), int(bonds(11,1,4)), int(bonds(12,1,4)), int(bonds(13,1,4))      ! JMcD DIAGNOS
!
!------------------------- Blank the angle_list array --------------------------
!
  do i = 1, limit
     angle_list(i) = ''
  end do
!
!----------------------- Find the unique dihedrls ------------------------------
!
  write(200,'(20X)')
  write(200,'(A)') 'Unique atom Dihedral quadmers'
  i = 1
  j = 1
  hostatom = 1
  storepoint = 1
  reang = 0
  a1234 = 0
  do i = 1,atom_count

!        
     if(hostatom.ne.i) then
        hostatom = hostatom + 1                       
        storepoint = 1                                 
     end if
!
     write(200,'(20X)')
     write(200,'(A)')'Unique dihedrals quadmers'
     write(200,'(A)')'------------ Next Host Atom ------------------'
     write(200,'(A11,15X,A5,4X,A12)')'Bonding/Not','Atoms'
     write(200,'(A)')'----------------------------------------------'
!     
     do n = 1,3
        coord_1(n) = bonds(i,1,n)
     end do
!
     do j = 1,atom_count
!
        do m = 1,3
           coord_2(m) = bonds(j,1,m)
        end do
! 
        call bondlength(coord_1,coord_2,bond_length)
        bond_length1 = bond_length
!        
        n = int(bonds(i,1,4))
        m = int(bonds(j,1,4))
        if(i.eq.j) then
            cycle
        end if  
!
        do k = 1, atom_count
!           
           do x = 1,3
              coord_3(x) = bonds(k,1,x)
           end do
!           
           call bondlength(coord_2,coord_3,bond_length)
           o = int(bonds(k,1,4))
           hydrogen = index(atoms(n),'H')
           hydrogen2 = index(atoms(m),'H')
           hydrogen3 = index(atoms(o),'H')
!   
           if(hydrogen2.ne.0) then
              exit
           end if
!
           if(hydrogen3.ne.0) then
                cycle
           end if
           
           if(n.eq.o) then
               cycle
           end if
           
!
           if(hydrogen.ne.0) then
              if (bond_length1.ge.0.4.and.bond_length1.le.1.2.and.bond_length.ge.0.4.and.bond_length.le.2.0) then                   
!
                 storepoint = storepoint + 1      
                 do x = 2, atom_count
                    if(bonds(k,x,4).eq.0) then
                       exit
                    end if
!
                    label = int(bonds(k,x,4))
!
                    if(n.eq.label.or.m.eq.label) then
                       cycle
                    end if
!                 
                    do l = 1,limit                    
                       dihedrals = 'A ' // atoms(n) // atoms(m) // atoms(o) // atoms(label)
                       a1 = index(angle_list(l),atoms(n))
                       a2 = index(angle_list(l),atoms(m))
                       a3 = index(angle_list(l),atoms(o))
                       a4 = index(angle_list(l),atoms(label))
                       angl = verify(angle,angle_list(l))
                       hydrogen = index(atoms(label),'H')                   
!
                       if(a1.ne.0.and.a2.ne.0.and.a3.ne.0.and.a4.ne.0) then
                          a1234 = 1
                       end if
!
                       if(angl.eq.0) then
                          reang = 1
                       end if
                    end do
!
                    if(a1234.eq.0.and.reang.eq.1) then
                       write(200,'(38X,A103,4(1X,A3))') 'WARNING : module bonding : The code is unclear if this entry is a  &
                            replicate it has been included anyway',atoms(n), atoms(m), atoms(o), atoms(label)
                       reang = 0
                    end if
!
                    if(reang.eq.0.and.a1234.eq.0) then
                       angle_list_line = angle_list_line + 1
                       write(160,'(A1,1X,4(A3,1X))') 'D', atoms(n), atoms(m), atoms(o), atoms(label)
                       write(200,'(A29,1X,4(A3,1X))') 'Valid dihedral s term       ', atoms(n), atoms(m), atoms(o),atoms(label)
                       angle_list(angle_list_line) = 'D ' // atoms(n) // atoms(m) // atoms(o) // atoms(label)
                    end if
                    reang = 0
                    a1234 = 0
                 end do                 ! JMcD x loop
              else
                 write(200,'(A29,4(1XA3))')'Not valid dihedral s term   ', atoms(n), atoms(m), atoms(o), atoms(label) ! JMcD log file
              end if                    ! JMcD bond length
           else
              if(bond_length1.ge.0.4.and.bond_length1.le.2.0.and.bond_length.ge.0.4.and.bond_length.le.2.0) then
                 storepoint = storepoint + 1      
                 do x = 2, atom_count
                    if(bonds(k,x,4).eq.0) then
                       exit
                    end if
!
                    label = int(bonds(k,x,4))
                    if(n.eq.label.or.m.eq.label) then
                       exit
                    end if
!            
                    do l = 1,limit                    
                       dihedrals = 'A ' // atoms(n) // atoms(m) // atoms(o) // atoms(label)
                       a1 = index(angle_list(l),atoms(n))
                       a2 = index(angle_list(l),atoms(m))
                       a3 = index(angle_list(l),atoms(o))
                       a4 = index(angle_list(l),atoms(label))
                       angl = verify(angle,angle_list(l))
                       hydrogen = index(atoms(label),'H')                   
!
                       if(a1.ne.0.and.a2.ne.0.and.a3.ne.0.and.a4.ne.0) then
                          a1234 = 1
                       end if
!
                       if(angl.eq.0) then
                          reang = 1
                       end if
                    end do
!
                    if(a1234.eq.0.and.reang.eq.1) then
                       write(200,'(38X,A103,4(1X,A3))') 'WARNING : module bonding : The code is unclear if this entry is a  &
                           replicate it has been included anyway',atoms(n), atoms(m), atoms(o), atoms(label)
                       reang = 0
                    end if
!
                    if(reang.eq.0.and.a1234.eq.0) then
                       angle_list_line = angle_list_line + 1
                       write(160,'(A1,1X,4(A3,1X))') 'D', atoms(n), atoms(m), atoms(o), atoms(label)
                       write(200,'(A29,1X,4(A3,1X))') 'Valid dihedral s term       ', atoms(n), atoms(m), atoms(o),atoms(label)
                       angle_list(angle_list_line) = 'D ' // atoms(n) // atoms(m) // atoms(o) // atoms(label)
                    end if
                    reang = 0
                    a1234 = 0
                 end do
              else
                 write(200,'(A29,4(1XA3))')'Not valid dihedral s term   ', atoms(n), atoms(m), atoms(o), atoms(label) ! JMcD log file
              end if                   ! JMcD bond length
           end if                      ! JMcD first atom is hydrogen
        end do                         ! JMcD k loop
     end do                            ! JMcD j loop
  end do                               ! JMcD i loop
!- - - - - - - - - Log - - - - - - - - - - - - - - 
  write(200,'(A)') 'Bonded list'
  write(200,'(13(A5,1X))') atoms 
!  do var = 1, atom_count
!       write(200,'(13(I5,1X))') int(bonds(1,var,4)), int(bonds(2,var,4)), int(bonds(3,var,4)), &
!          int(bonds(4,var,4)), int(bonds(5,var,4)), int(bonds(6,var,4)) ,int(bonds(7,var,4)), &
!          int(bonds(8,var,4)), int(bonds(9,var,4)), int(bonds(10,var,4)), int(bonds(11,var,4)) &
!          , int(bonds(12,var,4)), int(bonds(13,var,4))
!  end do
  write(200,'(20X)')
!
!-------------------- Closing and deallocating memory ---------------------------
!
  write(160,'(A)')'END'
  deallocate(atoms, coords, bonds)
  deallocate(angle_list)
  close(50)
  close(160)
end subroutine bonding
!
end module module_bonding
