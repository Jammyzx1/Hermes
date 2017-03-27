!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Originally created by James McDonagh at the University of Manchester 2015, in the Popelier group        ! 
! Components of the dihedral subroutine are acknowledged to Dr Tanja van Mourik University of St Andrews  ! 
! Licensed under MIT License                                        !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Originally created by James.
! Version 1.1
! CHANGE LOG 
! Version 1  : Subroutines for calculating strutcural parameters
! Version 1.1: Modular format and incorporation in the Hermes program. Additional routines for jmol images
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module module_subroutines 

  use module_functions 
  use module_constants 

  implicit none
contains

  subroutine bondlength(atom_co_1,atom_co_2,bond_length)
    real, dimension(3), intent(in) :: atom_co_1, atom_co_2
    real, dimension(3) :: bond_vec
    real, intent(inout) :: bond_length
    
    bond_vec = bond_vector(atom_co_1,atom_co_2)
    bond_length = norm(bond_vec)
    
  end subroutine bondlength


!**************************************************************************************************************

  subroutine bondangle(atom_co_1,atom_co_2,atom_co_3,angle,angle_deg)
    
    real, dimension(3), intent(in) :: atom_co_1, atom_co_2, atom_co_3
    real, dimension(3) :: vec1, vec2
    real :: normvec1, normvec2, norms, dotvec1vec2, argument 
    real,intent(inout) :: angle, angle_deg
    
    vec1 = bond_vector(atom_co_1, atom_co_2)
    vec2 = bond_vector(atom_co_3, atom_co_2)
      
    normvec1 = norm(vec1)
    normvec2 = norm(vec2)
    norms = normvec1 * normvec2

    dotvec1vec2 = vector_dot(vec1,vec2)
!    print*,'Dot product',dotvec1vec2               ! JMcD DIAGNOSTIC PRINT
    
    argument = dotvec1vec2 / norms
!   print*, argument                                ! JMcD DIAGNOSTIC PRINT
    
    angle = acos(argument)
    angle_deg = angle * radians_to_degrees
    
  end subroutine bondangle


!**************************************************************************************************************

  subroutine dihedral(atom_co_1,atom_co_2,atom_co_3,atom_co_4,dihedangle,dihedangle_deg)
    
    integer :: i
    real, dimension(3), intent(in) :: atom_co_1, atom_co_2, atom_co_3, atom_co_4
    real, dimension(3) :: b1, b2, b3, cross1, cross2, m
    real :: norm_b1, norm_b2, norm_b3, y, x
    real, intent(inout) :: dihedangle, dihedangle_deg
    Character (len=1) :: debug="F" ! Note debug is define on = T or off = F here
      
    If(debug=="T")then
       Print*, "In subroutine dihedral debug is turned on output will be verbose"
       Print*, "atom 1 coordinates are ", atom_co_1
       Print*, "atom 2 coordinates are ", atom_co_2
       Print*, "atom 3 coordinates are ", atom_co_3
       Print*, "atom 4 coordinates are ", atom_co_4
    end If
      
    b1=bond_vector(atom_co_1, atom_co_2)
    b2=bond_vector(atom_co_2, atom_co_3)
    b3=bond_vector(atom_co_3, atom_co_4)
     
    norm_b1=norm(b1)
    norm_b2=norm(b2)
    norm_b3=norm(b3)
    
    b1=normalise(b1,norm_b1)
    b2=normalise(b2,norm_b2)
    b3=normalise(b3,norm_b3)
    
    If(debug=="T")then
       Print*, "normalised bond vectors are b1 ", b1
       Print*, "normalised bond vectors are b2 ", b2
       Print*, "normalised bond vectors are b3 ", b3
    end If

    cross1=vector_cross(b1,b2)
    cross2=vector_cross(b2,b3)
      
    If(debug=="T")then
       Print*, "cross product vector 1 ", cross1
       Print*, "cross product vector 2 ", cross2
    end If

    m=vector_cross(cross1,b2)
    
    y=vector_dot(m,cross2)
    x=vector_dot(cross1,cross2)
      
    dihedangle=-atan2(y,x)
    dihedangle_deg=dihedangle*radians_to_degrees

  end subroutine dihedral


!**************************************************************************************************************
  
  subroutine StripSpaces(string)                      ! JMcD Acknowldge Jauch 27 November 2014 http ://stackoverflow.com/questions/27179549/removing-whitespace-in-string 
    character(len=*) :: string
    integer :: stringLen 
    integer :: last, actual

    stringLen = len (string)
    last = 1
    actual = 1

    do while (actual.lt.stringLen)
        if (string(last:last) == ' ') then
            actual = actual + 1
            string(last:last) = string(actual:actual)
            string(actual:actual) = ' '
        else
            last = last + 1
            if (actual.lt.last) &
                actual = last
        endif
    end do

  end subroutine StripSpaces


!**************************************************************************************************************

  subroutine sort_asend(array,toend,sorted)              ! JMcD Acknowledge Rossetta code.org Sorting algorithms/Bubble sort  http://rosettacode.org/wiki/Sorting_algorithms/Bubble_sort#Fortran
    ! Dummy arguments
    integer,intent(in) :: toend
    real, intent(in) :: array(1:toend)
    real, intent(inout) :: sorted(1:toend)
    ! Variables
    real :: temp
    logical :: swap
    integer :: i, j
    
    temp = 0
    sorted = 0
    i = 0
    j = 0
    sorted = array
    
    do j = toend, 1, -1
       swap = .FALSE.
       do i = 1, j
          if (sorted(i).gt.sorted(i + 1)) then
             temp = sorted(i)
             sorted(i) = sorted(i + 1)
             sorted(i + 1) = temp
             swap = .TRUE.
          end if
       end do
       if (.not. swap) exit
    end do
    
  end subroutine sort_asend


!**************************************************************************************************************

  subroutine create_norm_xyz(atom_counts,fil_nam,write_norm_xyz,title,ichar_title)

    character (len=100), intent(inout) :: fil_nam, write_norm_xyz
    integer, intent(in) :: atom_counts
    ! local variables
    character (len=5) :: atom_label, atom_label1
    character (len=10) :: atom_xyz_1, atom_xyz_2, atom_xyz_3
    integer ::  ierr, line_number, start_count, end_count, stat, ichars, location, location1
    integer :: geom, label_length, i, j, stats,l, structures, atom_number, ichar_title, temp_count
    character (len=11) :: search = "-----------" 
    character (len=256) :: line_from_50, title
    character (len=20), dimension(:,:),allocatable :: input_coords
    real, dimension(3) :: atom_coord_1, atom_coord_2, norm_coord

! - - - - - - - - - - - Open files - - - - - - - - - - -
    
    label_length = len_trim(write_norm_xyz)
    atom_label = write_norm_xyz
    allocate(input_coords(4,atom_counts))
    open (status="old", unit=50, action="read",access="sequential", file="labelled.xyz", iostat=ierr)                       ! unit 50 labelled.xyz
    if ( ierr.ne.0 ) then 
       write(*,'(A)') 'Unable to open labelled.xyz for creating normalised  xyz, please check file'
    else  
       write(*,'(A)')' Labelled.xyz has been opened successfully.'
    end if

    stat = 0
    structures = 0
    do while (stat.eq.0)
       read(50,'(A)',iostat=stat) line_from_50
       if (stat.gt.0) then
          print*,'ERROR - something went wrong when reading labelled.xyz.'
          stop
       else if (stat.lt.0) then
          print*,'WARNING - end of labelled.xyz file reached while counting structures ensure &
               END ifthe last line of the file for labelled.xyz.'
          stop
       else
          ichars=len_trim(line_from_50)
       end if

       if (line_from_50.eq.'END') then
          print*,'All geometries cycled'
          exit
       end if
       
       location = index(line_from_50(1:ichars),search)

       if(location.ne.0) then
          structures = structures + 1
       end if
    end do
    rewind(50)
 
    open (status="replace", unit=120, action="write",access="sequential", file="NORMAL.xyz", iostat=ierr)                         ! unit 120 NORMAL.xyz
    if ( ierr.ne.0 ) then 
       write(*,'(A)') 'Unable to open NORMAL.xyz for creating normalised  xyz, please check file'
    else  
       write(*,'(A)')' NORMAL.xyz has been opened successfully, now creating normalised xyz file'
    end if

    open (status="replace", unit=140, action="write",access="sequential", file="MULTI.xyz", iostat=ierr)                         ! unit 140 MULTI.xyz
    if ( ierr.ne.0 ) then
       write(*,'(A)') 'Unable to open MULTI.xyz for creating normalised  xyz, please check file'
    else
       write(*,'(A)')' MULTI.xyz has been opened successfully, now creating normalised xyz file'
    end if

    open (status="replace", unit=145, action="write",access="sequential", file="NORM-MULTI.xyz", iostat=ierr)                    ! unit 144 NORM-MULTI.xyz
    if ( ierr.ne.0 ) then
       write(*,'(A)') 'Unable to open NORM-MULTI.xyz for creating normalised  xyz, please check file'
    else
       write(*,'(A)')' NORM-MULTI.xyz has been opened successfully, now creating normalised xyz file'
    end if
! - - - - - - - - - - - Initalisation - - - - - - - - - -

    line_number = 0
    geom = 0
    start_count = 0
    end_count = 0
    stats = 0
    stat = 0
    do i = 1, 3
       atom_coord_1(i) = 0.0
       atom_coord_2(i) = 0.0
       norm_coord(i) = 0.0
    end do
! - - - - - - - - - - - Read and store/write output - - -
    atom_number = atom_counts * structures
    write(120,'(I10)') atom_number
    write(120,'(A)'), title(1:ichar_title)
    
    do while (stat.eq.0)

       line_number = line_number + 1
       read(50,'(A)',iostat=stat) line_from_50
       
       if (stat.gt.0) then
          print*,'ERROR - something went wrong when reading labelled.xyz for NORMAL.xyz.'
          stop
       else if (stat.lt.0) then
          print*,'WARNING - end of labelled.xyz file reached ensure END is the last line of the file for labelled.xyz.'
          stop
       else
          ichars=len_trim(line_from_50)
       end if
       
       if (line_from_50(1:ichars).eq.'END') then
          print*,'All geometries cycled'
          exit
       end if
       
       location = index(line_from_50(1:ichars),search)

       if(location.ne.0) then
          start_count = line_number + 1
          end_count = start_count + atom_counts
          geom = geom + 1
          write(140,'(I6)') atom_counts
          write(140,'(A17,I6)') 'Snapshot number ', geom
          write(145,'(I6)') atom_counts
          write(145,'(A17,I6)') 'Snapshot number ', geom
          temp_count = start_count
          do i = 1, atom_counts
             temp_count = temp_count + 1
             read(50,'(A)',iostat=stat) line_from_50
             location1 = index(line_from_50(1:ichars),atom_label(1:label_length))
             if(location1.ne.0) then

                read(line_from_50,*) atom_label, atom_xyz_1, atom_xyz_2, atom_xyz_3
             
                if (index(atom_xyz_1,'E').ne.0) then
                   read(atom_xyz_1,'(e15.8)') atom_coord_1(1)
                else
                   read(atom_xyz_1,'(f15.8)') atom_coord_1(1)
                end if

                if (index(atom_xyz_2,'E').ne.0) then
                   read(atom_xyz_2,'(e15.8)') atom_coord_1(2) 
                else
                   read(atom_xyz_2,'(f15.8)') atom_coord_1(2)
                end if

                if (index(atom_xyz_3,'E').ne.0) then
                   read(atom_xyz_3,'(e15.8)') atom_coord_1(3)
                else
                   read(atom_xyz_3,'(f15.8)') atom_coord_1(3)
                end if
             end if
          end do
          
          do i = 1, atom_counts
             backspace(50)
             temp_count = temp_count - 1
          end do
          if(start_count.ne.temp_count) then
             print*, 'Subroutine create_norm_xyz : ERROR : in locating atom coordintes to normalise against, please use &
                  MULTI.xyz as your visualisation file'
          end if
       end if
       
       if(line_number.ge.start_count.and.line_number.le.end_count - 1) then

          l = l + 1
          if(l.gt.atom_counts) then
             print*, 'ERROR - array assignment out of range when making NORMAL.xyz'
          end if

          read(line_from_50,*) input_coords(1,l), input_coords(2,l), input_coords(3,l), input_coords(4,l)

          read(input_coords(1,l),'(A)') atom_label1
          read(input_coords(2,l),'(f15.8)') atom_coord_2(1)
          read(input_coords(3,l),'(f15.8)') atom_coord_2(2)
          read(input_coords(4,l),'(f15.8)') atom_coord_2(3)
                
          do j = 1,3
             norm_coord(j) = atom_coord_2(j) - atom_coord_1(j)
          end do

          write(120,'(A3,3x,3(f15.8,3X))')atom_label1(:1), norm_coord(1), norm_coord(2), norm_coord(3)
          write(140,'(A3,3x,3(f15.8,3X))')atom_label1(:1), atom_coord_2(1), atom_coord_2(2), atom_coord_2(3)
          write(145,'(A3,3x,3(f15.8,3X))')atom_label1(:1), norm_coord(1), norm_coord(2), norm_coord(3)
          
          do j = 1, 3
             norm_coord(j) = 0.0
             atom_coord_2(j) = 0.0
          end do
       end if
       if(line_number.ge.end_count - 1) then
          l = 0
       end if
    end do
    close(50)
    close(120)
    close(140)
    deallocate(input_coords)
  end subroutine create_norm_xyz


!**************************************************************************************************************

  subroutine rdfcreate(tstep, bl, total, shell1, nx, ny, nz, solvent_atom_type, solutes_atom_type)     ! Subrountine originally created by Rob Coates Univeristy of Manchester 2015, Adapted by J. L. McDonagh
                                                                                                       ! R.C. calculates a rdf from an atom to surrounding water oxygen atoms
    ! Dummy variables
    integer, intent(in) :: tstep
    integer, intent(inout) :: total
    real, intent(in) :: bl
    real, intent(inout) :: shell1
    real, dimension(tstep), intent(inout) :: nx, ny, nz
    character (len=6), intent(in) :: solvent_atom_type, solutes_atom_type
    
    ! Local variables
    integer :: ier, i, counter, atm, stepcount, k, ncount, icharlensolv, icharlensolu
    real :: rho, dr, distance, rmax, dis, perdis
    real, dimension(:), allocatable :: x1, y1, z1
    real, dimension(:, :), allocatable :: r
    character (len=256) :: lrd
    character (len=17) :: dummy0, dummy1, dummy2, dummy3   
!    
! - - - - - - - - - - - - - - - - - - - - - Initialisation - - - - - - - - - - - - - - - - - - - - - 
!    
    ncount = 0
    k = 0
    stepcount = 0
    rmax = 0
    dr = 0
    rho = 0
    i = 0
    ier = 0
    counter = 0
    atm = 0
    shell1 = 0.0
    k = 12*bl ! divides box area by k into shells
    k = int(k)
    icharlensolv = len(trim(adjustl(solvent_atom_type)))
    icharlensolu = len(trim(adjustl(solutes_atom_type)))
!    
    print*,',' , solvent_atom_type, ',' , solutes_atom_type, ','
    allocate(x1(total), y1(total), z1(total))	  ! allocates an array for storing water O coordinates at each frame

! - - - - - - - - - - - - - - - - - - - - - Radial distribution setup - - - - - - - - - - - - - - - - - -

    allocate(r(k, 4))				  ! a 2D array containing distances, probability density, 1st and second derivatives of prob density
    rmax = sqrt(3.0*((bl*0.5)**2.0))		  ! the max distance an atom can be from another in the box 		
    
    dr = rmax / k		                  ! The shell distances
  
    do i = 1, k					  ! Fill array with shell distances
       r(i, 1) = dr*(i-1)
       r(i, 2) = 0.0
       r(i, 3) = 0.0
       r(i, 4) = 0.0
    end do

    rho = ((total)*tstep) / bl**3	          ! Mean atom density in the box across all timesteps
  
  
! - - - - - - - - - - - Calculate the distance for each atom of type - - - - - - - - - - - - - - - - - -

    print*, '|*********** CREATING * RDF *************************************************|'
    print*, '|****************************************************************************|'
    do while (ier.eq.0)
       
       read(51, "(A)", iostat=ier) lrd
!       print*, ier
!       print*,lrd
     
       if (index(lrd, solvent_atom_type(1:icharlensolv)).gt.0) then! Store coordinate data for each water O in a frame in arrays, 
          counter = counter + 1		          ! the atom label searched for may need changing depending on model used	
          read(lrd,*) dummy0,x1(counter), y1(counter), z1(counter)
!          print*,counter
!          read(dummy1, "(ES17.10)") x1(counter)					
!          read(dummy2, "(ES17.10)") y1(counter)
!          read(dummy3, "(ES17.10)") z1(counter)		
       else if (index(lrd, solutes_atom_type(1:icharlensolu)).gt.0) then
          ncount = ncount + 1		          ! counts through the atoms, one for each frame
          read(lrd,*) dummy0, nx(ncount), ny(ncount), nz(ncount)
!          read(dummy1, "(ES17.10)") nx(ncount) ! store coordinate data for use in rdf calculation and for use in later cutting out other water molecules		
!          read(dummy2, "(ES17.10)") ny(ncount)
!          read(dummy3, "(ES17.10)") nz(ncount)
       else if (index(lrd, 'Timestep').gt.0.and.stepcount.eq.0) then	! This only applies to the first line of the file
          stepcount = stepcount + 1				        ! begins counting of frames
          counter = 0
       else if (index(lrd, 'Timestep').gt.0.and.stepcount.gt.0) then	! Occurs at the end of each step
     
          do counter = 1, total			                        ! Calculates the distance from current target atom to each Oxygen atom in the current frame, using created function
             dis = perdist(x1(counter), y1(counter), z1(counter), nx(ncount), ny(ncount), nz(ncount), bl)
  
             do i = 2, k
                                                                !************|
                                                                ! Counts the |
                if (dis.ge.r(i-1,1).and.dis.lt.r(i,1)) then     !   number   |
                  r(i-1,2) = r(i-1,2) + 1                      	!     of     |
                else						                    !   Oxygen   |
                   cycle					                    !  atoms in  |
                end if						                    ! each shell |
                    					                        !************|   
           end do
          end do
    
          counter = 0	                                                 ! reset oxygen counter for next frame
    
          cycle
       else
          cycle
       end if
    end do
		
    do i = 1, k
       r(i,2) = r(i,2) / (rho*4*pi*((r(i,1)+dr)**2)*dr)		         ! Converts the number of atoms in each shell to a probability
    end do
	
    do i = 2, (k-1)
       r(i,3) = r(i+1,2) - r(i-1,2) 					 ! differential
    end do
	
    do i = 3, k
       r(i,4) = r(i,3) - r(i-1,3) 					 ! 2nd differential
    end do
	
    do i = 3, k
       if (abs(r(i, 3)).lt.0.08.and.r(i-1, 3).LT.0.and.r(i, 4).GT.0) then	! this tries to find the first minimum on the graph 
          shell1 = r(i, 1)						        ! May need tweaking to get correct, so check by plotting the rdf
          print*, 'First solvation shell: ', r(i, 1), ' Angstroms for atom ', solutes_atom_type  ! I've tried to make it less sensitive to small fluctuations
          exit								
       end if
    end do
	
! - - - - - - - - - - - Write RDF data to file - - - - - - - - - - - - - - - - - - - - - - - - - -
    
    write(52, "(A5,12X,A5)") solutes_atom_type, solvent_atom_type
    write(52, "(A)") 'Distance    Prob. Density    Derivative    2nd Derivative'													
	
    do i = 1, k									
       write(52, "(ES17.10, 2X, ES17.10, 2X, ES17.10, 2X, ES17.10)") r(i, 1), r(i, 2), r(i, 3), r(i, 4)	
    end do
														
    print*, '|*********** RDF * PRODUCED *************************************************|'
    print*, '|****************************************************************************|'
	
	
! - - - - - - - - - - - Clean up - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    
    deallocate(x1, y1, z1)
	
  end subroutine rdfcreate

  
!**************************************************************************************************************

  subroutine rdfcreatemicro(tstep, bl, total, shell1, nx, ny, nz, solvent_atom_type, solutes_atom_type) ! Subrountine created by Rob Coates Univeristy of Manchester 2015
    integer :: ier, i, counter, atm, total, stepcount, k, ncount, tstep
    real :: bl, rho, dr, distance, rmax, dis, shell1
    real, dimension(tstep) :: nx, ny, nz
    real, dimension(:), allocatable :: x1, y1, z1
    real, dimension(:, :), allocatable :: r
    real, parameter :: Pi = 3.1415927
    character (len=256) :: lrd
    character (len=6) :: solvent_atom_type, solutes_atom_type
	
! - - - - - - - - - - - - - - - - - - - - - I N I T I A L I S A T I O N - - - - - - - - - - - - - - - - - - - - - 
	
    ncount = 0
    k = 0
    stepcount = 0
    rmax = 0
    dr = 0
    rho = 0
    i = 0
    ier = 0
    counter = 0
    atm = 0
    k = 12*int(bl) !divides box area by k into shells
  
! - - - - - - - - - - - - - - - - - - - - - C O U N T - A T O M - T Y P E - - - - - - - - - - - - - - - - - - - - 
	
    allocate(x1(total), y1(total), z1(total))	!allocates an array for storing water O coordinates at each frame

! - - - - - - - - - - - - - - - - - - - - - R A D I A L - D I S - S E T U P - - - - - - - - - - - - - - - - - -

    allocate(r(k, 4))				! a 2D array containing distances, probability density, 1st and second derivatives of prob density
    rmax = sqrt(3.0*((bl*0.5)**2.0))		! the max distance an atom can be from another in the box 		
    
    dr = rmax / k					! The shell distances
  
    do i = 1, k					! Fill array with shell distances
       r(i, 1) = dr*(i-1)
       r(i, 2) = 0.0
       r(i, 3) = 0.0
       r(i, 4) = 0.0
    end do

    rho = ((total)*tstep) / bl**3	! Mean atom density in the box across all timesteps
  
  
! - - - - - - - - - - - C A L C U L A T E - T H E - D I S T A N C E S - F O R - E A C H - A T O M - O F - T Y P E - - - - - - -

    print*, '|*********** CREATING * RDF *************************************************|'
    print*, '|****************************************************************************|'
    do while (ier.EQ.0)
     
       read(51, "(A)", iostat=ier) lrd
     
       if (index(lrd, solvent_atom_type).GT.0) then	! Store coordinate data for each water O in a frame in arrays, 
          counter = counter + 1		! the atom label searched for may need changing depending on model used		
        
          read(lrd(7:24), "(ES17.10)") x1(counter)					
          
          read(lrd(26:43), "(ES17.10)") y1(counter)
        
          read(lrd(45:62), "(ES17.10)") z1(counter)
			
       else if (index(lrd, solutes_atom_type).GT.0) then	!change NH3 to whichever atom is chosen as target, here Nitrogen was used
          ncount = ncount + 1		!counts through the nitrogen atoms, one for each frame
        
          read(lrd(7:24), "(ES17.10)") nx(ncount)		!store coordinate data for use in rdf calculation		
        !and for use in later cutting out other water molecules
          read(lrd(26:43), "(ES17.10)") ny(ncount)
        
          read(lrd(45:62), "(ES17.10)") nz(ncount)
        
        
       else if (index(lrd, 'Timestep').GT.0.AND.stepcount.EQ.0) then	! This only applies to the first line of the file
          stepcount = stepcount + 1				! begins counting of frames
          counter = 0
       else if (index(lrd, 'Timestep').GT.0.AND.stepcount.GT.0) then	! Occurs at the end of each step
        
          do counter = 1, total			! Calculates the distance from current target atom 
									! to each Oxygen atom in the current frame, using created function
             dis = distancemicro(x1(counter), y1(counter), z1(counter), nx(ncount), ny(ncount), nz(ncount), bl)
           
             do i = 1, k								       !************|
                                                               ! Counts the |
                if (dis.GE.r(i,1).AND.dis.LT.r(i+1,1)) then    !   number   |
                   r(i,2) = r(i,2) + 1                         !     of     |
                else								           !   Oxygen   |
                   cycle							           !  atoms in  |
                end if								           ! each shell |
             end do								               !************|
           
          end do
    
          counter = 0	!reset oxygen counter for next frame
    
          cycle
       else
          cycle
       end if
    end do
		
    do i = 1, k
       r(i,2) = r(i,2) / (rho*4*Pi*((r(i,1)+dr)**2)*dr)		! Converts the number of atoms in each shell to a probability
    end do
	
    do i = 2, (k-1)
       r(i,3) = r(i+1,2) - r(i-1,2) 					! differential
    end do
	
    do i = 3, k
       r(i,4) = r(i,3) - r(i-1,3) 					! 2nd differential
    end do
	
    do i = 3, k
       if (abs(r(i, 3)).LT.0.08.AND.r(i-1, 3).LT.0.AND.r(i, 4).GT.0) then	!this tries to find the first minimum on the graph 
          shell1 = r(i, 1)						!May need tweaking to get correct, so check by plotting the rdf
          print*, 'First solvation shell: ', r(i, 1), ' Angstroms' 	!I've tried to make it less sensitive to small fluctuations
          exit								
       end if
    end do
	
! - - - - - - - - - - - W R I T E - R D F - D A T A - T O - F I L E - - - - - - - - - - - - - - - - - - - - - - - - - -
	
    write(52, "(A)") "N_3			OW"								
    write(52, "(A)") 'Distance    Prob. Density    Derivative    2nd Derivative'													
	
    do i = 1, k									
       write(52, "(ES17.10, 2X, ES17.10, 2X, ES17.10, 2X, ES17.10)") r(i, 1), r(i, 2), r(i, 3), r(i, 4)	
    end do
														
    print*, '|*********** RDF * PRODUCED *************************************************|'
    print*, '|****************************************************************************|'
	
	
! - - - - - - - - - - - C L E A N - U P - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    deallocate(x1, y1, z1)
	
  end subroutine rdfcreatemicro


!**************************************************************************************************************

end module module_subroutines
  
  
  
