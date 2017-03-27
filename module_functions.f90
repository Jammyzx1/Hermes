!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Originally created by James McDonagh at the University of Manchester 2015, in the Popelier group      ! 
! Components of the module are acknowledged to Dr Tanja van Mourik University of St Andrews             ! 
! Licensed under MIT License                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Originally created by James.
! Version 1.1
! CHANGE LOG 
! Version 1  : Functions for calculating strutcural parameters
! Version 1.1: Modular format and incorporation in the Hermes program.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module module_functions 

use module_constants

  implicit none
  
contains
  
  function bond_vector(v1,v2)
    !Dummy arguements
    real, dimension(3) :: v1,v2
    !Function argument
    real, dimension(3) :: bond_vector
    !local
    integer :: i
    
    do i=1,3
       bond_vector(i)=v2(i)-v1(i)
    end do

  end function bond_vector


!************************************************************************

  function dist(x1, y1, z1, x2, y2, z2)	! Calculates the distance between two coordinates given as separate values
    implicit none                       ! Function originally created by Rob Coates Univeristy of Manchester 2015
	
    real, dimension(3) :: vector
    real :: distance2, dist, x1, y1, z1, x2, y2, z2
    integer :: i
  
    i = 0
    distance2 = 0
    dist = 0
    vector = 0
	
    vector(1) = x1 - x2
    vector(2) = y1 - y2
    vector(3) = z1 - z2
	
    do i = 1,3,1
       distance2 = distance2 + (vector(i)**2.0)
    end do
	
    dist = sqrt(distance2)

  end function dist


!************************************************************************

  Function norm(v)

    !dummy arguments
    real, dimension (3) :: v
    !function result
    real :: norm

    norm=sqrt(v(1)*v(1) + v(2)*v(2) + v(3)*v(3))
!    Print*, "norm is = ", norm                ! JMcD DIAGNOSTIC PRINT

  end Function norm


!************************************************************************

  function normalise(vector,normal)
    !Dummy arguments 
    real, dimension(3) :: vector
    real :: normal
    !fucntion argument
    real, dimension(3) ::normalise
    !local
    integer :: inc

    do inc=1,3,1
       normalise(inc)=vector(inc)/normal
    end do

  end function normalise


!************************************************************************

  function vector_cross(bvec1,bvec2)
    !Dummy arguments
    real, dimension(3) :: bvec1, bvec2
    !Function argument
    real, dimension(3) :: vector_cross

    vector_cross(1)=(bvec1(2)*bvec2(3)-bvec1(3)*bvec2(2))
    vector_cross(2)=(bvec1(3)*bvec2(1)-bvec1(1)*bvec2(3))
    vector_cross(3)=(bvec1(1)*bvec2(2)-bvec1(2)*bvec2(1))

  end function vector_cross


!************************************************************************

  function vector_dot(v1,v2)
    !Dummy arguments
    real, dimension(3) :: v1, v2
    !function argument
    real :: vector_dot
   
    vector_dot=v1(1)*v2(1)+v1(2)*v2(2)+v1(3)*v2(3)
    
  end function vector_dot


!************************************************************************

  function perdist(x1, y1, z1, x2, y2, z2, boxlen)	    ! Calculates the distance between two coordinates in periodic Boundary conditions
    implicit none                                       ! Function originally created by Rob Coates University of Manchester 2015
	
    real, dimension(3) :: vector
    real :: distance2, distance, x1, y1, z1, x2, y2, z2, boxlen, perdist, halfboxlen
    integer :: i
	
    i = 0
    distance2 = 0
    distance = 0
    vector = 0
  
    vector(1) = x1 - x2
    vector(2) = y1 - y2
    vector(3) = z1 - z2
    halfboxlen = boxlen*0.5
  
    do i = 1,3,1					! R.C. takes into account periodic boundry conditions
       if (vector(i).gt.halfboxlen) then
           vector(i) = vector(i) - boxlen
!       else if(vector(i).lt.-halfboxlen) then
!           vector(i) = vector(i) + boxlen
       else 
          cycle
       end if
    end do
	
    do i = 1,3,1
       distance2 = distance2 + (vector(i)**2.0)
    end do
	
    distance = sqrt(distance2)
    perdist = distance

  end function perdist
  
  
!************************************************************************ 

function distancemicro(x1, y1, z1, x2, y2, z2, bl)	! Calculates the distance between two coordinates
  implicit none
	
  real, dimension(3) :: vector
  real :: distance2, distancemicro, x1, y1, z1, x2, y2, z2, bl
  integer :: i
	
  i = 0
  distance2 = 0
  distancemicro = 0
  vector = 0
  
  vector(1) = x1 - x2
  vector(2) = y1 - y2
  vector(3) = z1 - z2
  
  do i = 1,3,1					!takes into account periodic boundry conditions
     if (vector(i).GT.(bl*0.5)) then
     vector(i) = vector(i) - bl
  else 
     cycle
  end if
end do
	
do i = 1,3,1
   distance2 = distance2 + (vector(i)**2.0)
end do
	
distancemicro = sqrt(distance2)

end function distancemicro


!************************************************************************

function distmicro(x1, y1, z1, x2, y2, z2)	! Calculates the distance between two coordinates
  implicit none
	
  real, dimension(3) :: vector
  real :: distance2, distmicro, x1, y1, z1, x2, y2, z2
  integer :: i
  
  i = 0
  distance2 = 0
  distmicro = 0
  vector = 0
	
  vector(1) = x1 - x2
  vector(2) = y1 - y2
  vector(3) = z1 - z2
	
  do i = 1,3,1
     distance2 = distance2 + (vector(i)**2.0)
  end do
	
  distmicro = sqrt(distance2)

end function distmicro


!************************************************************************

end module module_functions
