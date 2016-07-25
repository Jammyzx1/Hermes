!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Originally created by James McDonagh at the University of Manchester 2015, in the Popelier group      ! 
! Components of the module are acknowledged to Dr Tanja van Mourik University of St Andrews             ! 
! Licensed under Attribution-ShareAlike 2.5 Generic (CC BY-SA 2.5)                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Originally created by James.
! Version 1.1
! CHANGE LOG 
! Version 1  : Paramters for calculating structural parameters
! Version 1.1: Modular format and incorporation in the Hermes program. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module module_constants 

  implicit none

  real, parameter :: pi = 3.1415926536
  real, parameter :: e = 2.7182818285
  real, parameter :: degrees_to_radians = 0.0174532925
  Real, parameter :: radians_to_degrees = 57.2957795

  contains
    subroutine show_contants()
      Print*, "Constants used are "
      Write(*,'(A,f12.10)') "pi = ",pi
      Write(*,'(A,f12.10)') "e = ",e
      Write(*,'(A,f12.10)') "Degrees to radians = ",degrees_to_radians
      Write(*,'(A,f12.10)') "Radians to degrees = ",radians_to_degrees
    end subroutine show_contants
  end module module_constants
