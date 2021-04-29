MODULE matrix_class
  use region_class_1d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Filename: matrix_class.f90                                               !!
!!                                                                           !!
!!  Dependant files: egion_class_1d.f90                                      !!
!!                                                                           !!
!!  Author: Oliver Conway                             Start date: 28/04/2021 !!
!!                                                                           !!
!!  Purpose: Base class for matrix hierarchy.                                !!
!!                                                                           !!
!!  Revisions:                                                               !!
!!                                                                           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IMPLICIT NONE
! Type definition
TYPE,PUBLIC :: matrix ! This will be the name we instantiate
! Instance variables.
CONTAINS
! Bound procedures
PROCEDURE,PUBLIC :: solve => solve_fn ! Performs a matrix equation solve
END TYPE matrix
! Restrict access to the actual procedure names
private :: solve_fn
! Now add methods
CONTAINS

FUNCTION solve_fn(this) result(solution)
  class(matrix),INTENT(INOUT) :: this
  real(dp),allocatable,dimension(:) :: solution
END FUNCTION solve_fn

END MODULE matrix_class
