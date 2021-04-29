MODULE populate_class
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Filename: populate_class.f90                                             !!
!!                                                                           !!
!!  Dependant files:                                                         !!
!!                                                                           !!
!!  Author: Oliver Conway                             Start date: 28/04/2021 !!
!!                                                                           !!
!!  Purpose: Base class for matrix populator hierarchy.                      !!
!!                                                                           !!
!!  Revisions:                                                               !!
!!                                                                           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IMPLICIT NONE
! Type definition
TYPE,PUBLIC :: populate ! This will be the name we instantiate
! Instance variables.
CONTAINS
! Bound procedures
PROCEDURE,PUBLIC :: populate => populate_sub ! Populates matrix
END TYPE populate
! Restrict access to the actual procedure names
private :: populate_sub
! Now add methods
CONTAINS

subroutine populate_sub(this)
    class(populate),INTENT(IN) :: this
END subroutine populate_sub

END MODULE populate_class
