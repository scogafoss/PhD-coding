MODULE populate_class
    use compressed_matrix_class
    use tridiagonal_matrix_class
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
! PROCEDURE,PUBLIC :: populate_matrix => populate_matrix_sub ! Populates matrix
END TYPE populate
! Restrict access to the actual procedure names
! private :: populate_matrix_sub
! Now add methods
CONTAINS

! SUBROUTINE populate_matrix_sub(this,in_matrix,regions,group)
!     !
!     ! Subroutine to perform discretisation.
!     !
!     implicit none
!     class(populate), intent(inout) :: this ! Nuclear_matrix object
!     class(matrix),INTENT(INOUT) :: in_matrix
!     type(region_1d), intent(in), dimension(:) :: regions ! Region object
!     INTEGER,INTENT(IN) :: group ! Which group's data is required
! END subroutine populate_matrix_sub

END MODULE populate_class
