MODULE tridiagonal_matrix_class
  use matrix_class
  use region_class_1d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Filename: matrix_class.f90                                               !!
!!                                                                           !!
!!  Dependant files: precision_set.f90, line_class.f90, material_class.f90,  !!
!!    region_class_1d.f90                                                    !!
!!                                                                           !!
!!  Author: Oliver Conway                             Start date: 14/01/2021 !!
!!                                                                           !!
!!  Purpose: Class to store tridiagonal matrix as three vectors:a, b and c.  !!
!!    The class will calculate the a b and c diagonals and contains a        !!
!!    tridiagonal matrix equation solver solver using the Thomas algorithm,  !!
!!    which solves equations of form Ax=B, where A and B are known.          !!
!!                                                                           !!
!!  Revisions:                                                               !!
!!    17/02/2021: Added power iteration                                      !!
!!    28/04/2021: Moved to part of matrix hierarchy                          !!
!!                                                                           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IMPLICIT NONE
! Type definition
TYPE,PUBLIC,extends(matrix) :: tridiagonal_matrix ! This will be the name we instantiate
! Instance variables.
real(dp), allocatable, dimension(:) :: a ! Bottom (leftmost) diagonal in matrix |b1 c1 0  0 |
real(dp), allocatable, dimension(:) :: b ! Middle (main) diagonal in matrix     |a2 b2 c2 0 |
real(dp), allocatable, dimension(:) :: c ! Top (rightmost) diagonal             |0  a3 b3 c3|
                                         !                                      |0  0  a4 b4|
CONTAINS
! Bound procedures
PROCEDURE,PUBLIC :: set_variables => set_variables_sub ! Allows user to input desired matrix
PROCEDURE,PUBLIC :: solve => solve_fn ! Performs a tridiagonal solve using Thomas Algorithm
PROCEDURE,PUBLIC :: get_a => get_a_fn ! Returns array a
PROCEDURE,PUBLIC :: get_b => get_b_fn ! Returns array b
PROCEDURE,PUBLIC :: get_c => get_c_fn ! Returns array c
END TYPE tridiagonal_matrix
! Restrict access to the actual procedure names
PRIVATE :: set_variables_sub
private :: solve_fn
private :: get_a_fn
private :: get_b_fn
private :: get_c_fn
! Now add methods
CONTAINS

SUBROUTINE set_variables_sub(this, a, b, c)
  !
  ! Subroutine to set the variables
  !
  IMPLICIT NONE
  ! Declare calling arguments
  CLASS(matrix) :: this ! Matrix object
  real(dp),INTENT(IN),dimension(:) :: a
  real(dp),INTENT(IN),dimension(:) :: b
  real(dp),INTENT(IN),dimension(:) :: c
  ! Save data
  this%a = a
  this%b = b
  this%c = c
END SUBROUTINE set_variables_sub

FUNCTION solve_fn(this, source_flux) result(solution)
  !
  ! Function to return solution to Ax=B matrix equation, where A is tridiagonal
  !
  IMPLICIT NONE
  ! Declare calling arguments
  CLASS(matrix),INTENT(IN) :: this ! Matrix object
  real(dp), INTENT(IN), dimension(:) :: source_flux
  REAL(dp), DIMENSION( SIZE(this%a) ) :: solution
  REAL(dp), DIMENSION( SIZE(this%b) ) :: btemp
  REAL(dp), DIMENSION( SIZE(source_flux) ) :: dtemp
  REAL(dp) :: w
  INTEGER :: j ! do loop
  INTEGER :: k ! second do loop
  solution = 0.0 ! clear whole vector
  ! Evaluate expression.
  btemp=this%b
  dtemp=source_flux
  DO j = 2,SIZE(this%b)
    w=this%a(j)/btemp(j-1)
    btemp(j)=btemp(j)-w*this%c(j-1)
    dtemp(j)=dtemp(j)-w*dtemp(j-1)
  END DO
  DO k = SIZE(this%b),1,-1
    IF (k == SIZE(this%a)) THEN ! different for last term
      solution(k)=dtemp(k)/btemp(k)
    ELSE ! calculation for other terms
      solution(k)=(dtemp(k)-(this%c(k)*solution(k+1)))/btemp(k)
    END IF
  END DO
END FUNCTION solve_fn

function get_a_fn(this) result(get_a)
  !
  ! Function to return a
  !
  implicit none
  ! Declare calling arguments
  class(matrix),intent(in) :: this ! Matrix object
  real(dp), allocatable, dimension(:) :: get_a
  get_a = this%a
end function get_a_fn

function get_b_fn(this) result(get_b)
  !
  ! Function to return b
  !
  implicit none
  ! Declare calling arguments
  class(matrix),intent(in) :: this ! Matrix object
  real(dp), allocatable, dimension(:) :: get_b
  get_b = this%b
end function get_b_fn

function get_c_fn(this) result(get_c)
  !
  ! Function to return c
  !
  implicit none
  ! Declare calling arguments
  class(matrix),intent(in) :: this ! Matrix object
  real(dp), allocatable, dimension(:) :: get_c
  get_c = this%c
end function get_c_fn

END MODULE tridiagonal_matrix_class
