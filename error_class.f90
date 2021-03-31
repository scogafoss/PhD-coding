MODULE error_class
  use maths_class
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Filename: error_class.f90                                                !!
!!                                                                           !!
!!  Dependant files: precision_set.f90                                       !!
!!                                                                           !!
!!  Author: Oliver Conway                             Start date: 16/12/2020 !!
!!                                                                           !!
!!  Purpose: Class to calculate error for two input arrays using defined     !!
!!    functions. The class will calculate the L2 error for the two arrays to !!
!!    test the second order accuracy of descretised solutions.               !!
!!                                                                           !!
!!  Revisions:                                                               !!
!!                                                                           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IMPLICIT NONE
! Type definition
TYPE,PUBLIC :: error ! This will be the name we instantiate
! Instance variables.
PRIVATE
real(dp), allocatable, dimension(:) :: variable1 ! First variable array
real(dp), allocatable, dimension(:) :: variable2 ! Second variable array
CONTAINS
! Bound procedures
PROCEDURE,PUBLIC :: set_variables => set_variables_sub
PROCEDURE,PUBLIC :: get_l2 => get_l2_fn
procedure,public :: l2_from_fluxes => l2_from_fluxes_fn
END TYPE error
! Restrict access to the actual procedure names
PRIVATE :: set_variables_sub
private :: get_l2_fn
PRIVATE :: l2_from_fluxes_fn
! Now add methods
CONTAINS

SUBROUTINE set_variables_sub(this, variable1, variable2)
  !
  ! Subroutine to set the variables
  !
  IMPLICIT NONE
  ! Declare calling arguments
  CLASS(error) :: this ! Error object
  real(dp),INTENT(IN),allocatable, dimension(:) :: variable1
  real(dp),INTENT(IN),allocatable,dimension(:) :: variable2
  ! Save data
  this%variable1 = variable1
  this%variable2 = variable2
END SUBROUTINE set_variables_sub

real(dp) FUNCTION get_l2_fn(this)
  !
  ! Function to return the l2 error from the two variable arrays
  !
  IMPLICIT NONE
  ! Declare calling arguments
  CLASS(error),INTENT(IN) :: this ! Error object
  real(dp), dimension (size(this%variable1)) :: l2_vector ! Tracks error at each element
  ! Get l2 error
  integer :: l
  if (size(this%variable1)==size(this%variable2)) then
    do l = 1 , size(this%variable1)
      l2_vector(l) = (this%variable1(l)-this%variable2(l))**2
    end do
    get_l2_fn =  (sum(l2_vector)/size(l2_vector))**0.5
  else
    write(*,*) 'Arrays must be of the same dimensions to calculate error.'
  end if
END FUNCTION get_l2_fn

real(dp) function l2_from_fluxes_fn(this,phi,x_phi,phi_gem,x_gem)
  !
  ! Function to return the l2 error from the gem flux and code flux
  !
  IMPLICIT NONE
  ! Declare calling arguments
  CLASS(error),INTENT(INOUT) :: this ! Error object
  real(dp),INTENT(IN),DIMENSION(:) :: phi
  real(dp),INTENT(IN),DIMENSION(:) :: phi_gem
  real(dp),INTENT(IN),DIMENSION(:) :: x_gem
  real(dp),INTENT(IN),DIMENSION(:) :: x_phi
  real(dp),DIMENSION(size(phi_gem)) :: phi_interpolated
  type(maths) :: m1
  call m1%set_variables(x_gem,x_phi,phi)
  phi_interpolated=m1%get_interpolation()
  this%variable1=phi_interpolated
  this%variable2=phi_gem
  l2_from_fluxes_fn = this%get_l2()
end function l2_from_fluxes_fn

END MODULE error_class
