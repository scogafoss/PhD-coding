MODULE nuclear_matrix_class
  use precision_set
  use line_class
  use material_class
  use matrix_class
  use region_class_1d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Filename: nuclear_matrix_class.f90                                       !!
!!                                                                           !!
!!  Dependant files: precision_set.f90, line_class.f90, material_class.f90   !!
!!    matrix_class.f90, region_class_1d.f90                                  !!
!!                                                                           !!
!!  Author: Oliver Conway                             Start date: 18/01/2021 !!
!!                                                                           !!
!!  Purpose: Class to store tridiagonal matrix as three vectors:a, b and c.  !!
!!    a, b and c can be read in using line and material classes. The class   !!
!!    calculates discretisation based on input values.                       !!
!!                                                                           !!
!!  Revisions:                                                               !!
!!    15/02/2021: Updated to allow use of array of region class objects      !!
!!    09/03/2021: Changed all instances of get_absorption() to get_removal() !!
!!    10/03/2021: Added index to check which group is required and added     !!
!!      scatter cross section to the diffusion coefficient                   !!
!!                                                                           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IMPLICIT NONE
! Type definition
TYPE,PUBLIC :: nuclear_matrix ! This will be the name we instantiate
! Instance variables.
PRIVATE
real(dp), allocatable, dimension(:) :: a ! Bottom (leftmost) diagonal in matrix |b1 c1 0  0 |
real(dp), allocatable, dimension(:) :: b ! Middle (main) diagonal in matrix     |a2 b2 c2 0 |
real(dp), allocatable, dimension(:) :: c ! Top (rightmost) diagonal             |0  a3 b3 c3|
                                         !                                      |0  0  a4 b4|
CONTAINS
! Bound procedures
PROCEDURE,PUBLIC :: set_variables => set_variables_sub ! Allows user to input desired matrix
procedure,public :: discretise => discretise_sub ! Calculates the variables based on input values
procedure,public :: discretise_regions => discretise_regions_sub ! Discretises input region array.
PROCEDURE,PUBLIC :: get_a => get_a_fn ! Returns array a
PROCEDURE,PUBLIC :: get_b => get_b_fn ! Returns array b
PROCEDURE,PUBLIC :: get_c => get_c_fn ! Returns array c
END TYPE nuclear_matrix
! Restrict access to the actual procedure names
PRIVATE :: set_variables_sub
private :: discretise_sub
private :: discretise_regions_sub
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
  CLASS(nuclear_matrix) :: this ! Matrix object
  real(dp),INTENT(IN),allocatable,dimension(:) :: a
  real(dp),INTENT(IN),allocatable,dimension(:) :: b
  real(dp),INTENT(IN),allocatable,dimension(:) :: c
  ! Save data
  this%a = a
  this%b = b
  this%c = c
END SUBROUTINE set_variables_sub
SUBROUTINE discretise_sub(this, lines, materials)
  !
  ! Subroutine to perform discretisation.
  !
  implicit none
  class(nuclear_matrix), intent(out) :: this ! Nuclear_matrix object
  class(line), intent(in) :: lines ! Line object
  class(material), intent(in) :: materials ! Material object
  INTEGER :: i ! integer for do loop
  real(dp) :: i_double ! double equal to i for do loop
  real(dp) :: delta ! Node spacing
  real(dp) :: D ! Diffusion coefficient
  ! Allocate the class array sizes.
  allocate (this%a(1:lines%get_steps()+1))
  allocate (this%b(1:size(this%a)))
  allocate (this%c(1:size(this%a)))
  delta = lines%get_length()/lines%get_steps()
  D = 1 / (3 * materials%get_removal(1))
  DO i = 1 , SIZE(this%b)
    i_double=i
  !------------------------------------------------------------------------------
    IF (i == 1 .AND. materials%get_left_boundary() == 'z') THEN ! Zero flux left boundary
      this%a(1) = 0 ! set to zero, even though technically not present
      this%b(1) = 1.0E30 ! very large number
      this%c(1) = -(D) / (delta**2) ! a0,1
    ELSE IF ( i == 1 .and. materials%get_left_boundary() == 'v') THEN ! Vacuum left boundary
      this%a(1) = 0 ! set to zero even though not present
      ! this%b(1) = materials%get_absorption()+(2*D/(delta**2))+(D/delta) ! a0,0 - (delta/D) a0,-1
      this%b(1) = materials%get_removal(1)+(2*D/(delta**2))+(1/delta) ! a0,0 - (delta/D) a0,-1
      this%c(1) = -2*D/(delta**2) !a0,1 + a0,-1
    else if (i == 1 .and. materials%get_left_boundary() == 'r') then ! Reflective left boundary
      this%a(1) = 0 ! set to zero, even though technically not present
      this%b(1) = (materials%get_removal(1)+2*(D)/(delta**2)) ! a0,0
      this%c(1) = -2 * D / (delta**2) ! a0,1 + a0,-1
    else if (i ==1 .and. materials%get_left_boundary() == 's') then ! Surface left boundary
      print *, 'Have you adjusted source flux?'
      this%a(1) = 0
      this%b(1) = materials%get_removal(1)+(2*D)/(delta**2)+(1/delta) !  a0,0 - (delta/D) a0,-1 note: -1/delta = delta/D * -0.5*(D+D)/delta^2
      this%c(1) = ((-2*D)/(delta**2)) !a0,1 + a0,-1
    else if (i ==1 .and. materials%get_left_boundary() == 'a') then ! Albedo left boundary
      this%a(1) = 0
      this%b(1) = materials%get_removal(1)+((2*D)/(delta**2)) - &
      ((1/delta)*(materials%get_left_albedo()-1)/(materials%get_left_albedo()+1)) ! aii + aii-1(delta/D)(alpha-1)/(alpha+1) note: -1/delta = delta/D * -0.5*(D+D)/delta^2
      this%c(1) = ((-2*D)/(delta**2)) ! aii+1 + aii-1
  ELSE IF (i == SIZE(this%b) .and. materials%get_right_boundary() == 'z') THEN ! Zero flux right boundary
      this%c(i) = 0 ! set to zero, even though technically not present
      this%b(i) = 1.0E30 ! very large number?
      this%a(i) = (-D / (delta**2)) *((1-(1/(2*(i_double-1))))**lines%get_geomtype())
    ELSE IF (i == SIZE(this%b) .and. materials%get_right_boundary() == 'v') THEN ! Vacuum right boundary
      this%a(i) = ((-D/(delta**2))*(1-(1/(2*(i_double-1))))**lines%get_geomtype())+((-D/(delta**2))*(1+(1/(2*(i_double-1))))&
      **lines%get_geomtype()) ! aII-1 + aII+1
      ! this%b(i) = materials%get_absorption()+((D/(delta**2))*(1-(1/(2*(i_double-1))))**lines%get_geomtype()) + &
      ! ((D/(delta**2))*(1+(1/(2*(i_double-1))))**lines%get_geomtype())-(1/D)*((1-(1/(2*(i_double-1))))**lines%get_geomtype()) ! aII - delta/D aII+1
      this%b(i) = materials%get_removal(1)+((D/(delta**2))*(1-(1/(2*(i_double-1))))**lines%get_geomtype()) + &
      ((D/(delta**2))*(1+(1/(2*(i_double-1))))**lines%get_geomtype())+(1/delta)*((1+(1/(2*(i_double-1))))**lines%get_geomtype()) ! aII - delta/D aII+1
      this%c(i) = 0 ! set to zero, even though not present
    else if (i == SIZE(this%b) .and. materials%get_right_boundary() == 'r') then ! Reflective right boundary
      this%a(i) = -2 * D / (delta**2) ! a0,1 + a0,-1
      this%b(i) = materials%get_removal(1)+(((D)/(delta**2))*(1-(1/(2*(i_double-1))))**lines%get_geomtype()) + &
       ((D/(delta**2))*(1+(1/(2*(i_double-1))))**lines%get_geomtype())! a0,0
      this%c(i) = 0 ! set to zero, even though technically not present
    else if (i == SIZE(this%b) .and. materials%get_right_boundary() == 's') then ! Surface right boundary
      print *, 'Have you adjusted source flux?'
      this%a(i) = ((-D/(delta**2))*(1-(1/(2*(i_double-1))))**lines%get_geomtype())+((-D/(delta**2))*(1+(1/(2*(i_double-1))))&
      **lines%get_geomtype()) ! aII-1 + aII+1
      this%b(i) = materials%get_removal(1)+((D/(delta**2))*(1-(1/(2*(i_double-1))))**lines%get_geomtype()) + &
      ((D/(delta**2))*(1+(1/(2*(i_double-1))))**lines%get_geomtype())-(1/D)*((1-(1/(2*(i_double-1))))**lines%get_geomtype()) ! aII - delta/D aII+1
      this%c(i) = 0 ! set to zero, even though not present
    else if (i == SIZE(this%b) .and. materials%get_right_boundary() == 'a') then ! Albedo right boundary
      this%a(i) = ((-D/(delta**2))*(1-(1/(2*(i_double-1))))**lines%get_geomtype())+((-D/(delta**2))*(1+(1/(2*(i_double-1))))&
      **lines%get_geomtype()) ! aII-1 + aII+1
      this%b(i) = materials%get_removal(1)+((D/(delta**2))*(1-(1/(2*(i_double-1))))**lines%get_geomtype()) + &
      ((D/(delta**2))*(1+(1/(2*(i_double-1))))**lines%get_geomtype()) - &
      ((1/delta)*((1+(1/(2*(i_double-1))))**lines%get_geomtype())*(materials%get_right_albedo()-1)/(materials%get_right_albedo()+1)) ! aII - (delta/D)(alpha-1)/(alpha+1) aII+1

      this%c(i) = 0 ! set to zero, even though not present
  !------------------------------------------------------------------------------
    ELSE
      this%a(i) = (-D / (delta**2))*((1-(1/(2*(i_double-1))))**lines%get_geomtype())
      this%b(i) = materials%get_removal(1) + (D/(delta**2))*((1-(1/(2*(i_double-1))))**lines%get_geomtype()) + &
      ((D/(delta**2))*(1+(1/(2*(i_double-1))))**lines%get_geomtype())
      this%c(i) = (-D / (delta**2))*((1+(1/(2*(i_double-1))))**lines%get_geomtype())
    END IF
  END DO
  !print *, tridiagonal
  !------------------------------------------------------------------------------
  ! Check if this is correct
  !tridiagonal(1,:) = a
  !tridiagonal(2,:) = b
  !tridiagonal(3,:) = c
END SUBROUTINE discretise_sub
SUBROUTINE discretise_regions_sub(this, regions,group)
  !
  ! Subroutine to perform discretisation.
  !
  implicit none
  class(nuclear_matrix), intent(out) :: this ! Nuclear_matrix object
  type(region_1d), intent(in), dimension(:) :: regions ! Region object
  integer,dimension(size(regions)) :: boundary_tracker ! Labels the values of i where boundaries between regions are
  INTEGER :: i ! integer for do loop
  real(dp),allocatable,dimension(:) :: x_coordinate ! Tracks x position in the system
  integer :: region_iterator ! Iterates over the array of regions.
  real(dp) :: i_double ! double equal to i for do loop
  real(dp) :: delta ! Node spacing
  real(dp) :: D ! Diffusion coefficient
  real(dp) :: absorption ! Macroscopic absorption coefficient.
  integer :: total_steps ! Total steps across regions.
  INTEGER :: group ! Which group's data is required
  ! Allocate the class array sizes.
  total_steps = 0
  do region_iterator =1,size(regions)
    total_steps = total_steps + regions(region_iterator)%get_steps()
    boundary_tracker(region_iterator) = total_steps+1 ! tracks the x values where there is a boundary, also the last boundary
  end do
  allocate(x_coordinate(1:total_steps+1))
  region_iterator = 1
  ! Track the x coordinate for each region
  do i = 1, total_steps+1
    if (i==1) then
      x_coordinate(i) = regions(region_iterator)%get_start()
    ! At a boundnary
    else if (i == boundary_tracker(region_iterator) .and. i /= size(this%b)) then
      x_coordinate(i) = x_coordinate(i-1) + regions(region_iterator)%get_length()/regions(region_iterator)%get_steps()
      region_iterator = region_iterator + 1
    else
      x_coordinate(i) = x_coordinate(i-1) + regions(region_iterator)%get_length()/regions(region_iterator)%get_steps()
    end if
  end do
  ! print *, 'x',x_coordinate
  allocate (this%a(1:total_steps+1))
  allocate (this%b(1:size(this%a)))
  allocate (this%c(1:size(this%a)))
  delta = regions(1)%get_length()/regions(1)%get_steps()
  absorption = regions(1)%get_removal(group)
  D = 1 / (3 * (regions(1)%get_absorption()(group)+sum(regions(1)%get_scatter(group,:))))
  region_iterator = 1
  DO i = 1 , SIZE(this%b)
    i_double=i
    ! Check if at a boundary, where average values are needed. No need for last value of i.
    if (i == boundary_tracker(region_iterator) .and. i /= size(this%b)) then
      ! D = ((1/(3*absorption))*delta+(1/(3*(regions(region_iterator+1)%get_absorption()))*(regions(region_iterator+1)%get_length()&
      ! /(regions(region_iterator+1)%get_steps()))))/(delta+(regions(region_iterator+1)%get_length()&
      ! /(regions(region_iterator+1)%get_steps())))
      absorption = ((absorption*delta)+(regions(region_iterator+1)%get_removal(group)*regions(region_iterator+1)%get_length()&
      /(regions(region_iterator+1)%get_steps())))/(delta+(regions(region_iterator+1)%get_length()&
      /(regions(region_iterator+1)%get_steps())))
      delta = (delta + (regions(region_iterator+1)%get_length()/(regions(region_iterator+1)%get_steps())))/2
      ! region_iterator=region_iterator+1
    end if
    ! Check if just after a boudary, where values are updated to new region
    if (i==boundary_tracker(region_iterator)+1 .and. i /= size(this%b)) then
      region_iterator=region_iterator+1
      delta = regions(region_iterator)%get_length()/regions(region_iterator)%get_steps()
      absorption = regions(region_iterator)%get_removal(group)
      D = 1 / (3 * (regions(region_iterator)%get_absorption()(group)+sum(regions(region_iterator)%get_scatter(group,:))))
    end if
  !
  !------------------------------------------------------------------------------
  !
    ! At the left boundary, only basic coefficients needed.
    IF (i == 1 .AND. regions(region_iterator)%get_left_boundary() == 'z') THEN ! Zero flux left boundary
      this%a(1) = 0 ! set to zero, even though technically not present
      this%b(1) = 1.0E30 ! very large number
      this%c(1) = -(D) / (delta**2)
    ELSE IF ( i == 1 .and. regions(region_iterator)%get_left_boundary() == 'v') THEN ! Vacuum left boundary
      this%a(1) = 0 ! set to zero even though not present
      ! this%b(1) = absorption+(2*D/(delta**2))+(D/delta)
      this%b(1) = absorption+(2*D/(delta**2))+(1/delta)
      this%c(1) = -2*D/(delta**2)
    else if (i == 1 .and. regions(region_iterator)%get_left_boundary() == 'r') then ! Reflective left boundary
      this%a(1) = 0 ! set to zero, even though technically not present
      this%b(1) = (absorption+2*(D)/(delta**2)) ! a0,0
      this%c(1) = -2 * D / (delta**2) ! a0,1 + a0,-1
    else if (i ==1 .and. regions(region_iterator)%get_left_boundary() == 's') then ! Surface left boundary
      print *, 'Have you adjusted source flux?'
      this%a(1) = 0
      this%b(1) = absorption+(2*D)/(delta**2)+(1/delta) ! note: -1/delta = delta/D * -0.5*(D+D)/delta^2
      this%c(1) = ((-2*D)/(delta**2))
    else if (i ==1 .and. regions(region_iterator)%get_left_boundary() == 'a') then ! Albedo left boundary
      this%a(1) = 0
      this%b(1) = absorption+((2*D)/(delta**2)) - &
      ((1/delta)*(regions(region_iterator)%get_left_albedo()-1)/(regions(region_iterator)%get_left_albedo()+1)) ! note: -1/delta = delta/D * -0.5*(D+D)/delta^2
      this%c(1) = ((-2*D)/(delta**2))
  !
  !------------------------------------------------------------------------------
  !
    ! At the right boundary, need to correct for geometry and position.
    ELSE IF (i == SIZE(this%b) .and. regions(region_iterator)%get_right_boundary() == 'z') THEN ! Zero flux right boundary
      this%c(i) = 0 ! set to zero, even though technically not present
      this%b(i) = 1.0E30 ! very large number?
      ! aii-1
      this%a(i) = -(D/(delta**2))*(1-(delta/(2*x_coordinate(i))))**regions(region_iterator)%get_geomtype()
    ELSE IF (i == SIZE(this%b) .and. regions(region_iterator)%get_right_boundary() == 'v') THEN ! Vacuum right boundary
      ! aii-1 + aii+1
      this%a(i) = (-(D/(delta**2))*(1-(delta/(2*x_coordinate(i))))**regions(region_iterator)%get_geomtype())+(-(D/(delta**2))*(1+(delta/(2*x_coordinate(i))))**regions(region_iterator)%get_geomtype())
      ! this%b(i) = absorption+((D/(delta**2))*(1-(1/(2*(i_double-1))))**regions(region_iterator)%get_geomtype()) + &
      ! ((D/(delta**2))*(1+(1/(2*(i_double-1))))**regions(region_iterator)%get_geomtype())-(1/D)*((1-(1/(2*(i_double-1))))**regions(region_iterator)%get_geomtype()) ! aII - delta/D aII+1
      ! aii - aii+1(delta/D)
      this%b(i) = (absorption + (D/(delta**2))*(((1-(delta/(2*x_coordinate(i))))**regions(region_iterator)%get_geomtype())+((1+(delta/(2*x_coordinate(i))))&
      **regions(region_iterator)%get_geomtype())))-((-(D/(delta**2))*(1+(delta/(2*x_coordinate(i))))**regions(region_iterator)%get_geomtype())*(delta/D))
      this%c(i) = 0 ! set to zero, even though not present
    else if (i == SIZE(this%b) .and. regions(region_iterator)%get_right_boundary() == 'r') then ! Reflective right boundary
      !aii-1+aii+1
      this%a(i) = (-(D/(delta**2))*(1-(delta/(2*x_coordinate(i))))**regions(region_iterator)%get_geomtype())+(-(D/(delta**2))*(1+(delta/(2*x_coordinate(i))))**regions(region_iterator)%get_geomtype())
      !aii
      this%b(i) = absorption + (D/(delta**2))*(((1-(delta/(2*x_coordinate(i))))**regions(region_iterator)%get_geomtype())+((1+(delta/(2*x_coordinate(i))))&
      **regions(region_iterator)%get_geomtype()))
      this%c(i) = 0 ! set to zero, even though technically not present
    else if (i == SIZE(this%b) .and. regions(region_iterator)%get_right_boundary() == 's') then ! Surface right boundary
      print *, 'Have you adjusted source flux?'
      ! aii-1 + aii+1
      this%a(i) = (-(D/(delta**2))*(1-(delta/(2*x_coordinate(i))))**regions(region_iterator)%get_geomtype())+(-(D/(delta**2))*(1+(delta/(2*x_coordinate(i))))**regions(region_iterator)%get_geomtype())
      ! this%b(i) = absorption+((D/(delta**2))*(1-(1/(2*(i_double-1))))**regions(region_iterator)%get_geomtype()) + &
      ! ((D/(delta**2))*(1+(1/(2*(i_double-1))))**regions(region_iterator)%get_geomtype())-(1/D)*((1-(1/(2*(i_double-1))))**regions(region_iterator)%get_geomtype()) ! aII - delta/D aII+1
      ! aii - aii+1(delta/D)
      this%b(i) = (absorption + (D/(delta**2))*(((1-(delta/(2*x_coordinate(i))))**regions(region_iterator)%get_geomtype())+((1+(delta/(2*x_coordinate(i))))&
      **regions(region_iterator)%get_geomtype())))-((-(D/(delta**2))*(1+(delta/(2*x_coordinate(i))))**regions(region_iterator)%get_geomtype())*(delta/D))
      this%c(i) = 0 ! set to zero, even though not present
    else if (i == SIZE(this%b) .and. regions(region_iterator)%get_right_boundary() == 'a') then ! Albedo right boundary
      !aii-1 + aii+1
      this%a(i) = (-(D/(delta**2))*(1-(delta/(2*x_coordinate(i))))**regions(region_iterator)%get_geomtype())+(-(D/(delta**2))*(1+(delta/(2*x_coordinate(i))))**regions(region_iterator)%get_geomtype())
      !aii + aii+1(delta/D)(alpha-1)/(alpha+1)
      this%b(i) = (absorption + (D/(delta**2))*(((1-(delta/(2*x_coordinate(i))))**regions(region_iterator)%get_geomtype())+((1+(delta/(2*x_coordinate(i))))&
      **regions(region_iterator)%get_geomtype())))+((-(D/(delta**2))*(1+(delta/(2*x_coordinate(i))))**regions(region_iterator)%get_geomtype())*(delta/D)*((regions(region_iterator)%get_right_albedo()-1)&
      /(regions(region_iterator)%get_right_albedo()+1)))
      this%c(i) = 0 ! set to zero, even though not present
  !
  !------------------------------------------------------------------------------
  !
    ! At each boundary need more terms
  else if(i == boundary_tracker(region_iterator) .and. i /= size(this%b)) then
      ! ai,i-1
      this%a(i) = -(((1/(3*(regions(region_iterator)%get_absorption()(group)+sum(regions(region_iterator)%get_scatter()(group,:)))))/(delta*regions(region_iterator)%get_length()/&
      regions(region_iterator)%get_steps()))*(1-((regions(region_iterator)%get_length()/regions(region_iterator)%get_steps())/(2*x_coordinate(i))))**regions(region_iterator)%get_geomtype())
      ! ai,i
      this%b(i) = absorption + (((1/(3*(regions(region_iterator)%get_absorption()(group)+sum(regions(region_iterator)%get_scatter()(group,:)))))/(delta*regions(region_iterator)%get_length()/&
      regions(region_iterator)%get_steps()))*(1-((regions(region_iterator)%get_length()/regions(region_iterator)%get_steps())/(2*x_coordinate(i))))**regions(region_iterator)%get_geomtype()) +&
       (((1/(3*(regions(region_iterator+1)%get_absorption()(group)+sum(regions(region_iterator+1)%get_scatter()(group,:)))))/(delta*regions(region_iterator+1)%get_length()/&
       regions(region_iterator+1)%get_steps()))*(1+((regions(region_iterator+1)%get_length()/regions(region_iterator+1)%get_steps())/(2*x_coordinate(i))))**regions(region_iterator+1)%get_geomtype())
      ! ai,i+1
      this%c(i) = -(((1/(3*(regions(region_iterator+1)%get_absorption()(group)+sum(regions(region_iterator+1)%get_scatter()(group,:)))))/&
      (delta*regions(region_iterator+1)%get_length()/regions(region_iterator+1)%get_steps()))*(1+((regions(region_iterator+1)%get_length()/&
      regions(region_iterator+1)%get_steps())/(2*x_coordinate(i))))**regions(region_iterator+1)%get_geomtype())
  !
  !------------------------------------------------------------------------------
  !
    ! Elsewhere need to correct for geometry and position
    ELSE
      ! ai,i-1
      this%a(i) = -(D/(delta**2))*(1-(delta/(2*x_coordinate(i))))**regions(region_iterator)%get_geomtype()
      ! ai,i
      this%b(i) = absorption + (D/(delta**2))*(((1-(delta/(2*x_coordinate(i))))**regions(region_iterator)%get_geomtype())+((1+(delta/(2*x_coordinate(i))))&
      **regions(region_iterator)%get_geomtype()))
      ! ai,i+1
      this%c(i) = -(D/(delta**2))*(1+(delta/(2*x_coordinate(i))))**regions(region_iterator)%get_geomtype()
    END IF
  END DO
  !print *, tridiagonal
  !------------------------------------------------------------------------------
  ! Check if this is correct
  !tridiagonal(1,:) = a
  !tridiagonal(2,:) = b
  !tridiagonal(3,:) = c
END SUBROUTINE discretise_regions_sub
function get_a_fn(this) result(get_a)
  !
  ! Function to return a
  !
  implicit none
  ! Declare calling arguments
  class(nuclear_matrix),intent(in) :: this ! Matrix object
  real(dp), allocatable, dimension(:) :: get_a
  get_a = this%a
end function get_a_fn
function get_b_fn(this) result(get_b)
  !
  ! Function to return b
  !
  implicit none
  ! Declare calling arguments
  class(nuclear_matrix),intent(in) :: this ! Matrix object
  real(dp), allocatable, dimension(:) :: get_b
  get_b = this%b
end function get_b_fn
function get_c_fn(this) result(get_c)
  !
  ! Function to return c
  !
  implicit none
  ! Declare calling arguments
  class(nuclear_matrix),intent(in) :: this ! Matrix object
  real(dp), allocatable, dimension(:) :: get_c
  get_c = this%c
end function get_c_fn
END MODULE nuclear_matrix_class
