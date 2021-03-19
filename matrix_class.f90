MODULE matrix_class
  use precision_set
  use line_class
  use material_class
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
!!                                                                           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IMPLICIT NONE
! Type definition
TYPE,PUBLIC :: matrix ! This will be the name we instantiate
! Instance variables.
PRIVATE
real(dp), allocatable, dimension(:) :: a ! Bottom (leftmost) diagonal in matrix |b1 c1 0  0 |
real(dp), allocatable, dimension(:) :: b ! Middle (main) diagonal in matrix     |a2 b2 c2 0 |
real(dp), allocatable, dimension(:) :: c ! Top (rightmost) diagonal             |0  a3 b3 c3|
                                         !                                      |0  0  a4 b4|
CONTAINS
! Bound procedures
PROCEDURE,PUBLIC :: set_variables => set_variables_sub ! Allows user to input desired matrix
PROCEDURE,PUBLIC :: thomas_solve => thomas_solve_fn ! Performs a tridiagonal solve
procedure,public :: power_iteration => power_iteration_sub ! Performs the power iteration for given nu*sigma_f (should still work if zero)
PROCEDURE,PUBLIC :: get_a => get_a_fn ! Returns array a
PROCEDURE,PUBLIC :: get_b => get_b_fn ! Returns array b
PROCEDURE,PUBLIC :: get_c => get_c_fn ! Returns array c
END TYPE matrix
! Restrict access to the actual procedure names
PRIVATE :: set_variables_sub
private :: thomas_solve_fn
private :: power_iteration_sub
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
  real(dp),INTENT(IN),allocatable,dimension(:) :: a
  real(dp),INTENT(IN),allocatable,dimension(:) :: b
  real(dp),INTENT(IN),allocatable,dimension(:) :: c
  ! Save data
  this%a = a
  this%b = b
  this%c = c
END SUBROUTINE set_variables_sub
FUNCTION thomas_solve_fn(this, source_flux) result(solution)
  !
  ! Function to return solution to Ax=B matrix equation, whre A is tridiagonal
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
END FUNCTION thomas_solve_fn
SUBROUTINE power_iteration_sub(this, phi, keff, regions)
  !
  ! Subroutine to perfrom power iteration (assumes no volumetric_source)
  !
  IMPLICIT NONE
  ! Declare calling arguments
  CLASS(matrix) :: this ! Matrix object
  integer :: iterations
  real(dp),INTENT(inout),allocatable,dimension(:) :: phi
  real(dp),allocatable,dimension(:) :: phi_temp ! phi from the previous iteration
  real(dp),INTENT(inout) :: keff
  real(dp) :: keff_temp !from the previous iteration
  real(dp),allocatable,dimension(:) :: source
  type(region_1d),INTENT(in),allocatable,dimension(:) :: regions
  integer,dimension(size(regions)) :: boundary_tracker ! Labels the values of i where boundaries between regions are
  real(dp) :: convergence_criterion
  integer :: region_iterator
  integer :: total_steps
  integer :: i
  integer :: j
  real(dp) :: numerator
  real(dp) :: denominator
  real(dp) :: normalisation
  real(dp),allocatable,dimension(:) :: x_coordinate ! Tracks x position in the system
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TEMPORARY FOR JACK
  open(90, file ='x_y_iterations.txt')
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Convergence condition
  convergence_criterion = 1e-5
  ! Record where boundaries are, so can use correct values of fission and delta
  total_steps = 0
  do region_iterator =1,size(regions)
    total_steps = total_steps + regions(region_iterator)%get_steps()
    boundary_tracker(region_iterator) = total_steps+1 ! tracks the x values where there is a boundary, also the last boundary
  end do
  allocate(source(1:total_steps+1))
  ! Initial guesses
  keff = 1
  allocate(phi(1:total_steps+1))
  phi = 1
  allocate(x_coordinate(1:total_steps+1))
  !
  ! Do loop to perform power iteration
  !
  iterations=0
  do
    !
    ! Correct source flux for fission
    !
    region_iterator=1
    do i =1,size(source)
      ! If at a boundary need to average the source flux
      if (i == boundary_tracker(region_iterator) .and. i /= size(source)) then
        source(i) = (1/keff)*((((regions(region_iterator)%get_fission(1)*(regions(region_iterator)%get_length()/&
        regions(region_iterator)%get_steps()))+(regions(region_iterator+1)%get_fission(1)*&
        (regions(region_iterator+1)%get_length()/regions(region_iterator+1)%get_steps())))/&
        ((regions(region_iterator)%get_length()/regions(region_iterator)%get_steps())+&
        (regions(region_iterator+1)%get_length()/regions(region_iterator+1)%get_steps())))*phi(i))
        region_iterator=region_iterator+1
      ! If not at boundary
      else
        source(i)=(1/keff)*((regions(region_iterator)%get_fission(1))*phi(i))
      end if
    end do
    !
    ! Calculate the next phi with thomas algorithm and store the previous
    !
    phi_temp=phi
    phi=this%thomas_solve(source)

    !
    ! Calculate next keff and store the previous
    !
    keff_temp=keff
    numerator = 0
    denominator = 0
    region_iterator=1
    do i = 1,size(phi)
      ! First find numerator and denominator
      ! If at boundary use average values
      if (i == boundary_tracker(region_iterator) .and. i /= size(source)) then
        ! Numerator for source
        numerator = numerator + (((((regions(region_iterator)%get_fission(1)*(regions(region_iterator)%get_length()/&
        regions(region_iterator)%get_steps()))+(regions(region_iterator+1)%get_fission(1)*&
        (regions(region_iterator+1)%get_length()/regions(region_iterator+1)%get_steps())))/&
        ((regions(region_iterator)%get_length()/regions(region_iterator)%get_steps())+&
        (regions(region_iterator+1)%get_length()/regions(region_iterator+1)%get_steps())))*phi(i))*&
        (((regions(region_iterator)%get_length()/regions(region_iterator)%get_steps())+(regions(region_iterator+1)%get_length()/&
        regions(region_iterator+1)%get_steps()))/2))
        ! denominator for source
        denominator = denominator + (((((regions(region_iterator)%get_fission(1)*(regions(region_iterator)%get_length()/&
        regions(region_iterator)%get_steps()))+(regions(region_iterator+1)%get_fission(1)*(regions(region_iterator+1)%get_length()/&
        regions(region_iterator+1)%get_steps())))/((regions(region_iterator)%get_length()/regions(region_iterator)%get_steps())+&
        (regions(region_iterator+1)%get_length()/regions(region_iterator+1)%get_steps())))*phi_temp(i))*&
        (((regions(region_iterator)%get_length()/regions(region_iterator)%get_steps())+(regions(region_iterator+1)%get_length()/&
        regions(region_iterator+1)%get_steps()))/2))
        region_iterator=region_iterator+1
      ! IF not at boundary
      else
        numerator = numerator + (regions(region_iterator)%get_fission(1)*phi(i)*(regions(region_iterator)%get_length()/&
        regions(region_iterator)%get_steps()))
        denominator = denominator + (regions(region_iterator)%get_fission(1)*phi_temp(i)*(regions(region_iterator)%get_length()/&
        regions(region_iterator)%get_steps()))
      end if
    end do
    iterations = iterations +1
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TEMPORARY FOR JACK
    do j=1,size(phi)
      write(90,*) phi(j)
    end do
    write(90,*) 'Iteration',iterations
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TEMPORARY FOR JACK
    ! Calculate the new keff
    keff=keff_temp * (numerator/denominator)
    ! If the convergence_criterion has been met, exit the loop
    ! print *,'test converge k1-k0/k0 = ',abs((keff-keff_temp)/keff_temp)
    if (abs((keff-keff_temp)/keff_temp)<convergence_criterion) exit
  end do
  !
  ! Normalise the flux values
  !
  ! First track the x coordinate for each region
  region_iterator=1
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
  ! Now calculate the normalisation
  do i =1,total_steps ! note ther are total_steps+1 nodes so this goes to second to last node, which is fine as it integrates over steps.
    ! If at boundary
    if (i == boundary_tracker(region_iterator) .and. i /= size(this%b)) then
      region_iterator = region_iterator + 1
      normalisation=normalisation+(0.5*(regions(region_iterator)%get_length()/regions(region_iterator)%get_steps())*&
      ((source(i)*(x_coordinate(i)**regions(region_iterator)%get_geomtype()))+source(i+1)*(x_coordinate(i+1)&
      **regions(region_iterator)%get_geomtype())))
    else
      normalisation=normalisation+(0.5*(regions(region_iterator)%get_length()/regions(region_iterator)%get_steps())*&
      ((source(i)*x_coordinate(i)**regions(region_iterator)%get_geomtype())+source(i+1)*x_coordinate(i+1)&
      **regions(region_iterator)%get_geomtype()))
    end if
  end do
  ! Make correction for geometry
  if(regions(region_iterator)%get_geomtype()==1) then! Cylindrical
    normalisation=2*pi_dp*normalisation
  else if (regions(region_iterator)%get_geomtype()==2) then! spherical
    normalisation=4*pi_dp*normalisation
  end if
  ! Now normalise flux
  phi = phi / normalisation
  print *, 'Number of iterations = ',iterations
  print *,'normalisation = ',normalisation,'sum(s) = ', sum(source)
END SUBROUTINE power_iteration_sub
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
END MODULE matrix_class