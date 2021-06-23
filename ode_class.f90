MODULE ode_class
    use precision_set
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!  Filename: ode_class.f90                                                  !!
  !!                                                                           !!
  !!  Dependant files: precision_set.f90                                       !!
  !!                                                                           !!
  !!  Author: Oliver Conway                             Start date: 23/06/2021 !!
  !!                                                                           !!
  !!  Purpose: Class to perfrom ODE solves in time                             !!
  !!                                                                           !!
  !!  Revisions:                                                               !!
  !!                                                                           !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  IMPLICIT NONE
  ! Type definition
  TYPE,PUBLIC :: ode ! This will be the name we instantiate
  ! Instance variables.
  PRIVATE
  real(dp),dimension(:),allocatable :: solutions
  CONTAINS
  ! Bound procedures
  PROCEDURE,PUBLIC :: explicit_euler => explicit_euler_sub ! Performs the explicit Euler method
  PROCEDURE,PUBLIC :: implicit_euler => implicit_euler_sub
  PROCEDURE,PUBLIC :: get_solution => get_solution_fn ! Returns the calculated solution
  PROCEDURE,PUBLIC :: runge_kutta => runge_kutta_sub
  END TYPE ode
  ! Restrict access to the actual procedure names
  PRIVATE :: explicit_euler_sub
  PRIVATE :: implicit_euler_sub
  PRIVATE :: get_solution_fn
  PRIVATE :: runge_kutta_sub
  ! Now add methods
  CONTAINS

  SUBROUTINE explicit_euler_sub(this,derivative,initial,delta,steps)
    !
    ! Subroutine to set the variables
    !
    IMPLICIT NONE
    ! Declare calling arguments
    CLASS(ode),INTENT(OUT) :: this ! ODE object
    CHARACTER(len=*),INTENT(IN) :: derivative
    real(dp),INTENT(IN) :: initial,delta
    INTEGER,INTENT(IN) :: steps ! Number of solutions (time steps) required
    INTEGER :: i
    real(dp) :: f,initial_new ! Derivative funciton
    !
    ! Allocate the number of solutions
    !
    allocate(this%solutions(1:steps))
    !
    ! Fill the solutions array
    !
    initial_new=initial
    do i=1,steps
        this%solutions(i) = initial_new
        if(derivative=='radiation') then
            f=radiation_derivative(initial_new)
        endif
        initial_new=initial_new + (delta*f)        
    end do
  END SUBROUTINE explicit_euler_sub

  SUBROUTINE implicit_euler_sub(this,derivative,initial,delta,steps)
    !
    ! Subroutine to set the variables
    !
    IMPLICIT NONE
    ! Declare calling arguments
    CLASS(ode),INTENT(OUT) :: this ! ODE object
    CHARACTER(len=*),INTENT(IN) :: derivative
    real(dp),INTENT(IN) :: initial,delta
    INTEGER,INTENT(IN) :: steps ! Number of solutions (time steps) required
    INTEGER :: i
    real(dp) :: f,initial_new ! Derivative funciton
    !
    ! Allocate the number of solutions
    !
    allocate(this%solutions(1:steps))
    !
    ! Fill the solutions array
    !
    ! Fill the first immediately
    this%solutions(1) = initial
    ! Set the first solution outside of loop
    initial_new=newton_raphson(derivative,initial,delta)
    do i=2,steps
        ! Add each new solution to the solutions array
        this%solutions(i) = initial_new
        ! Perform newton raphson to calculate the next value
        initial_new=newton_raphson(derivative,initial_new,delta)        
    end do
  END SUBROUTINE implicit_euler_sub

  function get_solution_fn(this) result(solution)
    !
    ! Subroutine to set the variables
    !
    IMPLICIT NONE
    ! Declare calling arguments
    CLASS(ode),INTENT(IN) :: this ! ODE object
    real(dp),DIMENSION(:),ALLOCATABLE :: solution ! The solution found by this method
    !
    ! Check if allocated
    !
    if(.not.allocated(this%solutions)) stop 'Error, ODE solutions not yet calculated.'
    ALLOCATE(solution(size(this%solutions)))
    !
    ! Return the solutions
    !
    solution = this%solutions
  END function get_solution_fn

  real(dp) function newton_raphson(derivative,initial,delta)
    !
    ! Function to calculate newton raphson
    !
    IMPLICIT NONE
    ! Declare calling arguments
    CHARACTER(len=*),INTENT(IN) :: derivative
    real(dp),INTENT(IN) :: initial,delta
    INTEGER :: iterations,max
    real(dp) :: f,fprime,convergence,previous ! Derivative funciton
    !
    ! Perform the Newton raphson
    !
    convergence=1e-8_dp
    newton_raphson=initial
    iterations=0
    max=10
    previous=newton_raphson
    ! Loop until convergence
    do
        if(derivative=='radiation') then
            f=newton_raphson-initial+(delta*(-1.0_dp)*radiation_derivative(newton_raphson))
            fprime=1.0_dp+(delta*4.0_dp*4.0e-12_dp*(newton_raphson**3))
        endif
        previous=newton_raphson
        ! print*,'newton =',newton_raphson
        newton_raphson=previous-(f/fprime)
        iterations=iterations+1
        if(abs(newton_raphson-previous)/newton_raphson<=convergence .or. iterations==max) exit
    end do
    if(iterations==max) print*,'Max Newton-Raphson iterations met'
  end function newton_raphson

  SUBROUTINE runge_kutta_sub(this,derivative,initial,delta,steps)
    !
    ! Subroutine to set the variables
    !
    IMPLICIT NONE
    ! Declare calling arguments
    CLASS(ode),INTENT(OUT) :: this ! ODE object
    CHARACTER(len=*),INTENT(IN) :: derivative
    real(dp),INTENT(IN) :: initial,delta
    INTEGER,INTENT(IN) :: steps ! Number of solutions (time steps) required
    INTEGER :: i
    real(dp) :: initial_new,dt1,dt2,dt3,dt4 ! Derivative funciton
    !
    ! Allocate the number of solutions
    !
    allocate(this%solutions(1:steps))
    !
    ! Fill the solutions array
    !
    initial_new=initial
    do i=1,steps
        this%solutions(i) = initial_new
        if(derivative=='radiation') then
            dt1=radiation_derivative(initial_new)*delta
            dt2=radiation_derivative(initial_new+(dt1/2))*delta
            dt3=radiation_derivative(initial_new+(dt2/2))*delta
            dt4=radiation_derivative(initial_new+dt3)*delta
        endif
        initial_new=initial_new + (1.0_dp/6.0_dp)*(dt1+(2*dt2)+(2*dt3)+dt4)        
    end do
  END SUBROUTINE runge_kutta_sub

  real(dp) function radiation_derivative(value)
    !
    ! Derivative function for radiation problem
    !
    IMPLICIT NONE
    ! Declare calling arguments
    real(dp),INTENT(IN) :: value
    radiation_derivative=-4.0e-12_dp*((value**4)-(250.0_dp)**4)    
END function radiation_derivative

END MODULE ode_class
