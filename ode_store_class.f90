MODULE ode_store_class
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!  Filename: ode_store_class.f90                                            !!
  !!                                                                           !!
  !!  Dependant files:                                                         !!
  !!                                                                           !!
  !!  Author: Oliver Conway                             Start date: 20/06/2021 !!
  !!                                                                           !!
  !!  Purpose: Class to store ODEs                                             !!
  !!                                                                           !!
  !!  Revisions:                                                               !!
  !!                                                                           !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  IMPLICIT NONE
  ! Type definition
  TYPE,PUBLIC :: odes ! This will be the name we instantiate
  ! Instance variables.
  CONTAINS
  ! Bound procedures
  PROCEDURE,PUBLIC :: lorenz => lorenz_sub ! Lorenz system of equations
  
  END TYPE odes
  ! Restrict access to the actual procedure names
  PRIVATE :: lorenz_sub
  ! Now add methods
  CONTAINS

  SUBROUTINE lorenz_sub(this,time,values,differentials)
    !
    ! Subroutine containing the lorenz system of equations
    !
    IMPLICIT NONE
    ! Declare calling arguments
    CLASS(odes) :: this ! ODEs object
    real(kind=8),INTENT(INOUT),dimension(3) :: differentials,values
    real(kind=8),INTENT(IN) :: time
    real(kind=8) :: sigma,rho,beta
    ! Set constants
    sigma=1._8
    rho=1._8
    beta=1._8
    differentials(1)=sigma*(values(2)-values(1))
    differentials(2)=values(1)*(rho-values(3))-values(2)
    differentials(3)=(values(1)*values(2))-(beta*values(3))
  END SUBROUTINE lorenz_sub

  END MODULE ode_store_class