SUBROUTINE lorenz(time,values,differentials)
    !
    ! Subroutine containing the lorenz system of equations
    !
    IMPLICIT NONE
    ! Declare calling arguments
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
  END SUBROUTINE lorenz
program test_ode
    use ShampGordon_wrapper
    ! use ode_store_class
    use ShampGordon_oli
    implicit none
    ! procedure(ode%lorenz),POINTER :: f => null()
    external lorenz
    real(kind=8),DIMENSION(3) :: values
    real(kind=8) :: t,tout
    ! f=>ode%lorenz
    values=1.0_8
    t=0.0_8
    tout=1.0_8
    call wrapper(lorenz,size(values),values,t,tout)
    print*,'new values'
    print*,'x',values(1),'y',values(2),'z',values(3)
end program test_ode