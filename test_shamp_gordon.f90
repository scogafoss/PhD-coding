program test_ode
    use ShampGordon_wrapper
    use ode_store_class
    use ShampGordon_oli
    implicit none
    type(odes) :: ode
    real(kind=8),DIMENSION(3) :: values
    real(kind=8) :: t,tout
    values=1.0_8
    t=0.0_8
    tout=1.0_8
    call wrapper(ode%lorenz,ode,size(values),values,t,tout)
    print*,'new values'
    print*,'x',values(1),'y',values(2),'z',values(3)
end program test_ode